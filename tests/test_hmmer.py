import logging

import pytest

from homoppi.config import HmmscanParams
from homoppi.hmmer import filter_domains, parse_domtblout, parse_pfam_dat, resolve_clan_overlaps

logger = logging.getLogger("test")


def domtblout_line(
    target, acc, query, cevalue,
    ali_from=10, ali_to=90, seq_evalue="1.2e-40", tlen=264, hmm_from=5, hmm_to=180,
):
    fields = [
        target, acc, str(tlen), query, "-", "500",
        seq_evalue, "140.0", "0.1", "1", "1",
        str(cevalue), "1.2e-16", "62.0", "0.0",
        str(hmm_from), str(hmm_to), str(ali_from), str(ali_to), "8", "200", "0.90",
        "some description here",
    ]
    return "  ".join(fields)


def write_domtblout(tmp_path, lines):
    path = tmp_path / "raw.domtblout"
    path.write_text("# comment line\n" + "\n".join(lines) + "\n")
    return path


def test_parse_domtblout(tmp_path):
    path = write_domtblout(
        tmp_path,
        [
            domtblout_line("Pkinase", "PF00069.28", "PROT1", "3.4e-20"),
            domtblout_line("CustomDom", "-", "PROT2", "1e-15"),  # no accession -> target name
        ],
    )
    df = parse_domtblout(path)
    assert df.iloc[0].pfam_acc == "PF00069"  # version stripped
    assert df.iloc[0].query_id == "PROT1"
    assert df.iloc[0].ali_from == 10
    assert df.iloc[0].hmm_cov == pytest.approx(100 * (180 - 5 + 1) / 264, abs=0.1)
    assert df.iloc[0].seq_evalue == pytest.approx(1.2e-40)
    assert df.iloc[1].pfam_acc == "CustomDom"


def test_filter_domains_cevalue_and_dedup(tmp_path):
    path = write_domtblout(
        tmp_path,
        [
            domtblout_line("Pkinase", "PF00069.28", "PROT1", "1e-20", ali_from=10),
            domtblout_line("Pkinase", "PF00069.28", "PROT1", "1e-30", ali_from=200),  # better instance
            domtblout_line("SH3_1", "PF00018.28", "PROT1", "1e-05"),  # above cutoff -> dropped
        ],
    )
    df = filter_domains(parse_domtblout(path), HmmscanParams(cevalue=1e-10), 1, logger)
    assert len(df) == 1
    assert df.iloc[0].ali_from == 200  # best c-Evalue instance kept


def test_filter_domains_emulates_hmmscan_full_seq_gate(tmp_path):
    # Re-filtering a cached raw file with a tighter threshold must also apply
    # the full-sequence E-value gate that a fresh `hmmscan -E` run would.
    path = write_domtblout(
        tmp_path,
        [
            domtblout_line("DomA", "PF00001.1", "PROT1", "1e-25", seq_evalue="1e-30"),  # kept
            domtblout_line("DomB", "PF00002.1", "PROT1", "1e-25", seq_evalue="1e-12"),  # seq gate -> dropped
        ],
    )
    df = filter_domains(parse_domtblout(path), HmmscanParams(cevalue=1e-20), 1, logger)
    assert df.pfam_acc.tolist() == ["PF00001"]


def test_filter_domains_cut_tc_skips_cevalue(tmp_path):
    path = write_domtblout(tmp_path, [domtblout_line("SH3_1", "PF00018.28", "PROT1", "1e-05")])
    df = filter_domains(parse_domtblout(path), HmmscanParams(cut_tc=True), 1, logger)
    assert len(df) == 1


def test_min_hmm_coverage_filter(tmp_path):
    path = write_domtblout(
        tmp_path,
        [
            domtblout_line("DomA", "PF00001.1", "PROT1", "1e-20", hmm_from=5, hmm_to=180),  # 66.7%
            domtblout_line("DomB", "PF00002.1", "PROT1", "1e-20", hmm_from=5, hmm_to=20),  # 6.1%
        ],
    )
    df = filter_domains(parse_domtblout(path), HmmscanParams(min_hmm_coverage=50.0), 1, logger)
    assert df.pfam_acc.tolist() == ["PF00001"]


CLANS = {"PF04851": "CL0023", "PF00270": "CL0023", "PF99999": "CL9999"}


def test_clan_overlap_resolution(tmp_path):
    path = write_domtblout(
        tmp_path,
        [
            # ResIII vs DEAD: same clan, overlapping -> only the better (DEAD) survives
            domtblout_line("ResIII", "PF04851.1", "PROT1", "1e-20", ali_from=135, ali_to=300),
            domtblout_line("DEAD", "PF00270.1", "PROT1", "9e-21", ali_from=137, ali_to=303),
            # different clan overlapping DEAD -> kept
            domtblout_line("OtherClan", "PF99999.1", "PROT1", "1e-15", ali_from=140, ali_to=280),
            # family without a clan overlapping DEAD -> kept
            domtblout_line("NoClan", "PF88888.1", "PROT1", "1e-12", ali_from=150, ali_to=250),
            # non-overlapping same-clan hit -> kept
            domtblout_line("ResIII", "PF04851.1", "PROT2", "1e-20", ali_from=1, ali_to=100),
        ],
    )
    params = HmmscanParams(resolve_clan_overlap=True)
    df = filter_domains(parse_domtblout(path), params, 2, logger, clans=CLANS)
    prot1 = set(df[df.query_id == "PROT1"].pfam_acc)
    assert prot1 == {"PF00270", "PF99999", "PF88888"}
    assert df[df.query_id == "PROT2"].pfam_acc.tolist() == ["PF04851"]


def test_clan_resolution_requires_mapping(tmp_path):
    path = write_domtblout(tmp_path, [domtblout_line("DEAD", "PF00270.1", "PROT1", "1e-20")])
    with pytest.raises(ValueError, match="clan mapping"):
        filter_domains(parse_domtblout(path), HmmscanParams(resolve_clan_overlap=True), 1, logger, clans=None)


def test_resolve_clan_overlaps_direct(tmp_path):
    path = write_domtblout(
        tmp_path,
        [
            domtblout_line("ResIII", "PF04851.1", "PROT1", "1e-20", ali_from=135, ali_to=300),
            domtblout_line("DEAD", "PF00270.1", "PROT1", "9e-21", ali_from=137, ali_to=303),
        ],
    )
    kept = resolve_clan_overlaps(parse_domtblout(path), CLANS, logger)
    assert kept.pfam_acc.tolist() == ["PF00270"]


def test_parse_pfam_dat(tmp_path):
    dat = tmp_path / "Pfam-A.hmm.dat"
    dat.write_text(
        "# STOCKHOLM 1.0\n"
        "#=GF ID   1-cysPrx_C\n"
        "#=GF AC   PF10417.14\n"
        "#=GF CL   CL0172\n"
        "//\n"
        "# STOCKHOLM 1.0\n"
        "#=GF ID   NoClanFam\n"
        "#=GF AC   PF00001.20\n"
        "//\n"
    )
    assert parse_pfam_dat(dat) == {"PF10417": "CL0172"}
