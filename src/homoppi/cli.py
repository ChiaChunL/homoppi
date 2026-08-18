"""homoppi command-line interface."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from . import __version__
from .config import DEFAULT_TEMPLATE_SCORE, BlastParams, HmmscanParams, parse_taxids
from .log import get_logger

app = typer.Typer(
    name="homoppi",
    help=(
        "Homology-based protein-protein interaction prediction: "
        "interolog mapping (IM) and domain-domain interaction (DDI) inference."
    ),
    no_args_is_help=True,
    add_completion=False,
)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"homoppi {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool | None,
        typer.Option("--version", callback=_version_callback, is_eager=True, help="Show version and exit."),
    ] = None,
) -> None:
    """Stages: makedb -> blast/domainanno -> interolog/ddi (or `run` for the whole pipeline).

    Omit --pairs for proteome-wide all-vs-all mode.
    """


# --------------------------------------------------------------------------- makedb

@app.command()
def makedb(
    out: Annotated[Path, typer.Option(help="Output database directory.")],
    ppi: Annotated[Path | None, typer.Option(exists=True, dir_okay=False, help="PPI template library TSV: protein_a, protein_b, taxid[, score].")] = None,
    fasta: Annotated[list[str] | None, typer.Option(help="Template proteome as TAXID=PATH (repeatable), e.g. --fasta 9606=human.fasta.")] = None,
    ddi_3did: Annotated[Path | None, typer.Option("--ddi-3did", exists=True, dir_okay=False, help="3did flat file (3did_flat) to build the DDI library from.")] = None,
    ddi_em: Annotated[Path | None, typer.Option("--ddi-em", exists=True, dir_okay=False, help="EM score TSV produced by `homoppi ddi-em` (pfam_a, pfam_b, em_score).")] = None,
    ddi: Annotated[Path | None, typer.Option(exists=True, dir_okay=False, help="Pre-scored DDI library TSV: pfam_a, pfam_b[, score]. Mutually exclusive with --ddi-3did/--ddi-em.")] = None,
    pfam_hmm: Annotated[Path | None, typer.Option("--pfam-hmm", exists=True, dir_okay=False, help="Pfam-A.hmm file; symlinked into the database and pressed with hmmpress.")] = None,
    pfam_dat: Annotated[Path | None, typer.Option("--pfam-dat", exists=True, dir_okay=False, help="Pfam-A.hmm.dat file; provides the family->clan mapping needed by --resolve-clan-overlap.")] = None,
    id_regex: Annotated[str | None, typer.Option(help="Regex with one capture group to extract protein IDs from FASTA headers (UniProt-style headers are auto-detected).")] = None,
    keep_all_proteins: Annotated[bool, typer.Option("--keep-all-proteins", help="Keep FASTA proteins absent from the PPI library (default: drop them; they cannot support templates).")] = False,
    makeblastdb_bin: Annotated[Path | None, typer.Option(help="Path to makeblastdb (default: search PATH).")] = None,
    hmmpress_bin: Annotated[Path | None, typer.Option(help="Path to hmmpress (default: search PATH).")] = None,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
) -> None:
    """Build database components: PPI library + BLAST DB, DDI library, pressed Pfam HMMs.

    Components are independent; rebuild any subset without touching the others.
    """
    from .makedb import build_database

    logger = get_logger(verbose=verbose)
    build_database(
        db_dir=out, ppi_tsv=ppi, fasta_specs=fasta, logger=logger,
        id_regex=id_regex, keep_all_proteins=keep_all_proteins, makeblastdb_bin=makeblastdb_bin,
        ddi_3did=ddi_3did, ddi_em=ddi_em, ddi_plain=ddi,
        pfam_hmm=pfam_hmm, pfam_dat=pfam_dat, hmmpress_bin=hmmpress_bin,
    )


# --------------------------------------------------------------------------- evidence preparation

_DB_OPT = typer.Option(exists=True, file_okay=False, help="Database directory built by `homoppi makedb`.")
_WD_OPT = typer.Option(help="Working directory for this run (created if missing).")


@app.command()
def blast(
    db: Annotated[Path, _DB_OPT],
    workdir: Annotated[Path, _WD_OPT],
    fasta: Annotated[Path, typer.Option(exists=True, dir_okay=False, help="Query proteome/protein FASTA.")],
    evalue: Annotated[float, typer.Option(help="Maximum E-value.")] = 1e-10,
    min_identity: Annotated[float, typer.Option(help="Minimum sequence identity (%).")] = 30.0,
    min_coverage: Annotated[float, typer.Option(help="Minimum query coverage (%).")] = 40.0,
    min_subject_coverage: Annotated[float, typer.Option(help="Minimum template (subject) coverage (%); 0 disables.")] = 0.0,
    max_target_seqs: Annotated[int, typer.Option(help="blastp -max_target_seqs.")] = 500,
    threads: Annotated[int, typer.Option("--threads", "-t")] = 4,
    blastp_bin: Annotated[Path | None, typer.Option(help="Path to blastp (default: search PATH).")] = None,
    force: Annotated[bool, typer.Option("--force", help="Rerun even if the stage is up to date.")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
) -> None:
    """Find template homologs of the query proteins (blastp + threshold filtering)."""
    from .blast import ensure_blast_stage
    from .workdir import Workdir

    wd = Workdir(workdir)
    logger = get_logger(log_file=wd.logs_dir / "homoppi.log", verbose=verbose)
    params = BlastParams(
        evalue=evalue, min_identity=min_identity, min_coverage=min_coverage,
        min_subject_coverage=min_subject_coverage, max_target_seqs=max_target_seqs, threads=threads,
    )
    ensure_blast_stage(wd, db, fasta, params, logger, force=force, blastp_bin=blastp_bin)


@app.command()
def domainanno(
    db: Annotated[Path, _DB_OPT],
    workdir: Annotated[Path, _WD_OPT],
    fasta: Annotated[Path, typer.Option(exists=True, dir_okay=False, help="Query proteome/protein FASTA.")],
    cevalue: Annotated[float, typer.Option(help="Maximum conditional (domain) E-value.")] = 1e-10,
    cut_tc: Annotated[bool, typer.Option("--cut-tc", help="Use Pfam trusted cutoffs instead of E-value filtering.")] = False,
    min_hmm_coverage: Annotated[float, typer.Option(help="Minimum fraction of the HMM model matched (%); 0 disables.")] = 0.0,
    resolve_clan_overlap: Annotated[bool, typer.Option("--resolve-clan-overlap", help="Drop overlapping hits from the same Pfam clan, keeping the best (needs `makedb --pfam-dat`).")] = False,
    threads: Annotated[int, typer.Option("--threads", "-t")] = 4,
    hmmscan_bin: Annotated[Path | None, typer.Option(help="Path to hmmscan (default: search PATH).")] = None,
    force: Annotated[bool, typer.Option("--force", help="Rerun even if the stage is up to date.")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
) -> None:
    """Annotate Pfam domains on the query proteins (hmmscan + c-Evalue filtering)."""
    from .hmmer import ensure_hmmscan_stage
    from .workdir import Workdir

    wd = Workdir(workdir)
    logger = get_logger(log_file=wd.logs_dir / "homoppi.log", verbose=verbose)
    params = HmmscanParams(
        cevalue=cevalue, cut_tc=cut_tc, min_hmm_coverage=min_hmm_coverage,
        resolve_clan_overlap=resolve_clan_overlap, threads=threads,
    )
    ensure_hmmscan_stage(wd, db, fasta, params, logger, force=force, hmmscan_bin=hmmscan_bin)


# --------------------------------------------------------------------------- inference

@app.command()
def interolog(
    db: Annotated[Path, _DB_OPT],
    workdir: Annotated[Path, _WD_OPT],
    pairs: Annotated[Path | None, typer.Option(exists=True, dir_okay=False, help="Query pairs TSV (2 columns). Omit for proteome-wide all-vs-all mode.")] = None,
    fasta: Annotated[Path | None, typer.Option(exists=True, dir_okay=False, help="Query FASTA; triggers the blast stage automatically when its results are missing or stale.")] = None,
    taxids: Annotated[str | None, typer.Option(help="Restrict template evidence to these taxids, e.g. '9606,10090'.")] = None,
    include_self: Annotated[bool, typer.Option("--include-self", help="Score self pairs (A == A); discarded by default.")] = False,
    default_template_score: Annotated[float, typer.Option(help="Score assumed for templates without a library score.")] = DEFAULT_TEMPLATE_SCORE,
    min_score: Annotated[float, typer.Option(help="Proteome-wide mode only: drop predictions below this score.")] = 0.0,
    no_evidence: Annotated[bool, typer.Option("--no-evidence", help="Skip the per-template evidence table (smaller output).")] = False,
    evalue: Annotated[float, typer.Option(help="blast stage: maximum E-value.")] = 1e-10,
    min_identity: Annotated[float, typer.Option(help="blast stage: minimum identity (%).")] = 30.0,
    min_coverage: Annotated[float, typer.Option(help="blast stage: minimum query coverage (%).")] = 40.0,
    min_subject_coverage: Annotated[float, typer.Option(help="blast stage: minimum template (subject) coverage (%); 0 disables.")] = 0.0,
    max_target_seqs: Annotated[int, typer.Option(help="blast stage: blastp -max_target_seqs.")] = 500,
    threads: Annotated[int, typer.Option("--threads", "-t")] = 4,
    blastp_bin: Annotated[Path | None, typer.Option(help="Path to blastp (default: search PATH).")] = None,
    force: Annotated[bool, typer.Option("--force", help="Rerun the blast stage even if up to date.")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
) -> None:
    """Predict interactions by interolog mapping (batch pairs or proteome-wide)."""
    from .pipeline import interolog_stage
    from .workdir import Workdir

    wd = Workdir(workdir)
    logger = get_logger(log_file=wd.logs_dir / "homoppi.log", verbose=verbose)
    interolog_stage(
        wd, db, logger,
        pairs_path=pairs, fasta=fasta, taxid_set=parse_taxids(taxids),
        include_self=include_self, default_score=default_template_score, min_score=min_score,
        no_evidence=no_evidence,
        blast_params=BlastParams(
            evalue=evalue, min_identity=min_identity, min_coverage=min_coverage,
            min_subject_coverage=min_subject_coverage, max_target_seqs=max_target_seqs, threads=threads,
        ),
        blastp_bin=blastp_bin, force=force,
    )


@app.command()
def ddi(
    db: Annotated[Path, _DB_OPT],
    workdir: Annotated[Path, _WD_OPT],
    pairs: Annotated[Path | None, typer.Option(exists=True, dir_okay=False, help="Query pairs TSV (2 columns). Omit for proteome-wide all-vs-all mode.")] = None,
    fasta: Annotated[Path | None, typer.Option(exists=True, dir_okay=False, help="Query FASTA; triggers the hmmscan stage automatically when its results are missing or stale.")] = None,
    include_self: Annotated[bool, typer.Option("--include-self", help="Score self pairs (A == A); discarded by default.")] = False,
    default_template_score: Annotated[float, typer.Option(help="Score assumed for templates without a library score.")] = DEFAULT_TEMPLATE_SCORE,
    min_score: Annotated[float, typer.Option(help="Proteome-wide mode only: drop predictions below this score.")] = 0.0,
    no_evidence: Annotated[bool, typer.Option("--no-evidence", help="Skip the per-template evidence table (smaller output).")] = False,
    cevalue: Annotated[float, typer.Option(help="hmmscan stage: maximum conditional E-value.")] = 1e-10,
    cut_tc: Annotated[bool, typer.Option("--cut-tc", help="hmmscan stage: use Pfam trusted cutoffs.")] = False,
    min_hmm_coverage: Annotated[float, typer.Option(help="hmmscan stage: minimum fraction of the HMM model matched (%); 0 disables.")] = 0.0,
    resolve_clan_overlap: Annotated[bool, typer.Option("--resolve-clan-overlap", help="hmmscan stage: drop overlapping hits from the same Pfam clan (needs `makedb --pfam-dat`).")] = False,
    threads: Annotated[int, typer.Option("--threads", "-t")] = 4,
    hmmscan_bin: Annotated[Path | None, typer.Option(help="Path to hmmscan (default: search PATH).")] = None,
    force: Annotated[bool, typer.Option("--force", help="Rerun the hmmscan stage even if up to date.")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
) -> None:
    """Predict interactions by domain-domain interaction inference (batch pairs or proteome-wide)."""
    from .pipeline import ddi_stage
    from .workdir import Workdir

    wd = Workdir(workdir)
    logger = get_logger(log_file=wd.logs_dir / "homoppi.log", verbose=verbose)
    ddi_stage(
        wd, db, logger,
        pairs_path=pairs, fasta=fasta,
        include_self=include_self, default_score=default_template_score, min_score=min_score,
        no_evidence=no_evidence,
        hmm_params=HmmscanParams(
            cevalue=cevalue, cut_tc=cut_tc, min_hmm_coverage=min_hmm_coverage,
            resolve_clan_overlap=resolve_clan_overlap, threads=threads,
        ),
        hmmscan_bin=hmmscan_bin, force=force,
    )


# --------------------------------------------------------------------------- ddi-em

@app.command(name="ddi-em")
def ddi_em(
    db: Annotated[Path, _DB_OPT],
    domains: Annotated[list[Path], typer.Option(exists=True, dir_okay=False, help="Domain annotation TSV(s) for the template proteins (domainanno output; repeatable).")],
    out: Annotated[Path, typer.Option(help="Output TSV of EM scores (pfam_a, pfam_b, em_score, ...).")],
    false_neg_rate: Annotated[float, typer.Option(help="Assumed false-negative rate of the PPI library.")] = 0.8,
    false_pos_rate: Annotated[float | None, typer.Option(help="False-positive rate (default: derived from library counts).")] = None,
    max_iter: Annotated[int, typer.Option(help="Maximum EM iterations.")] = 50,
    tol: Annotated[float, typer.Option(help="Convergence threshold on max |delta lambda|.")] = 1e-4,
    max_proteins_per_domain: Annotated[int, typer.Option(help="Exclude domains present in more proteins than this.")] = 1000,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
) -> None:
    """Expand the DDI library: EM-score candidate domain pairs from the known PPIs.

    Feed the output into `homoppi makedb --ddi-3did ... --ddi-em <out>`.
    """
    from .blast import load_protein2taxid
    from .em import EMParams, load_domain_annotations, run_ddi_em
    from .ppidb import PPIIndex

    logger = get_logger(verbose=verbose)
    index = PPIIndex.load(db)
    protein_taxid = load_protein2taxid(db)
    annotations = load_domain_annotations(domains)
    logger.info(
        "Inputs: %s known PPIs, %s proteins with taxid, %s proteins with domains",
        f"{len(index.pairs):,}", f"{len(protein_taxid):,}", f"{len(annotations):,}",
    )

    params = EMParams(
        false_neg_rate=false_neg_rate, false_pos_rate=false_pos_rate,
        max_iter=max_iter, tol=tol, max_proteins_per_domain=max_proteins_per_domain,
    )
    result = run_ddi_em(set(index.pairs), protein_taxid, annotations, params, logger)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, sep="\t", index=False)
    logger.info("EM scores for %s domain pairs written to %s", f"{len(result):,}", out)


# --------------------------------------------------------------------------- run

@app.command()
def run(
    db: Annotated[Path, _DB_OPT],
    workdir: Annotated[Path, _WD_OPT],
    fasta: Annotated[Path, typer.Option(exists=True, dir_okay=False, help="Query proteome/protein FASTA.")],
    pairs: Annotated[Path | None, typer.Option(exists=True, dir_okay=False, help="Query pairs TSV. Omit for proteome-wide mode.")] = None,
    skip_im: Annotated[bool, typer.Option("--skip-im", help="Skip the interolog method.")] = False,
    skip_ddi: Annotated[bool, typer.Option("--skip-ddi", help="Skip the DDI method.")] = False,
    fused: Annotated[bool, typer.Option("--fused", help="Add an s_fused column: 1-(1-s_im)(1-s_ddi).")] = False,
    taxids: Annotated[str | None, typer.Option(help="IM only: restrict template evidence to these taxids.")] = None,
    include_self: Annotated[bool, typer.Option("--include-self", help="Score self pairs (A == A); discarded by default.")] = False,
    default_template_score: Annotated[float, typer.Option(help="Score assumed for templates without a library score.")] = DEFAULT_TEMPLATE_SCORE,
    min_score: Annotated[float, typer.Option(help="Proteome-wide mode only: drop predictions below this score.")] = 0.0,
    no_evidence: Annotated[bool, typer.Option("--no-evidence", help="Skip the per-template evidence tables.")] = False,
    evalue: Annotated[float, typer.Option(help="blast stage: maximum E-value.")] = 1e-10,
    min_identity: Annotated[float, typer.Option(help="blast stage: minimum identity (%).")] = 30.0,
    min_coverage: Annotated[float, typer.Option(help="blast stage: minimum query coverage (%).")] = 40.0,
    min_subject_coverage: Annotated[float, typer.Option(help="blast stage: minimum template (subject) coverage (%); 0 disables.")] = 0.0,
    max_target_seqs: Annotated[int, typer.Option(help="blast stage: blastp -max_target_seqs.")] = 500,
    cevalue: Annotated[float, typer.Option(help="hmmscan stage: maximum conditional E-value.")] = 1e-10,
    cut_tc: Annotated[bool, typer.Option("--cut-tc", help="hmmscan stage: use Pfam trusted cutoffs.")] = False,
    min_hmm_coverage: Annotated[float, typer.Option(help="hmmscan stage: minimum fraction of the HMM model matched (%); 0 disables.")] = 0.0,
    resolve_clan_overlap: Annotated[bool, typer.Option("--resolve-clan-overlap", help="hmmscan stage: drop overlapping hits from the same Pfam clan (needs `makedb --pfam-dat`).")] = False,
    threads: Annotated[int, typer.Option("--threads", "-t")] = 4,
    blastp_bin: Annotated[Path | None, typer.Option(help="Path to blastp (default: search PATH).")] = None,
    hmmscan_bin: Annotated[Path | None, typer.Option(help="Path to hmmscan (default: search PATH).")] = None,
    force: Annotated[bool, typer.Option("--force", help="Rerun cached stages.")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v")] = False,
) -> None:
    """Run the full pipeline: blast -> interolog, hmmscan -> ddi, then merge.

    Completed stages are cached in the workdir; an interrupted run resumes
    where it stopped.
    """
    from .pipeline import ddi_stage, interolog_stage, merge_stage
    from .workdir import Workdir

    if skip_im and skip_ddi:
        raise typer.BadParameter("--skip-im and --skip-ddi cannot both be set.")

    wd = Workdir(workdir)
    logger = get_logger(log_file=wd.logs_dir / "homoppi.log", verbose=verbose)

    if not skip_im:
        interolog_stage(
            wd, db, logger,
            pairs_path=pairs, fasta=fasta, taxid_set=parse_taxids(taxids),
            include_self=include_self, default_score=default_template_score, min_score=min_score,
            no_evidence=no_evidence,
            blast_params=BlastParams(
                evalue=evalue, min_identity=min_identity, min_coverage=min_coverage,
                min_subject_coverage=min_subject_coverage, max_target_seqs=max_target_seqs, threads=threads,
            ),
            blastp_bin=blastp_bin, force=force,
        )
    if not skip_ddi:
        ddi_stage(
            wd, db, logger,
            pairs_path=pairs, fasta=fasta,
            include_self=include_self, default_score=default_template_score, min_score=min_score,
            no_evidence=no_evidence,
            hmm_params=HmmscanParams(
                cevalue=cevalue, cut_tc=cut_tc, min_hmm_coverage=min_hmm_coverage,
                resolve_clan_overlap=resolve_clan_overlap, threads=threads,
            ),
            hmmscan_bin=hmmscan_bin, force=force,
        )
    if not skip_im and not skip_ddi:
        merge_stage(wd, logger, fused=fused)


if __name__ == "__main__":
    app()
