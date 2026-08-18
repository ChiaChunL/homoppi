"""Database building: normalize template FASTAs, build the BLAST DB, write the manifest."""

from __future__ import annotations

import json
import logging
import re
from datetime import UTC, datetime
from pathlib import Path

from tqdm import tqdm

from . import __version__
from .blast import PROTEIN2TAXID_TSV, blastdb_prefix
from .external import find_binary, run_command
from .fasta import extract_id, read_fasta
from .ppidb import build_ppi_component, lib_proteins_by_taxid


def parse_fasta_specs(specs: list[str]) -> list[tuple[int, Path]]:
    """Parse repeated --fasta TAXID=PATH options."""
    parsed: list[tuple[int, Path]] = []
    for spec in specs:
        taxid_str, sep, path_str = spec.partition("=")
        if not sep or not taxid_str.strip().isdigit():
            raise ValueError(f"--fasta expects TAXID=PATH (e.g. 9606=human.fasta), got: {spec!r}")
        path = Path(path_str).expanduser()
        if not path.exists():
            raise FileNotFoundError(f"FASTA not found: {path}")
        parsed.append((int(taxid_str), path))
    return parsed


def build_blast_component(
    fastas: list[tuple[int, Path]],
    db_dir: Path,
    lib_proteins: dict[int, set[str]],
    logger: logging.Logger,
    id_regex: str | None = None,
    keep_all_proteins: bool = False,
    makeblastdb_bin: Path | None = None,
) -> dict:
    """Normalize headers, restrict to library proteins, and run makeblastdb.

    By default only proteins that occur in the PPI template library are kept:
    hits outside the library can never support a template, and the smaller DB
    makes blastp much faster. Disable with keep_all_proteins.
    """
    regex = re.compile(id_regex) if id_regex else None
    blast_dir = db_dir / "blastdb"
    blast_dir.mkdir(parents=True, exist_ok=True)

    wanted_total = {p for proteins in lib_proteins.values() for p in proteins}
    protein2taxid: dict[str, int] = {}
    n_written = 0
    n_dupes = 0

    combined_fasta = blast_dir / "templates.fasta"
    with open(combined_fasta, "w") as out:
        for taxid, fasta_path in fastas:
            n_taxid = 0
            for header, seq in tqdm(read_fasta(fasta_path), desc=f"Reading {fasta_path.name}", unit="seq"):
                acc = extract_id(header, regex)
                if not keep_all_proteins and acc not in wanted_total:
                    continue
                if acc in protein2taxid:
                    n_dupes += 1
                    continue
                protein2taxid[acc] = taxid
                out.write(f">{acc}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i : i + 60] + "\n")
                n_written += 1
                n_taxid += 1
            logger.info("taxid %s: kept %s sequences from %s", taxid, f"{n_taxid:,}", fasta_path.name)
    if n_dupes:
        logger.warning(
            "%s duplicate accessions across FASTA files were skipped (first occurrence kept).", f"{n_dupes:,}"
        )
    if n_written == 0:
        raise ValueError(
            "No template sequences were kept. Check that FASTA IDs match the protein IDs used in the PPI library "
            "(UniProt-style 'sp|ACC|NAME' headers are handled automatically; otherwise pass --id-regex)."
        )

    with open(blast_dir / PROTEIN2TAXID_TSV, "w") as fh:
        fh.write("protein_id\ttaxid\n")
        for protein, taxid in protein2taxid.items():
            fh.write(f"{protein}\t{taxid}\n")

    makeblastdb = find_binary("makeblastdb", makeblastdb_bin)
    run_command(
        [
            makeblastdb,
            "-in", str(combined_fasta),
            "-dbtype", "prot",
            "-out", str(blastdb_prefix(db_dir)),
            "-title", "homoppi_templates",
        ],
        db_dir / "makeblastdb.log",
        logger,
    )

    # Coverage report: how many library proteins actually have a sequence.
    coverage: dict[int, dict] = {}
    for taxid, proteins in sorted(lib_proteins.items()):
        found = sum(1 for p in proteins if p in protein2taxid)
        pct = 100 * found / len(proteins) if proteins else 0.0
        coverage[taxid] = {"library_proteins": len(proteins), "with_sequence": found, "pct": round(pct, 1)}
        level = logging.WARNING if pct < 90 else logging.INFO
        logger.log(level, "taxid %-8s %s/%s library proteins have a sequence (%.1f%%)",
                   taxid, f"{found:,}", f"{len(proteins):,}", pct)

    return {"sequences": n_written, "coverage_by_taxid": coverage}


def build_pfam_component(
    pfam_hmm: Path,
    db_dir: Path,
    logger: logging.Logger,
    hmmpress_bin: Path | None = None,
) -> dict:
    """Link the Pfam HMM into the database and press it for hmmscan.

    hmmpress compiles the text HMM file into binary indices (.h3m/.h3i/.h3f/.h3p)
    that hmmscan requires. The source file is symlinked, not copied; the pressed
    indices live inside the database directory.
    """
    from .hmmer import PFAM_HMM

    pfam_dir = db_dir / "pfam"
    pfam_dir.mkdir(parents=True, exist_ok=True)
    link = pfam_dir / PFAM_HMM
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(pfam_hmm.resolve())

    hmmpress = find_binary("hmmpress", hmmpress_bin)
    run_command([hmmpress, "-f", str(link)], db_dir / "hmmpress.log", logger)
    logger.info("Pfam HMM pressed into %s", pfam_dir)
    return {"source": str(pfam_hmm.resolve()), "size_bytes": pfam_hmm.stat().st_size}


def build_clan_component(pfam_dat: Path, db_dir: Path, logger: logging.Logger) -> dict:
    """Extract the family -> clan mapping (for --resolve-clan-overlap) from Pfam-A.hmm.dat."""
    from .hmmer import clans_path, parse_pfam_dat

    clans = parse_pfam_dat(pfam_dat)
    if not clans:
        raise ValueError(f"no family->clan records found in {pfam_dat}; is this a Pfam-A.hmm.dat file?")
    out_path = clans_path(db_dir)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("pfam_acc\tclan\n")
        for acc, clan in sorted(clans.items()):
            fh.write(f"{acc}\t{clan}\n")
    logger.info("Clan mapping for %s Pfam families written to %s", f"{len(clans):,}", out_path)
    return {"source": str(pfam_dat.resolve()), "families_with_clan": len(clans)}


def write_manifest(db_dir: Path, components: dict, inputs: dict) -> None:
    manifest_path = db_dir / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        manifest.get("components", {}).update(components)
        manifest.get("inputs", {}).update(inputs)
        components, inputs = manifest["components"], manifest["inputs"]
    manifest = {
        "homoppi_version": __version__,
        "updated_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "inputs": inputs,
        "components": components,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")


def build_database(
    db_dir: Path,
    ppi_tsv: Path | None,
    fasta_specs: list[str] | None,
    logger: logging.Logger,
    id_regex: str | None = None,
    keep_all_proteins: bool = False,
    makeblastdb_bin: Path | None = None,
    ddi_3did: Path | None = None,
    ddi_em: Path | None = None,
    ddi_plain: Path | None = None,
    pfam_hmm: Path | None = None,
    pfam_dat: Path | None = None,
    hmmpress_bin: Path | None = None,
) -> None:
    """Build the requested database components; existing ones are left in place."""
    from .ddidb import build_ddi_component

    if (ppi_tsv is None) != (not fasta_specs):
        raise ValueError("--ppi and --fasta must be given together (the BLAST DB is built from library proteins).")
    wants_ddi = any(x is not None for x in (ddi_3did, ddi_em, ddi_plain))
    if ppi_tsv is None and not wants_ddi and pfam_hmm is None and pfam_dat is None:
        raise ValueError("nothing to build: provide --ppi/--fasta, DDI inputs, --pfam-hmm, and/or --pfam-dat.")

    db_dir.mkdir(parents=True, exist_ok=True)
    components: dict = {}
    inputs: dict = {}

    if ppi_tsv is not None:
        fastas = parse_fasta_specs(fasta_specs or [])
        components["ppi"] = build_ppi_component(ppi_tsv, db_dir / "ppi", logger)
        lib_proteins = lib_proteins_by_taxid(db_dir / "ppi")
        components["blast"] = build_blast_component(
            fastas, db_dir, lib_proteins, logger,
            id_regex=id_regex, keep_all_proteins=keep_all_proteins, makeblastdb_bin=makeblastdb_bin,
        )
        inputs["ppi"] = {
            "ppi_tsv": str(ppi_tsv.resolve()),
            "fastas": [{"taxid": t, "path": str(p.resolve())} for t, p in fastas],
            "id_regex": id_regex,
            "keep_all_proteins": keep_all_proteins,
        }

    if wants_ddi:
        components["ddi"] = build_ddi_component(
            db_dir / "ddi", logger, three_did=ddi_3did, em_scores=ddi_em, plain=ddi_plain
        )
        inputs["ddi"] = {
            "ddi_3did": str(ddi_3did.resolve()) if ddi_3did else None,
            "ddi_em": str(ddi_em.resolve()) if ddi_em else None,
            "ddi_plain": str(ddi_plain.resolve()) if ddi_plain else None,
        }

    if pfam_hmm is not None:
        components["pfam"] = build_pfam_component(pfam_hmm, db_dir, logger, hmmpress_bin=hmmpress_bin)
        inputs["pfam"] = {"pfam_hmm": str(pfam_hmm.resolve())}

    if pfam_dat is not None:
        components["clans"] = build_clan_component(pfam_dat, db_dir, logger)
        inputs["clans"] = {"pfam_dat": str(pfam_dat.resolve())}

    write_manifest(db_dir, components=components, inputs=inputs)
    logger.info("Database ready at %s", db_dir)
