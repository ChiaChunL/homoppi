# homoppi

<img src="https://raw.githubusercontent.com/ChiaChunL/homoppi/main/docs/assets/homoppi_banner.png" alt="homoppi" width="100%">

| Testing | [![CI](https://github.com/ChiaChunL/homoppi/actions/workflows/ci.yml/badge.svg)](https://github.com/ChiaChunL/homoppi/actions/workflows/ci.yml) |
|---|---|
| Package | [![PyPI Latest Release](https://img.shields.io/pypi/v/homoppi)](https://pypi.org/project/homoppi/) [![Python versions](https://img.shields.io/pypi/pyversions/homoppi)](https://pypi.org/project/homoppi/) [![PyPI Downloads](https://img.shields.io/pypi/dm/homoppi)](https://pypi.org/project/homoppi/) |
| Meta | [![License - MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) |

## 🧬 What is it?

`homoppi` predicts protein-protein interactions (PPIs) for any species with a
proteome by transferring experimentally known interactions through homology:

- **Interolog mapping (IM)** — a query pair (A, B) is predicted to interact
  when homologs (A′, B′), found by blastp, form a known PPI in a template
  library.
- **Domain-domain interaction (DDI)** — a query pair is predicted to interact
  when their Pfam domains (annotated by hmmscan) form a known domain-domain
  interaction (3did, optionally expanded with an EM algorithm from known PPIs).

Each method integrates all matched templates into a probability,
`S = 1 − ∏(1 − s_template)`, and reports the full supporting evidence
(template pairs, species, alignment statistics) so every prediction is
traceable. The two methods complement each other: IM needs conserved
full-length homologs, while DDI still fires on fast-evolving proteins that
keep their domain architecture.

<img src="https://raw.githubusercontent.com/ChiaChunL/homoppi/main/docs/assets/homoppi_workflow.png" alt="homoppi workflow" width="100%">

## 📦 Installation

```bash
pip install homoppi
conda install -c bioconda blast hmmer   # external engines, found on PATH
```

Requires Python >= 3.11. Binaries can also be pointed to explicitly with
`--blastp-bin` / `--hmmscan-bin` / etc.

## ⚡ Quickstart

A tiny self-contained demo ships in [examples/](examples/) — from that
directory, every command below runs as-is:

```bash
# 1. build the reference database once
homoppi makedb --out db/ --ppi data/ppis.tsv --fasta 9606=data/templates_9606.fasta \
    --ddi data/ddis.tsv --pfam-hmm data/toy_pfam.hmm

# 2. score query pairs with both methods and merge
homoppi run --db db/ --workdir runs/demo --fasta data/query.fasta \
    --pairs data/pairs.tsv --fused
```

For real use, swap in your own libraries (see Input formats) and query
proteome. Batch mode scales to tens of thousands of pairs — unique proteins
are searched once, then all pairs are scored from the cached tables. Omitting
`--pairs` switches to proteome-wide all-vs-all mode. Each stage is also its
own command (`blast`, `domainanno`, `interolog`, `ddi`) with cached, resumable
results.

## 📄 Input formats

| File | Columns | Notes |
|---|---|---|
| PPI library (`--ppi`) | `protein_a  protein_b  taxid  score` | `score` optional, must be in [0, 1]; unscored templates fall back to `--default-template-score` (default 0) |
| Template proteomes (`--fasta TAXID=PATH`, repeatable) | FASTA | UniProt-style headers auto-parsed; else first token, or `--id-regex`; proteins absent from the PPI library are dropped |
| DDI library | `--ddi-3did` 3did flat file, `--ddi-em` EM scores, or `--ddi` pre-scored TSV `pfam_a  pfam_b  score` | 3did+EM scored as `½(S_known + S_EM)` |
| Pfam (`--pfam-hmm`, `--pfam-dat`) | `Pfam-A.hmm`, `Pfam-A.hmm.dat` | pressed into the db; `.dat` provides clans for `--resolve-clan-overlap` |
| Query pairs (`--pairs`) | 2 columns, header optional | IDs must match the query FASTA headers (first token) |

**Where the data comes from.** PPI templates: [IntAct](https://www.ebi.ac.uk/intact/)
ships ready-made MI confidence scores, or aggregate several databases
(BioGRID, HitPredict, ...) and score them with the
[HIPPIE](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/) scheme — either
way, normalize scores into [0, 1]. DDI templates:
[3did](https://3did.irbbarcelona.org/) flat file, optionally expanded with
`homoppi ddi-em`. Domain models: `Pfam-A.hmm` (+ `.hmm.dat`) from the
[EBI FTP](https://ftp.ebi.ac.uk/pub/databases/Pfam/).

## 📊 Outputs

`<workdir>/results/` per method (`interolog.*` / `ddi.*`), plus
`combined.summary.tsv` from `run`:

| Table | Columns |
|---|---|
| `*.summary.tsv` (one row per query pair) | `query_a, query_b, n_templates, s_im`/`s_ddi`, `best_template_a/_b`, `best_template_taxid`/`_source`, `best_template_score` |
| `interolog.evidence.tsv` (one row per template) | template pair, `taxid`, `template_score`, per-side `pident, qcov, scov, evalue` |
| `ddi.evidence.tsv` (one row per template) | domain pair, `source` (3did/em), `template_score`, per-side `cevalue` |
| `combined.summary.tsv` | `n_im_templates, s_im, n_ddi_templates, s_ddi` (+ `s_fused = 1−(1−s_im)(1−s_ddi)` with `--fused`) |
| `*.params.json` | snapshot of every effective parameter of the run |

Workdir layout: `blast/` and `hmmscan/` keep both raw and filtered search
results, `logs/` the external-tool logs, `state.json` the stage cache that
powers resume.

## 🎛 Thresholds and options

Defaults follow the published method: blastp homologs at identity ≥ 30%,
query coverage ≥ 40%, E-value ≤ 1e-10; hmmscan domains at conditional
E-value ≤ 1e-10 (or `--cut-tc`). Self pairs are discarded unless
`--include-self`; `--taxids` restricts IM evidence to chosen template species.

Optional strictness knobs, all off by default:

- `--min-subject-coverage` — also require the template side of a blastp hit
  to be covered (HSP-union fraction of the template length).
- `--min-hmm-coverage` — require a domain hit to span this fraction of the
  Pfam HMM model.
- `--resolve-clan-overlap` — overlapping hits from the same Pfam clan keep
  only the best one, preventing double-counted DDI evidence (needs
  `makedb --pfam-dat`).

Tightening thresholds on an existing workdir only re-filters the cached raw
output; loosening them (or changing the query) reruns the external search
automatically.

## 🧮 Expanding the DDI library with EM

Domain pairs co-occurring across known interacting proteins are scored by
expectation-maximization (Deng et al. 2002) and merged with 3did:

```bash
homoppi domainanno --db db/ --workdir runs/templates --fasta human.fasta
homoppi ddi-em --db db/ --domains runs/templates/hmmscan/domains.tsv --out em_scores.tsv
homoppi makedb --out db/ --ddi-3did 3did_flat --ddi-em em_scores.tsv
```

## ⏱ Performance

- Unique query proteins are searched once; scoring 1,000 pairs from cached
  tables takes seconds.
- Measured on 8 threads against a 7-species library (1.06 M scored PPIs,
  49 k proteins): 300 query proteins × 1,000 pairs through both methods in
  ~5 minutes (hmmscan dominates).
- Proteome-wide mode streams evidence to disk; memory scales with the number
  of predicted pairs, not the evidence volume.

## 📚 References

- Schaefer MH, Fontaine JF, Vinayagam A, Porras P, Wanker EE,
  Andrade-Navarro MA. [*HIPPIE: Integrating protein interaction networks with
  experiment based quality scores.*](https://doi.org/10.1371/journal.pone.0031826)
  PLoS ONE 7, e31826 (2012). — template PPI confidence scoring
- Orchard S et al. [*The MIntAct project — IntAct as a common curation
  platform for 11 molecular interaction
  databases.*](https://doi.org/10.1093/nar/gkt1115) Nucleic Acids Res 42,
  D358–D363 (2014). — scored template PPIs
- Yu H et al. [*Annotation transfer between genomes: protein-protein
  interologs and protein-DNA regulogs.*](https://doi.org/10.1101/gr.1774904)
  Genome Res 14, 1107–1118 (2004). — interolog mapping
- Mosca R, Céol A, Stein A, Olivella R, Aloy P. [*3did: a catalog of
  domain-based interactions of known three-dimensional
  structure.*](https://doi.org/10.1093/nar/gkt887) Nucleic Acids Res 42,
  D374–D379 (2014). — DDI templates
- Deng M, Mehta S, Sun F, Chen T. [*Inferring domain-domain interactions from
  protein-protein interactions.*](https://doi.org/10.1101/gr.153002) Genome
  Res 12, 1540–1548 (2002). — EM expansion of the DDI library
- Mistry J et al. [*Pfam: The protein families database in
  2021.*](https://doi.org/10.1093/nar/gkaa913) Nucleic Acids Res 49,
  D412–D419 (2021). — domain models

## 📄 License

[MIT](LICENSE)
