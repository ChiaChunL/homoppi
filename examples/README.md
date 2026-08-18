# homoppi example

A tiny self-contained demo (toy sequences, toy libraries). Requires BLAST+ and
HMMER on PATH. From this directory:

```bash
homoppi makedb --out db/ --ppi data/ppis.tsv --fasta 9606=data/templates_9606.fasta \
    --ddi data/ddis.tsv --pfam-hmm data/toy_pfam.hmm

homoppi run --db db/ --workdir runs/demo --fasta data/query.fasta --pairs data/pairs.tsv --fused

column -t runs/demo/results/combined.summary.tsv
```

Expected result: the query pair A-B is supported by the template PPI T1-T2
(s_im = 0.8) and by the domain pair DOMA-DOMB (s_ddi = 0.7), fused to 0.94.
