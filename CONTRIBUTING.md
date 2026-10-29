# Contributing

## Changing a rule

1. Keep one rule per tool. If a tool needs a different input form (annotation,
   index, container), express it in that rule rather than in a shared script.
2. Any parameter that changes a count must be measurable. Add the measurement
   (sample, input size, runtime, peak memory) to the rule's comment or to
   [docs/parameters.md](docs/parameters.md).
3. Do not add an expression filter to the quantification or merge steps.
   Filtering belongs to analysis and differs per method.
4. Anything that touches the annotation must keep the audit in
   `rules/annotation.smk` passing: 693,078 ERV loci, 550 subfamilies, zero
   row-wise differences against the SQuIRE ERV subset.

## Before opening a change

```bash
.venv/bin/snakemake -n -p --cores 8 all     # must plan the expected job count
```

If the change affects a count, run one sample end to end and compare
`qc_per_sample.tsv` with the previous value.

## Style

- English in code, comments and documentation.
- Comments explain *why* a non-obvious step exists, including the failure it
  prevents. Do not narrate what the code already says.
- Python scripts in `workflow/scripts/` are standalone: they take positional
  arguments and must run without Snakemake.
