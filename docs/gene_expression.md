# Gene expression

Gene-level counts from the same BAMs used by the ERV tools.

## Flow

```
BAM (coordinate sorted)
  -> featureCounts            one <sample>_featureCounts.txt per sample
  -> merge_gene_counts.py     one wide matrix per cohort + one loader table
```

## Parameters

Identical to the TCASIA pipeline (`workflows/01_alignment`, rule `featurecounts`):

```bash
featureCounts -T 4 -p -t exon -g gene_id \
    -a gencode.v34.annotation.gtf -o <sample>_featureCounts.txt <sample>.bam
```

- `-p`: count fragments (read pairs), not reads
- `-t exon`, `-g gene_id`: gene level, using the GTF `gene_id`
- no `-s`: unstranded (`-s 0`). The historical matrices were built this way, so
  new samples stay comparable with them
- `-T 4`: four threads per sample; measured peak memory below 1 GB

## Outputs

| Path | Content |
|---|---|
| `gene_counts/<cohort>/<sample>/<sample>_featureCounts.txt` | raw featureCounts output |
| `gene_expression/matrices/<cohort>.counts.tsv.gz` | wide matrix, genes as rows (ENSG with version), samples of that cohort as columns |
| `gene_expression/gene_expression_standard.csv.gz` | loader table: `cohort,sample,gene_id,gene_symbol,length,count`, `count > 0` rows only |

## No expression filter

The workflow drops zero-count rows only, matching the ERV standard tables. It
applies no expression threshold, because filtering is an analysis-time choice:

- DESeq2 performs independent filtering internally; its vignette recommends
  `rowSums(counts) >= 10` as the minimal pre-filter
- edgeR / limma-voom expect `filterByExpr()` before modelling

Hard-coding a threshold in the platform would make those analyses
irreproducible, so the wide matrices keep all genes.

## Relation to the existing matrices

`/data/data2/Projects/TCASIA/Count/*_merged_counts.csv` was produced with the same
parameters (60,669 genes, versioned ENSG identifiers). New samples can be merged
with them directly: row names, column names and counting rules match.
