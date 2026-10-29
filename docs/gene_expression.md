# 基因表达

与 ERV 定量共用同一批 BAM，产出基因层面的 count 表。

## 流程

```
BAM (coordinate sorted)
  -> featureCounts           每样本一个 _featureCounts.txt
  -> merge_gene_counts.py    每队列一个宽矩阵 + 一张入库长表
```

## 参数

与 TCASIA 流水线（`workflows/01_alignment`，rule `featurecounts`）完全一致：

```bash
featureCounts -T 4 -p -t exon -g gene_id \
    -a gencode.v34.annotation.gtf -o <sample>_featureCounts.txt <sample>.bam
```

- `-p`：按 fragment（read pair）计数，不是按 read
- `-t exon`、`-g gene_id`：基因层面，用 GTF 的 `gene_id`
- 未指定 `-s`，即非链特异性（`-s 0`）。历史上的矩阵就是这样生成的，保持可比
- `-T 4`：每样本 4 线程；实测内存 < 1 GB

## 输出

| 文件 | 内容 |
|---|---|
| `gene_counts/<cohort>/<sample>/<sample>_featureCounts.txt` | featureCounts 原始输出 |
| `gene_expression/matrices/<cohort>.counts.tsv.gz` | 宽矩阵，行是基因（带版本号的 ENSG），列是该队列样本 |
| `gene_expression/gene_expression_standard.csv.gz` | 入库长表：`cohort,sample,gene_id,gene_symbol,length,count`，只含 count>0 的行 |

## 不做表达过滤

workflow 只去掉 count=0 的行（与 ERV 标准表一致），不做表达量阈值过滤。
过滤是分析阶段的选择：

- DESeq2：内部做 independent filtering；官方 vignette 建议的最小预过滤是 `rowSums(counts) >= 10`
- edgeR / limma-voom：分析前调用 `filterByExpr()`

在平台上固化一个阈值会让上述方法的结果无法复现，因此宽矩阵保留全部基因，
需要过滤的分析自行处理。

## 与已有矩阵的关系

`/data/data2/Projects/TCASIA/Count/*_merged_counts.csv` 是同一套参数生成的历史矩阵
（60,669 个基因，ENSG 带版本号）。本流程产出的新样本可以直接与它们合并，
行名、列名和计数口径都一致。
