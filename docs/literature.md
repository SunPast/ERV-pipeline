# Literature notes: annotation handling in ERV/TE quantification tools

> Purpose: confirm "should cross-tool comparisons use a unified annotation" and "does changing the annotation change the algorithm or affect detection".
> Every conclusion is backed by the original sentence for easy verification. Search date 2026-10-27.

## One-sentence conclusion

**A unified annotation is the standard practice in the field for cross-tool comparisons, and has been explicitly implemented by at least two benchmarks**; changing the annotation does not change the algorithm, only the assignable space.

---

## 1. Bendall et al. 2019, PLoS Comput Biol 15(9):e1006453 (Telescope paper)

- Link: https://doi.org/10.1371/journal.pcbi.1006453
- Annotation source: the authors built 14,968 HERV loci / 60 subfamilies, repository `mlbendall/telescope_annotation_db`

Key original text (simulation benchmark):

> "To ensure equal comparisons, **all approaches use the same annotation** (S1 File), and modifications to the annotation were made to allow locus-specific quantification (instead of family-level quantification) for RepEnrich, TEtranscripts, and SalmonTE."

Methods section:

> "To ensure a fair comparison among approaches, **the same annotation (S1 File) was used as input for all approaches**... For RepEnrich, TEtranscripts, and SalmonTE, **the locus identifier was used in place of the family name** in order to generate locus-specific estimates."

Real data section (note this is the opposite of the above):

> "RepEnrich, TEtranscripts, and SalmonTE were all run according to author instructions, **with author-provided annotations** and default parameters."

S4 figure caption: `using author provided annotations and default arguments`

**Implication**: the original authors unified the annotation in simulation but not in real data; the latter can only be judged by "whether clustering recovers cell types". This shows that "unified annotation" is a precondition for comparability, while "each uses its own" makes the values not directly comparable.

---

## 2. Schwarz et al. 2022, Brief Bioinform 23(1):bbab417

- Link: https://doi.org/10.1093/bib/bbab417
- Tools: SalmonTE\*, TEtools\*, TEtranscripts\*, SQuIRE, Telescope (five tools at locus level)

Annotation: **all unified** to the same RepeatMasker `.align` (hg38, Repeat Library 20140131)

> "We used the repeat annotation of RepeatMasker of M. musculus (mm10 ...), H. sapiens (**hg38, based on Repeat Library 20140131** ...)"

How family-level tools were turned into locus-level tools (only the ids changed, not the code):

> "TEtools aligns reads to the instance-specific reference sequences and aggregates ... To suppress the aggregation step, **we substituted the ids of TE families with ids of TE instances**. Similarly, TEtranscripts uses an annotation file mapping TEs to their respective families. **Again, we substituted the family names by TE ids** to avoid the aggregation process."

SQuIRE's annotation entry point (key: external annotation is usable):

> "SQuIRE requires RepeatMasker's '.out' file format. To provide such a file, **we translated the downloaded '.align' files into the '.out' format via an in-house script. This mapping is bijective** ... From this, SQuIRE generates its own annotation file in BED format with SQuIRE-specific TE ids."

Explicit statement that "changing the annotation is not changing the algorithm":

> "As we did **not change the algorithm** of the tools, which are responsible for the assignment of multi-mapping reads, we do not expect an interfering of the outcomes."

On Telescope's read filtering (we need to verify this empirically):

> "By default, reads and read-pairs mapping to **more than 100 possible loci** are removed."
> — Note: this filter implementation was not found in the local Telescope 1.0.3 source (`telescope_assign.py` has no such parameter); trust the empirical measurement.

Tool ranking (Table 2):

| Tool | Detection | Quantification | Differential expression |
|---|---|---|---|
| SalmonTE\* | 1 | 1 | 1 |
| Telescope | 2 | 4 | 2 |
| SQuIRE | 3 | 2 | 5 |
| TEtools\* | 5 | 3 | 3 |
| TEtranscripts\* | 4 | 5 | 4 |

---

## 3. Savytska et al. 2022, Front Genet 13:1026847

- Link: https://doi.org/10.3389/fgene.2022.1026847
- Tools: SQuIRE, TElocal, SalmonTE, featureCounts (unique/fraction/random)

Locus-level false-positive conclusion:

> "Our findings suggest a **high false discovery number that exceeds the total number of correctly recovered active loci** for all the quantification strategies, including the best performing tool TElocal."

Specific numbers:

> "While expression was simulated for only 491 SVA loci for stranded data, **TElocal_UM detected 1455 SVA loci** as putatively active, including 67 loci for SVA_B that had no simulated loci and 911 loci for SVA_D that had only 4 loci simulated."

**Implication**: locus-level "positive detection" is inherently very noisy and must be paired with thresholds; this also means locus-level positive counts cannot be used directly as "expression quantities".

---

## 4. TE-Seq 2025 (PMC11507816)

- Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC11507816/

Why Telescope was chosen:

> "Our TE-Seq pipeline utilizes the Telescope tool since **it accepts user-provided gene-transfer format (GTF) formatted TE annotations**, making it more flexible than either TElocal or L1EM."

Annotation:

> "Repetitive element count estimates are produced with the Telescope tool and the **custom RepeatMasker derived annotation** resulting from the AREF module."

**Implication**: feeding a custom rmsk annotation to Telescope is a published standard usage, not our ad-hoc invention. It also outputs both "conservative (unique-only)" and "permissive (including multi-mapping)" count sets, which is worth borrowing.

---

## 5. NSCLC 2022 (PMC9497127)

- Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC9497127/

> "HERV quantification was performed using Telescope **with its GTF annotation file**."
> "Of the **14,174 HERV genomic loci** analyzed, **4227 (29.8%)** and 4439 (31.3%) were found expressed (average RPM > 5) in TCGA-LUAD and TCGA-LUSC."

**Implication**: application papers mostly use the official annotation; this can serve as a reference for "about 30% locus positivity under the official annotation".

---

## 6. Pan-tissue/pan-cancer retroelement database 2025 (PMC12495673)

- Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC12495673/

Uses Telescope, but **does not state which annotation** (only says it did not deviate from the authors' pipeline):

> "we did not deviate from the authors' published pipeline, which uses bowtie2"
> "11,388 retroelements were found to be actively transcribed"

Two directly reusable lessons:

1. Normalization must exclude TEs:
   > "to prevent **the very large number of lowly-expressed repetitive elements** in a RepeatMasker-derived annotation from adversely biasing median-based count normalization."
2. Gene/TE double counting is a recognized flaw:
   > "exonic reads from the GDC pipeline, which overlap with retroelements detected by the Telescope pipeline, **are double-counted**."

---

## 7. Tokuyama et al. 2018, PNAS 115(50):12565 (ERVmap paper)

- Link: https://doi.org/10.1073/pnas.1814589115

> "based on a **newly annotated database of 3,220 autonomous ERVs** that mirror a full-length provirus."

**Implication**: ERVmap's annotation system is independent of RepeatMasker; changing the annotation means it is no longer ERVmap. Keep it native.

---

## Five direct implications for our design

1. **Unifying the annotation is correct**, and converting TEtranscripts to locus level only requires swapping ids, not code (both the Telescope paper and Schwarz do this).
2. **SQuIRE can accept external annotation**: Schwarz used `.align → .out`, we use `te_locus_map.tsv → bed + copies.txt`; both routes bypass the auto-downloaded rmsk.
3. **Changing the annotation does not change the algorithm**, but it does change the assignable space; therefore we must measure "assignment rate, share of counts from short fragments" rather than reason about it.
4. **Size factors can only be computed from genes**; 720k TE loci must not be thrown into median-based normalization.
5. **Locus-level positive counts cannot be taken directly as biological conclusions** (Frontiers 2022: false positives exceed true positives); a `min_count` and a length threshold are required.

## Still unresolved / to be decided by us

- The threshold for calling a locus positive (the literature offers no directly copyable value; under the official annotation RPM>5 gives about 30%, but that is on a 14k-locus scale).
- Whether short fragments (<50 bp) should be filtered before counting, and the effect of such filtering on subfamily totals.
- The four-tool intersection (including ERVmap) has no precedent in the literature; the 2/3 consensus can only be our own choice.
