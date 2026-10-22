# Data dictionary

All tables are tab-separated text (plain or gzip/bgzip). Sample IDs are text
(AoU `research_id`s look numeric; they are never converted to numbers) and
must be identical across files. Chromosome names must match between the BEDs,
the GTF and the PLINK `.pvar` (AoU: `chr1` ... `chr22`).

## Inputs

| Input | Format | Notes |
| --- | --- | --- |
| `editing_ratios` | BED matrix: `#chr  start  end  site_id  <sample>...` | Bulk editing ratio G/(A+G) in [0, 1]; `NA` = no reads. One row per editing site. `site_id` is the caEditR id `chrom:pos:strand` (no `chr`, 1-based, e.g. `1:12345:+`); `start = pos - 1`, `end = pos`. Plain or gzipped, any order (the pipeline sorts and indexes it). |
| `editing_coverage` | Same layout and same sites as `editing_ratios` | Read depth (A+G reads) per site and sample; 0 = no reads. |
| `expression` | BED matrix: `#chr  start  end  gene_id  <sample>...` | Linear CPM (the same file given to CellTypeDeconvolution). Used for theta (host-gene expression per cell type). |
| `gtf` | GTF / GTF.GZ | `gene` records with `gene_id` (matching `expression`), `gene_name`, `gene_type`. Sites are mapped to host genes by gene-body overlap. |
| `proportions` | TSV: `sample_id` + one column per cell group | CellTypeDeconvolution output `proportions_tca_weights.tsv` (retained major groups, rows sum to 1), e.g. `Neutrophils`, `CD4 T cells`, `Monocyte/myeloid`. |
| `covariates` | TSV: `sample_id` + numeric columns | QTL covariates (genotype PCs, sex, age; categorical variables already encoded). Same layout as prepare_QTL `AdditionalCovariates`. Samples with missing values are dropped. |
| `plink_pgen/pvar/psam` | PLINK2 | From prepare_QTL `prepare_VCF` / `convertVCF2Plink` (variant IDs `chr1:123_A_G`). Effects are per ALT allele. |
| `sample_list` | Optional, one ID per line | Optional header `research_id` / `sample_id` / `ID`. |

The analysis cohort is the intersection of all inputs, in the column order of
`editing_ratios` (`sample_alignment_report.tsv` lists which input lacks which
sample).

## Main settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `min_coverage`, `min_covered_samples` | 20, 60 | A site is analysed when >= 60 samples have >= 20 reads (Li et al. 2022). `min_covered_fraction` adds a fraction rule. |
| `theta_floor` | 0.001 | NNLS floor: theta at the floor = host gene not expressed in that cell type (caEditR gating). |
| `sites_per_chunk` | 500 | Sites per scatter shard (100k sites -> 200 shards). |
| `n_editing_pcs` | 15 | Editing PCs of each phenotype matrix added to its tensorQTL covariates. |
| `max_missing_fraction` | 0.5 | Sites with more missing values in a phenotype matrix are dropped (caNRD-edit leaves non-identifiable cell types `NA`). |
| `cis_window` | 1000000 | +/- bp around each site (inclusive), the same for tensorQTL (nominal and permutations) and caNRD-editQTL; 1 Mb as in the GTEx edQTL analysis. |
| `canrd_n_editing_pcs` | 0 | Bulk editing PCs added to the caNRD-editQTL covariates. |
| `canrd_maf`, `canrd_min_minor_allele_samples` | 0.01, 10 | Variant filters of caNRD-editQTL. |
| `canrd_min_mean_phi` | 0.10 | A cell type is tested at a site only when its mean RNA share is >= 10 %. |

## Outputs (per stage)

| Output | Content |
| --- | --- |
| `site_to_gene.tsv`, `site_theta.tsv`, `gene_theta.tsv` | Host gene per site; per-cell-type theta (sites without a host gene get theta = 1, i.e. RNA share = proportion). |
| `<prefix>.site_qc.tsv.gz` | Per site: covered samples, mean / sd ratio, pass / reason. |
| `<prefix>.bulk.all.bed.gz` | QC-passing bulk ratios (input of the standard bulk edQTL scan). |
| `<prefix>.canrd_edit.<celltype>.bed.gz` | caNRD-edit per-sample cell-type editing (ratio scale; `NA` = not identifiable). |
| `<prefix>.tca.<celltype>.bed.gz` | TCA per-sample cell-type editing (ratio scale). |
| `<prefix>.canrd_diagnostics.tsv.gz` | caNRD-edit per-site status, condition number, excluded cell types. |
| `<label>.INT.bed.gz` + `<label>.tensorqtl_covariates.tsv` | tensorQTL phenotype (rank inverse-normal per site) and covariates (`ID` + samples; supplied covariates + editing PCs). |
| `<label>.cis_qtl_pairs.<chr>.parquet` | tensorQTL cis nominal (default): every site x variant pair (`pval_nominal`, `slope`, `slope_se`, `af`, `ma_samples`). |
| `<label>.cis_qtl.txt.gz` | Only with `run_permutations = true`: tensorQTL permutation result (lead variant, `pval_perm`, `pval_beta`). tensorQTL's own q-values are not used. |
| `<prefix>.tensorqtl_cis.all_analyses.tsv.gz` | All tensorQTL sets in one table, one row per site and set: lead variant, `n_variants_tested`, raw `pval_nominal`, `slope`, `slope_se` (`analysis` = bulk / canrd_edit / tca, `celltype_slug`). |
| `<chunk>.canrd_editqtl.all_pairs.tsv.gz` | caNRD-editQTL: every site x variant x cell type (`beta`, `se`, `p`, `p_site`, `status`, `vif`, `weakly_identifiable`, `refined`). |
| `<prefix>.canrd_editqtl.site_level.tsv` | Per site: lead variant, exact `p_site_lead` (K-df LRT), `n_variants_tested`. |
| `<prefix>.canrd_editqtl.lead_site_celltype.tsv` | Per site x cell type at the lead variant: exact LRT `p`, `beta`, `se`, VIF. |
| `<prefix>.edsite_calls_long.tsv.gz`, `<prefix>.edsite_overview.tsv` | The three approaches side by side. |

The q-value / "edsite" columns of the aggregated tables are a convenience
(lead p x number of cis variants of the site, then Storey q across sites); every
table keeps the raw p-values for your own multiple-testing correction.
