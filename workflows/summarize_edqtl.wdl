version 1.0

# ------------------------------------------------------------------------------
# workflows/summarize_edqtl.wdl
#
# Description: compares the three edQTL approaches after the separate entry
#   points have run (Dockstore name summarize_edqtl). The end-to-end pipeline
#   runs the same task itself.
# Version:     0.1.0
# Date:        2026-10-22
# Author:      Haim Krupkin+Claude
#
# Pseudo code:
#   SummarizeEdqtlResults(tensorQTL table, caNRD-editQTL lead table, cell types)
# Inputs:  outputs of tensorqtl_edqtl.wdl, canrd_editqtl.wdl,
#          prepare_editing_phenotypes.wdl (celltypes)
# Outputs: edsite_calls_long, edsite_overview
# ------------------------------------------------------------------------------

import "tasks/summary.wdl" as summary_tasks

workflow SummarizeEdqtl {
  input {
    File? tensorqtl_all_analyses
    File? canrd_lead_site_celltype
    File celltypes
    String output_prefix
    String docker_image = "ghcr.io/aou-multiomics-analysis/caeditr_edqtl:latest"
    Int preemptible_attempts = 2
    Int max_retries = 1
  }

  call summary_tasks.SummarizeEdqtlResults {
    input:
      tensorqtl_results = tensorqtl_all_analyses,
      canrd_lead_site_celltype = canrd_lead_site_celltype,
      celltypes = celltypes,
      output_prefix = output_prefix,
      docker_image = docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  output {
    File edsite_calls_long = SummarizeEdqtlResults.edsite_calls_long
    File edsite_overview = SummarizeEdqtlResults.edsite_overview
  }
}
