# Terra WDL maintenance

These pipelines target Terra's managed Cromwell environment (WDL 1.0). The
rules are the AoU-Multiomics-Analysis rules (see prepare_QTL
`docs/terra-file-paths.md`):

- Do not create files at workflow scope (`python3 tools/check_wdl_file_scope.py workflows`).
- Pass arguments as typed WDL inputs and safely quoted named CLI arguments; no
  JSON argument wrappers. Arrays of files: `write_lines()` directly in the
  command placeholder, copied to a task-local file.
- Keep files File-typed until localization; never wrap a generated File in `sub`.
- Link BED + `.tbi` and the PLINK trio into one task-local directory
  (Cromwell can localize them into different directories).
- Every command block logs `stage=`, `start_time=`, `completion_time=`,
  `dimensions=`, `outputs=` (`python3 tools/check_wdl_logging.py workflows`).
- Validate with `miniwdl check` and Cromwell `womtool validate` (CI does both).
- Tables are TSV/CSV; avoid JSON for parameters.
- State separately whether syntax checks, local miniwdl/Docker runs and a
  Terra run have passed. Local runs do not establish Terra success.
