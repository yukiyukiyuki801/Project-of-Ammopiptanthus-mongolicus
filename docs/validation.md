# Release validation

All1,615 stored assets and1,667 reconstructed input/result files passed SHA256 verification. Required inputs for all26 latest figure sources exist after restoration. All52 delivered figure PDF/PNG hashes match the latest manuscript package, and the workbook has28 tables. All24 Python and44 R sources passed syntax checks. Shell scripts passed bash syntax checks. The pinned ecospat/hierfstat fetches retrieved exact verified source bytes.

The final GAM fitter/postprocessor/input files retain the validated module's hashes, and its portable input-preparation check regenerated the production380-row input identically. Mapping configuration checks exercise217-sample identity and ten-shard498-contig coverage without processing study reads. Optional figure-QA control tests verify explicit NOT_RUN status when external helpers are absent.

`provenance/released_code_index.json` is the current code-hash index. Source-selection and source-copy manifests document earlier staging or original execution hashes; later portability edits are reflected in the current index. Original executed-source copies in mapping provenance are lineage evidence, not the portable CLI entry points. Source providers retain attribution.

Static checks found no raw FASTQs/BAMs, manuscript/reviewer documents, installed environments, private keys or recognizable credentials in the release. Code author and committer are configured as Pei-Wei Sun using the existing GitHub-linked noreply address. Public protein sequences can contain short patterns resembling token prefixes; these were checked as biological sequence data, not credentials.

No statistical fits, EEMS MCMC, mapping, calling or figure rendering were repeated while assembling this release. This is not a claim that every analysis has been rerun successfully in an independent fresh environment. See RELEASE_VALIDATION.json, INDEPENDENT_STATIC_AUDIT.json and the per-module validation records for the actual scope of checks.
