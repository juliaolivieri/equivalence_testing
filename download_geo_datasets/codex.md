Task: Search GEO for mouse bulk RNA-seq datasets suitable for later differential-expression and equivalence-analysis work. This is a dataset discovery and screening task only.

Scope:

* Search the full GEO window covering the last 2 years.
* Output one row per GSE series.
* Include eligible, candidate, and exclude rows.
* Do not download supplementary files unless absolutely necessary to confirm a candidate.

Inclusion criteria:

* Organism must be Mus musculus.
* Study must be bulk RNA-seq.
* Exclude single-cell RNA-seq, single-nucleus RNA-seq, spatial transcriptomics, microarrays, and other non-bulk assays unless the series clearly contains standard bulk RNA-seq gene-expression data.
* GEO series should be released within the last 2 years; if release date is unavailable, use submission date.
* There must be one clear baseline/control condition and at least one experimental/perturbation condition.
* Accept control labels such as wild type, WT, control, untreated, vehicle, sham, uninfected, or other clearly baseline conditions.
* At least one usable comparison must have at least 2 biological replicates per condition.
* A supplementary file should contain a gene-expression matrix that is either clearly raw counts or likely raw counts.
* There must be enough metadata to assign samples to conditions without guesswork.

Exclusion criteria:

* Non-mouse organism
* Not bulk RNA-seq
* No clear control-versus-experimental comparison
* Fewer than 2 replicates per condition for all usable comparisons
* No gene-level expression matrix available
* Only normalized values appear available, with no raw-count matrix
* Metadata too ambiguous to assign samples to groups
* Study design too complex to identify even one clean pairwise comparison

Important behavior rules:

* Be conservative. If uncertain, do not guess; mark the row as candidate or unclear.
* Save the CSV incrementally as you go, at least every 10 GSEs examined.
* Keep a short audit log or notes field explaining why each row was marked eligible, candidate, or exclude.
* Do not silently skip ambiguous studies.
* Stop only when the full last-2-years search space has been examined or no further relevant GSEs can be found through reasonable GEO search pagination.

For each GSE, produce a CSV row with these columns:

* gse_accession
* title
* geo_url
* organism
* assay_type
* release_date
* total_samples
* condition_names
* num_conditions
* usable_comparison_exists
* num_usable_comparisons
* comparison_description
* replicates_by_condition
* min_replicates_per_condition
* supplementary_file_names
* matrix_type_assessment
* metadata_sufficiency
* inclusion_status
* confidence
* notes
* reason_for_exclusion_or_uncertainty

Allowed values:

* matrix_type_assessment: raw_counts_confirmed, raw_counts_likely, not_raw_counts, unclear
* metadata_sufficiency: clear, partially_clear, ambiguous
* inclusion_status: eligible, candidate, exclude
* confidence: high, medium, low

At the end:

* Return the final CSV file
* Return a short summary with counts of eligible, candidate, and exclude rows
* Report any recurring ambiguity patterns that may require a second-pass manual review


Previous prompts, kept for historicity; do not implement these.
Below is a longer prompt that I eventually want done, but first I want you to try a pilot. So you can read the full spec, but then just try this first pilot before proceeding with the whole thing: Search GEO for up to 15 mouse bulk RNA-seq GSE series released in the last 2 years that may fit the screening criteria. Do not download supplementary files yet unless absolutely necessary. Output one CSV row per GSE, including eligible, candidate, and exclude rows, and save partial progress as you go. Prioritize simple control-versus-experimental designs and clear metadata.

Full prompt (do not do this yet):
Task: Search GEO for mouse bulk RNA-seq datasets suitable for later differential-expression and equivalence-analysis work. Do not download supplementary files yet. This task is only for dataset discovery and screening.

Inclusion criteria:
- Organism must be Mus musculus.
- Study must be bulk RNA-seq.
- Exclude single-cell RNA-seq, single-nucleus RNA-seq, spatial transcriptomics, microarrays, and other non-bulk assays unless the series clearly contains standard bulk RNA-seq gene-expression data.
- GEO series must have been deposited or released within the last 2 years.
- There must be one clear baseline/control condition and at least one experimental/perturbation condition.
- Accept control labels such as wild type, WT, control, untreated, vehicle, sham, uninfected, or other clearly baseline conditions.
- At least one usable comparison must have at least 2 biological replicates per condition.
- A supplementary file should contain a gene-expression matrix that is either clearly raw counts or likely raw counts.
- There must be enough metadata to assign samples to conditions without guesswork.

Exclusion criteria:
- Non-mouse organism
- Not bulk RNA-seq
- No clear control-versus-experimental comparison
- Fewer than 2 replicates per condition for all usable comparisons
- No gene-level expression matrix available
- Only normalized values appear available, with no raw-count matrix
- Metadata too ambiguous to assign samples to groups
- Study design too complex to identify even one clean pairwise comparison

Important behavior rules:
- Do not download or parse supplementary files yet unless absolutely necessary to confirm a candidate.
- Prefer information available from the GEO series page, GSM metadata, and supplementary-file listings.
- Be conservative. If uncertain, do not guess; mark the row as candidate or unclear.
- It is acceptable to include excluded rows with partial information so the output shows what was checked.
- Output one row per GSE series, not one row per comparison.
- Include the number of usable comparisons when possible.

For each GSE, produce a CSV row with these columns:
- gse_accession
- title
- geo_url
- organism
- assay_type
- release_date
- total_samples
- condition_names
- num_conditions
- usable_comparison_exists
- num_usable_comparisons
- comparison_description
- replicates_by_condition
- min_replicates_per_condition
- supplementary_file_names
- matrix_type_assessment
- metadata_sufficiency
- inclusion_status
- confidence
- notes
- reason_for_exclusion_or_uncertainty

Allowed values:
- matrix_type_assessment: raw_counts_confirmed, raw_counts_likely, not_raw_counts, unclear
- metadata_sufficiency: clear, partially_clear, ambiguous
- inclusion_status: eligible, candidate, exclude
- confidence: high, medium, low

Prioritization:
- Prefer simple designs with one control and one experimental group.
- Prefer datasets with clearer sample annotation.
- Prefer datasets with explicit raw-count matrices.
- Prefer stronger replication.

Return the final output as a CSV file.
