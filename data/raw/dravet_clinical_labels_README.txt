DeepGene Dravet clinical-label schema

Required columns:
  variant_key       SCN1A:p.Arg1648His or SCN1A:p.R1648H
  dravet_label      1 = clinically confirmed Dravet syndrome; 0 = explicitly non-Dravet SCN1A-related phenotype
  label_source      publication, registry, or curated clinical cohort
  source_accession  PMID, database accession, or de-identified cohort identifier
  cohort_id         patient/cohort grouping used to prevent related samples crossing splits
  reviewer_status   independently reviewed / unresolved

Do not infer dravet_label from ClinVar pathogenicity, LOF/GOF, phenotype text alone,
or the prediction output of another model. Unknown and unresolved cases must be
excluded rather than assigned 0.
