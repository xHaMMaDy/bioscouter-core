# Usability and Ethics Protocol

This protocol is an optional post-submission validation plan for BioScouter. It
is included to make any future usability claims auditable and to avoid implying
that unpublished user-study data already exist.

## Scope

The protocol evaluates whether bioinformatics users can complete common dataset
discovery tasks using the hosted BioScouter interface. It does not collect new
biological samples, clinical records, participant-level health data, or private
repository data.

## Participants

- Target sample: 5-10 bioinformatics, genomics, or computational biology users.
- Inclusion criteria: prior experience searching at least one public omics
  repository.
- Exclusion criteria: direct authorship on the BioScouter manuscript or direct
  involvement in BioScouter development.
- Recruitment: academic collaborators, students, or research staff invited
  without coercion.

## Consent and Privacy

Participants should receive a short information sheet explaining the study
purpose, voluntary nature, approximate duration, and data collected. Collect only
minimal study data: task completion, time on task, usability ratings, and
optional free-text comments. Do not collect passwords, OAuth credentials, private
dataset identifiers, patient data, IP addresses for analysis, or identifiable
biomedical information.

If institutional policy classifies this as human-subjects research, obtain the
required ethics/IRB determination before starting. If classified as service
evaluation or exempt usability testing, retain the determination note with the
study materials.

## Tasks

1. Run a keyword search for `TMT proteomics`.
2. Open one result and identify source, accession, organism, sample metadata,
   and metadata-readiness indicators.
3. Export search results as CSV.
4. Run a natural-language research question: `What genes are dysregulated in
   early-stage Alzheimer's disease?`
5. Use filters or sorting to identify a candidate dataset for reuse.
6. Locate citation or source-link information for one dataset.

## Measures

- Task completion: complete, partial, or not complete.
- Time on task: minutes and seconds.
- Error count: wrong route, failed query, wrong export, or incorrect accession.
- Usability rating: 1-5 Likert score after each task.
- Overall System Usability Scale (SUS), optional.
- Free-text comments: optional, de-identified.

## Analysis Plan

Report only aggregate descriptive statistics unless the sample size is expanded
and a statistical analysis plan is pre-registered. Do not claim clinical,
diagnostic, or biological validation from usability results. Report limitations,
including participant expertise, small sample size, and use of a live system
whose upstream repository responses can change.

## Data Retention

Store de-identified task sheets and summary outputs in the manuscript
reproducibility archive only if participants consent and institutional policy
allows sharing. Otherwise, report the availability restriction and keep raw
participant sheets private.
