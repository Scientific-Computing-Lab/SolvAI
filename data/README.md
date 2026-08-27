# Data layout

- `benchmark/`: the 85-solute reference table, including source provenance,
  physical reference values and frozen fold assignments.
- `processed/`: redistributable benchmark-disjoint teacher tables required for
  audit or documented full reconstruction.
- `manifests/`: source versions, licences, checksums, counts and filtering
  decisions.

The release does not redistribute a third-party source whose archive lacks
clear standalone redistribution terms. Such sources are represented by a
canonical URL, DOI and exact hash in the manifests. All headline predictions
and metrics are included independently of those source archives.

This policy applies to the processed CombiSolv-QM table and to the two merged
public-hydration tables containing CombiSolv-Exp records. Their exact frozen
metadata are in `manifests/training_source_manifest.json` and
`manifests/endpoint_label_manifest.json`, respectively.

Every supervised external source used by the final eligible model was screened
against benchmark canonical SMILES, full InChIKey and connectivity block before
training. See `../audits/leakage_audit.md`.
