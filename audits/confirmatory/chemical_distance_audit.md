# Confirmatory chemical-distance and leakage audit

- Benchmark molecules: **85**
- Endpoint experimental training rows: **1,280**
- Exact identity/connectivity firewall clean: **YES**
- Standardized parent/tautomer firewall clean: **NO**

| Source | Rows | Exact connectivity | Uncharged parent | Canonical tautomer | Shared scaffold rows | Maximum similarity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| endpoint_experimental | 1,280 | 0 | 0 | 0 | 439 | 1.000 |
| combisolv_qm | 3,963 | 0 | 0 | 2 | 874 | 1.000 |
| abraham | 8,098 | 0 | 0 | 0 | 3,000 | 1.000 |
| openff | 520 | 0 | 0 | 0 | 190 | 1.000 |
| gbn2 | 550 | 0 | 0 | 0 | 206 | 1.000 |
| molsolv_smd | 350,391 | 0 | 0 | 32 | 51,943 | 1.000 |
| confsolv | 39,878 | 0 | 22 | 22 | 7,008 | 1.000 |

The scaffold and similarity results are proximity diagnostics, not identity leakage. Every identity-equivalent match is listed in `chemical_identity_matches.csv`; nearest neighbours for every benchmark/source pair are listed in `chemical_distance_by_benchmark.csv`.
