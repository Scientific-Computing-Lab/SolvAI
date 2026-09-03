# Run and decision ledgers

`ledger.jsonl` is append-only and records every analysis, simulation attempt, failure,
retry, input/output hash, device and cost. `allocation_decisions.parquet` is generated
only after a policy is frozen and records every candidate utility, chosen action and
stop decision at each replay step.

