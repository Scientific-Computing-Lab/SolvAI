# Dense sentinel freeze amendment 001

**Date:** 2026-09-03 UTC  
**Original freeze commit:** `512aad3d46c25c28d590e76ebbd6cd2dc9be7186`

Before any prospective-sentinel dense window was generated, an implementation
review identified that acquisition cannot use the standard error of an unqueried
future window. Candidate observation noise is therefore fixed to the
lambda-specific median five-block SEM measured only on the four calibration
molecules. Once a window is actually acquired, its own measured SEM is used for
posterior conditioning.

This is a leakage correction and does not change the molecule panel, response
values, kernel candidates, acquisition objective, budgets, success criteria or
any other frozen analysis choice. At the time of this amendment, several
calibration windows had completed; no prospective-sentinel dense response had
been generated or read.

For consistent provenance, all run-ledger rows refer to the original scientific
freeze commit even when the runner implementation commit differs. Early
calibration rows recorded the then-current implementation commit because the
first runner version queried Git after each window; this metadata-only issue is
retained in the append-only ledger and does not affect configurations or results.
