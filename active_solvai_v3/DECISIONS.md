# Active SolvAI v3 decision log

This file is append-only after the first v3 protocol freeze. Every amendment must state
what evidence was already visible.

## 2026-09-03 — Isolate v3

- Branch: `active-solvai-v3-effort-allocation`.
- Parent: `770e13ce68ca80104557fd67a224abc5a2c44767`.
- Immutable v1 no-go: `8fb984c2eb26d016c6b81cf488f88dc667ca9cd3`.
- Immutable v2 diagnostic: `f2cc60fb73416e0e417c97c600d5abd00868031c`.
- Decision: preserve lambda-location optimization as a closed negative direction and
  test only allocation of additional sampling time over the fixed 15-point grid.
- No v3 policy result or new simulation was visible at this decision.

## 2026-09-03 — Gate 1 does not identify a SolvAI increment

- Freeze commits: `c19101b` and pre-result clarification `4608c21`.
- Generic observed diagnostics improved complementary log-difficulty MAE from
  1.547 to 1.427, but missed the frozen 10% materiality threshold.
- SolvAI conditioning worsened MAE from 1.427 to 1.475 and was
  molecule-shuffle-equivalent.
- Decision: do not fit or select a deployable v3 allocation policy from these
  data; proceed only to the predeclared quantitative power/reference gate.

## 2026-09-03 — Quantitative gate closes the campaign

- Power protocol commit: `9610813`.
- First adequate predeclared design: 256 molecules, 15 windows, two streams,
  3,000 ps per stream-window, 69,138 reserved RTX 3090 GPU-hours.
- Gate-1 provided no positive SolvAI signal, and the planning design is far
  outside the local exploratory envelope.
- Decision: launch no simulation, create no policy freeze and qualify no
  prospective sentinel. Preserve the generic diagnostic observation as a
  bounded result, not an Active SolvAI success.
