# ProteoBench PYE/plasma scoring

**Status:** open follow-up
**Extracted:** 2026-07-24 from the completed HYE implementation tracker

## Goal

Add PYE/plasma scoring to APB using the existing matrix-native ProteoBench
intermediate, storage, CLI, provenance, and APB Studio stage contracts.

## Requirements

- Define the PYE/plasma-specific expected ratios and score metrics without
  duplicating the HYE pipeline.
- Add golden intermediate and score fixtures from a pinned ProteoBench revision.
- Advertise only module/vendor combinations whose full golden parity passes.
- Preserve the existing score JSON, `varm`, `uns`, and Studio display contracts.
- Keep vendor parsing and value completion in APB conversion rules rather than
  adding vendor switches to the scorer.

## Done when

- PYE/plasma golden intermediate values and score JSON match ProteoBench.
- APB and APB Studio full quality gates pass.
- README, architecture, and supported-module documentation list only verified
  combinations.
