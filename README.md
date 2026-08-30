# Differential Pathlength Factor Models in CW-NIR Imaging

MATLAB analysis code accompanying the paper *Differential pathlength factor
model for near-infrared diffuse optical imaging*.

The project evaluates conventional and proposed differential pathlength factor
(DPF) definitions for continuous-wave near-infrared diffuse optical imaging.
Its principal model is

```text
DPF(d) = A + B/d + C log(d)/d,
```

where `d` is the source–detector separation.

## Repository organization

- `analysis/parameter_sweep/` — optical-property sweep and DPF comparison.
- `analysis/absorption_recovery/` — absolute and differential MBLL recovery.
- `experiments/` — phantom measurement validation.
- `reviewer_analysis/` — diffusion-validity and simulation checks.
- `archive/` — clustered development and alternative-model versions.
- `config/` — project-relative data configuration.
- `data/` — documented placeholders for non-included inputs.
- `docs/` — manuscript and reproducibility documentation.

## Getting started

1. Clone the repository and start MATLAB in its root directory.
2. Run `startup`.
3. Add the required databases described in `docs/DATA_REQUIREMENTS.md`, or set
   `DPF_SIMULATION_DATA` to an external data directory.
4. Run one of the four publication entry points listed in `docs/CODE_MAP.md`.

Example:

```matlab
startup
analyze_dpf_parameter_sweep
```

## Reproducibility status

The supplied source archive did not contain the Monte Carlo databases, raw
experimental tables, simulation-generation drivers, or three instrument
characterization functions. Consequently, the present package documents and
organizes the complete supplied analysis code but cannot yet reproduce every
paper figure from a fresh clone. Missing inputs are listed explicitly in
`docs/DATA_REQUIREMENTS.md`.

## Version history

Only the recommended paper-facing workflows are placed under `analysis/` and
`experiments/`. Near-duplicate development versions are retained in clustered
`archive/` folders so that scientific history is not lost or confused with the
recommended workflow.

## License

All rights reserved until the authors select a public software license. The
repository should remain private until the data, dependency, and licensing
items above are resolved.
