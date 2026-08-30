# Code-to-paper map

| Repository entry point | Scientific role |
|---|---|
| `analyze_dpf_parameter_sweep` | Compares DPF definitions across absorption, scattering, and source–detector separation; fits the inverse-distance model. |
| `validate_phantom_measurements` | Compares 750 nm SAP–TiO2 phantom measurements with Monte Carlo predictions. |
| `evaluate_absolute_recovery` | Assisted absolute-MBLL absorption recovery using `G ~= OD(d -> 0)`. |
| `evaluate_differential_recovery` | Non-assisted differential-MBLL absorption-change recovery using `G ~= OD(d -> 0)`. |

The `archive/` directory preserves prior and alternative versions for research
traceability. The `reviewer_analysis/` directory contains targeted checks that
are not primary paper entry points.
