# Simulation data

Place the Monte Carlo output databases here, or set the environment variable
`DPF_SIMULATION_DATA` to their location. Expected filename families are:

- `Photon_33_mua_*.mat` — parameter-sweep simulations.
- `Photon_81_mua_*.mat` — absorption-recovery simulations.
- `Photon_87_concentration_*.mat` — phantom-validation simulations.
- `Photon_58_mua_*.mat` — high-photon-count reviewer comparison.

The databases are not included in the supplied archive. See
`docs/DATA_REQUIREMENTS.md` for required variables.
