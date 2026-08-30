# Experimental data

The validation workflow expects these tab-delimited detector-grid recordings:

- `OD_grid_2025-12-29-11-40-40-750nm-OD-1gL.txt`
- `OD_grid_2025-12-30-12-15-15-750nm-OD-0_25gL.txt`

They were referenced by the source code but were not present in the supplied
archive. Each table is expected to contain time, detector voltage, grid row,
grid column, and recording-state fields in the column positions documented in
`validate_phantom_measurements.m`.
