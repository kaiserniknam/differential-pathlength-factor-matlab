# Data and dependency requirements

## Monte Carlo database variables

The parameter-sweep files (`Photon_33`) use `x`, `y`, `z`, `d`, `s`, `w`,
`c`, `a`, and `no_of_photons`. The recovery files (`Photon_81`) use photon
exit positions, pathlengths, weights, and photon count. The phantom files
(`Photon_87`) use `x_in`, `y_in`, `z_in`, `x_ot`, `y_ot`, `z_ot`, `s`, `w`,
and `no_of_photons`.

The full data arrays may be too large for ordinary Git. Publish them through a
versioned research-data archive and place checksums plus the DOI in this file.

## Missing instrument functions

Experimental conversion calls `do_FDS100`, `do_LED750L`, and `do_LED850LN`.
These functions were not in the supplied archive. Add documented, licensed
implementations under `utilities/characteristics/` before claiming that the
experimental workflow is reproducible.

## MATLAB products

- MATLAB
- Statistics and Machine Learning Toolbox (`fitlm`)

Some scripts use functions such as `nanmean` and `nanstd`; compatibility with
the intended MATLAB release should be verified before release.
