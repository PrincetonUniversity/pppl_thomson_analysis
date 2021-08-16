# PPPL Thomson Analysis

A collection of Python functions written to assist with data collection and
analysis for the real-time multi-point Thomson scattering evaulation system
(rt-MPTS) originally developed for NSTX-U.

**Developers:**    K. C. Hammond, F. M. Laggner, and R. Rozenblat   \
**Maintained by:** K. C. Hammond                                    \
**Contact:**       khammond@pppl.gov                                \
**Last updated:**  2021-07-12


## Software requirements

- Python 3+
- [MDSplus (Python module)](https://mdsplus.org)

## Modules and scripts

### generate_nstxu_rtmpts_tree.py

Contains functions for generating an MDSplus model tree file for storing
calibration and experimental data from the NSTX-U real-time multi-point 
Thomson scattering (rt-MPTS) system.

The tree structure will likely be different for other systems, but this module
may be instructive as a template for using the Python MDSplus module for 
constructing model trees.

To create the sample tree files with this script, it is first necessary
to define an environment variable indicating where the tree files should be
stored. The environment variable should have the name `[tree_name]_path`, where
`[tree_name]` is the name of the tree.

To create files for a tree with the name `nstxu_example` and to store the 
underlying files in the working directory, the following commands should be
entered:

    >> export nstxu_example_path=$(pwd)
    >> python generate_nstxu_rtmpts_tree.py nstxu_example

### polychromator.py

Defines the Polychromator object class. Instances of this class store 
experimental or idealized polychromator transmissivity spectra and can 
generate synthetic Thomson scattering signals based on the spectra.

A tutorial with example usage cases is available 
[here](tutorials/tutorial_polychromator.md).

### pulse_eval.py

Contains functions for extracting signal strength parameters from time-resolved
polychromator signals. The signal strength parameters are used as input for
the spectral fits.

- **processPulse():** Calculates a scalar strength parameter for a raw,
  time-resolved polychromator signal according to a user specified method.
  Currently supports peak detection, trapezoidal integration, and Gaussian
  curve-fitting.

  The keyword arguments (kwargs) generally correspond to properties specific
  to particular diagnostic system. They will probably not need to be changed
  from laser pulse to laser pulse or from shot to shot, but they must be
  chosen carefully any time the function is applied to a new system.
- **fitPulse():** Uses the Levenberg-Marquardt least-squares fitting algorithm
  to find best-fit parameters for a user-specified curve model to a dataset.
- **singleGaussian():** Computes the Gaussian function; for use in `fitPulse()`.
- **doubleGaussian():** Computes a double Gaussian function (the sum of two
  Gaussians); for use in `fitPulse()`.
- **jacobian_singleGaussian():** Computes the Jacobian matrix associated with
  `singleGaussian()`; for use in `fitPulse()`.
- **jacobian_doubleGaussian():** Computes the Jacobian matrix associated with
  `doubleGaussian()`; for use in `fitPulse()`.

### spectral_fitting.py

Contains functions for fitting polychromator data to the Thomson scattering
spectrum.

- **fitSpectrum():** Performs a fit of polychromator data to polychromator
  spectral calibration data with a Levenberg-Marquardt least-squares algorithm.
- **lookup():** Calculates expected polychromator signals and their derivatives
  for a given electron temperature and scaling factor by interpolating the
  calibration spectra.

### tree_queries.py

Contains miscellaneous functions that extract data from an MDSplus tree.

At present, the functions contained within this file are specifically tailored 
to the implementation of the real-time system on the LHD experiment. As such, 
they may not work for other implementations in which the tree structure may not 
be exactly the same.

- **getVoltageData():** Obtains raw polychromator signals in volts for a given 
  shot, spatial volume, and polychromator channel.

- **getRawData():** Obtains raw polychromator signals in terms of the unscaled 
  digitizer signals for a given shot, spatial volume, and polychromator channel.
  This query is independent of the digitizer calibration data.

## Relevant publications

- F. M. Laggner, A. Diallo, B. P. LeBlanc, R. Rozenblat, G. Tchilinguirian,
  E. Kolemen, and the NSTX-U Team, A scalable real-time framework for
  Thomson scattering analysis: Application to NSTX-U, 
  [*Review of Scientific Instruments* **90**, 043501 (2019)](https:/doi.org/10.1063/1.5088248)
- R. Rozenblat, E. Kolemen, F. M. Laggner, C. Freeman, G. Tchilinguirian, 
  P. Sichta, and G. Zimmer, Development of real-time software for Thomson
  scattering analysis at NSTX-U, [*Fusion Science and Technology* **75**, 835 (2019)](https://doi.org/10.1080/15361055.2019.1658037)
- K. C. Hammond, F. M. Laggner, A. Diallo, S. Doskoczynski, C. Freeman, 
  H. Funaba, D. A. Gates, R. Rozenblat, G. Tchilinguirian, Z. Xing, I. Yamada, 
  R. Yasuhara, G. Zimmer, and E. Kolemen, Initial operation and data processing
  on a system for real-time evaluation of Thomson scattering signals on the
  Large Helical Device, [*Review of Scientific Instruments* **92**, 063523 (2021)](https://doi.org/10.1063/5.0041507)


