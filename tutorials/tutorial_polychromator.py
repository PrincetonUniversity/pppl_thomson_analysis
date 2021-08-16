"""
tutorial_polychromator.py

Demonstrates the usage of the Polychromator class as well as the pulse
and spectral fitting functions.
"""

import numpy as np
import matplotlib.pylab as pl
import pulse_eval as pe
import spectral_fitting as sf
import polychromator as pc

# Colors for plotting
col_line = [(0.8, 0, 0), (1.0, 0.5, 0), (0, 0.6, 0), (0, 0, 0.8), (0.5, 0, 1.0)]
col_err  = [(1, .5, .5), (1, .75, .5), (.5, 1, .5), (.5, .5, 1), (.75, .5, 1)]
pl.rcParams['font.size']=8

if __name__ == '__main__':

    # Create the instance of the Polychromator class with a defined laser
    # wavelength and scattering angle associated with the polychromator's 
    # scattering volume
    wvl_laser = 1064.
    scat_ang = 167.*np.pi/180.
    pcTest = pc.Polychromator(wvl_laser=wvl_laser, scat_ang=scat_ang)

    # Define idealized transmissivity spectra for a 5-channel polychromator
    n_channels = 5
    transm_levels = [0.9, 1.1, 1.3, 1.2, 1.0]
    cutoff_low    = [1030., 1010., 950., 840., 690.]
    cutoff_high   = [1050., 1030., 1010., 950., 840.]
    transm_errs   = [0.05, 0.05, 0.05, 0.05, 0.05]
    wvl = np.linspace(600, 1100, 501) # units must be consistent with wvl_laser
    pcTest.set_transmissivity_ideal(transm_levels, cutoff_low, cutoff_high, \
                                    wvl, errs=transm_errs)

    # Extract the transmissivity spectra for plotting
    transm = [[]]*n_channels
    dtransm = [[]]*n_channels
    wvl = [[]]*n_channels
    for i in range(n_channels):
        transm[i], dtransm[i], wvl[i] = pcTest.transmissivity(i)

    # Plot and save a figure with the transmissivity spectra
    f_tr = pl.figure(figsize=(4,2.6))
    ax_tr = f_tr.add_subplot(1,1,1)
    ax_tr.set_position([0.15, 0.15, 0.8, 0.8])
    h_tr = []
    for i in range(n_channels):
        ax_tr.fill_between(wvl[i], transm[i]+dtransm[i], transm[i]-dtransm[i],\
                           color=col_err[i])
        h_tr.append(ax_tr.plot(wvl[i], transm[i], color=col_line[i], \
                               label='Channel %d'%(i)))
    ax_tr.set_xlabel('Wavelength (nm)')
    ax_tr.set_ylabel('Transmissivity (A.U.)')
    pl.legend([h_tr_plot[0].get_label() for h_tr_plot in h_tr])
    f_tr.savefig('pc_tutorial_transmissivity.png', dpi=150)
    pl.close(f_tr)

    # Output the lookup tables for expected signal versus electron temperature
    spec, dspec, spec_te = pcTest.expected_relative_signal_spectra()

    # Plot and save a figure with the expected signal spectra
    f_es = pl.figure(figsize=(4,2.6))
    ax_es = f_es.add_subplot(1,1,1)
    ax_es.set_position([0.15, 0.15, 0.8, 0.8])
    h_es = []
    for i in range(n_channels):
        ax_es.fill_between(spec_te, spec[:,i]+dspec[:,i], spec[:,i]-dspec[:,i],\
                           color=col_err[i])
        h_es.append(ax_es.plot(spec_te, spec[:,i], color=col_line[i], \
                               label='Channel %d'%(i)))
    ax_es.set_xlabel('Electron temperature (keV)')
    ax_es.set_ylabel('Transmissivity (A.U.)')
    pl.legend([h_es_plot[0].get_label() for h_es_plot in h_es])
    f_es.savefig('pc_tutorial_exp_signals.png', dpi=150)
    pl.close(f_es)

    # Generate synthetic signals for a temperature of 3.0 keV
    te_test = 3.0
    am_test = 1.0
    sigs, times = pcTest.synthetic_signals(te_test, am_test)

    # Plot and save a figure with the synthetic signals
    f_sg = pl.figure(figsize=(4,2.6))
    ax_sg = f_sg.add_subplot(1,1,1)
    ax_sg.set_position([0.15, 0.15, 0.8, 0.8])
    h_sg = []
    for i in range(n_channels):
        h_sg.append(ax_sg.plot(1.e9*times, sigs[i], color=col_line[i], \
                               label='Channel %d'%(i)))
    ax_sg.set_xlabel('Time (ns)')
    ax_sg.set_ylabel('Synthetic signal (A.U.)')
    pl.legend([h_sg_plot[0].get_label() for h_sg_plot in h_sg])
    f_sg.savefig('pc_tutorial_synth_signals.png', dpi=150)
    pl.close(f_sg)

    # Calculate signal strengths using the peak-detect method
    print('Calculating signal strengths using the peak-detect method')
    print('(strengths will be equal to the pulse amplitudes)')
    sigValsPk = [[]]*n_channels
    dSigValsPk = [[]]*n_channels
    info = [[]]*n_channels
    for i in range(n_channels):
        sigValsPk[i], dSigValsPk[i], info[i] = \
            pe.processPulse(sigs[i], method='peak')

    # Fit the signal strengths to the simulated signals in the look-up table
    te0 = 2.0
    s0  = 2.0
    sigOK = [True]*n_channels
    sFit, teFit, dsFit, dteFit, chi2, _, _, _ = \
        sf.fitSpectrum(sigOK, sigValsPk, dSigValsPk, s0, te0, \
                       spec_te, spec, dspec, verbose=True)

    # Calculate signal strengths using the Gaussian fitting method
    print('Calculating signal strengths using the curve-fitting method')
    print('(strengths will be equal to the integral of the pulse)')
    sigValsFit = [[]]*n_channels
    dSigValsFit = [[]]*n_channels
    info = [[]]*n_channels
    for i in range(n_channels):
        sigValsFit[i], dSigValsFit[i], info[i] = \
            pe.processPulse(sigs[i], method='fit')

    # Fit the signal strengths to the simulated signals in the look-up table
    te0 = 2.0
    s0  = 2.0
    sigOK = [True]*n_channels
    sFit, teFit, dsFit, dteFit, chi2, _, _, _ = \
        sf.fitSpectrum(sigOK, sigValsFit, dSigValsFit, s0, te0, \
                       spec_te, spec, dspec, verbose=True)

