"""
polychromator.py

Defines a class for storing polychromator calibration data, calculating 
expected signal spectra, and generating synthetic signals.
"""

import numpy as np

# Physical constants
me = 9.1093837015e-31
qe = 1.602176634e-19
c  = 2.99792458e8

class Polychromator(object):
    """
    Encapsulates a polychromator, containing data for the transmissivity of
    each channel and outputting relative signal spectra on request.
    """

    def __init__(self, wvl_laser=[], scat_ang=[], spec_te=[]):
        """
        Creates an instance of the Polychromator class. Input parameters (laser
        wavelength and scattering angle) are optional, but they must be set 
        prior to querying the class instance for spectral data or synthetic 
        signals.

        Parameters
        ----------
            wvl_laser: double (optional)
                Wavelength of Thomson laser. Units must be consistent with
                the those associated with the transmissivity spectra.
            scat_ang: double (optional)
                Scattering angle of the light received by the polychromator
                (radians)
            spec_te: double array (optional)
                Array of electron temperatures (keV) at which to evaluate
                expected signal strength when calculating spectral curves
        """

        self.n_channels = 0
        self.wvl = []
        self.transm = []
        self.dtransm = []
        self.spec_te = []
        self.spec = []
        self.dspec = []
        self.scat_ang = scat_ang
        self.wvl_laser = wvl_laser
        self.spec_calculated = False

        if not spec_te:
            self.spec_te = np.linspace(0.01, 10., 1000)
        else:
            self.spec_te = spec_te

    def set_transmissivity_ideal(self, levels, cutoff_low, cutoff_high, \
                                 wvl, errs=[]):
        """
        Set the transmissivity spectra of each channel as idealized filters
        with a constant level within the pass band and zero outside.

        Parameters
        ----------
            levels: double array
                Transmissivity level (arb. units) within the pass-band for 
                each channel
            cutoff_low: double array
                Lower cutoff wavelength for each channel. Units must be 
                consistent with those of the specified laser wavelength.
            cutoff_high: double array
                Upper cutoff wavelength for each channel. Units must be 
                consistent with those of the specified laser wavelength.
            wvl: double array 
                Array of wavelengths at which to store the transmissivity
                spectra of each channel for the purpose of calculating the 
                expected relative signals. Units must be consistent with those 
                of the specified laser wavelength.
            errs: double array (optional)
                Uncertainties of the values in `levels` (zero by default)
        """

        # Check inputs
        n_channels = len(levels)
        if len(cutoff_low) != n_channels or len(cutoff_low) != n_channels or \
           len(cutoff_high) != n_channels:
            raise ValueError('Input arrays `levels`, `cutoff_low`, and ' \
                             + '`cutoff_high` must have the same number of ' \
                             + 'elements')
        if not errs:
            errs = np.zeros(np.shape(levels))
        else:
            if len(errs) != n_channels:
                raise ValueError('Input `errs` must agree in length ' \
                                 + 'with `levels`')

        self.n_channels = n_channels

        self.wvl = []
        for i in range(self.n_channels):
            self.wvl.append(np.array(wvl))

        # Construct transmissivity spectral arrays and populate self.transm
        self.transm = []
        self.dtransm = []
        for i in range(self.n_channels):
            inds_passband = np.logical_and(self.wvl[i] >= cutoff_low[i], \
                                           self.wvl[i] <= cutoff_high[i])
            transm_i = np.zeros(len(self.wvl[i]))
            dtransm_i = np.zeros(len(self.wvl[i]))
            transm_i[inds_passband] = levels[i]
            dtransm_i[inds_passband] = errs[i]
            self.transm.append(transm_i)
            self.dtransm.append(dtransm_i)

        # (Re)calculate expected signal spectra
        self.calc_spectra()

    def set_transmissivity_experimental(self, transm, wvl, dtransm=[]):
        """
        Set the transmissivity spectra of each channel as user-input arrays,
        e.g. from experimental calibration curves.

        Parameters
        ----------
            transm: list of double arrays
                Transmissivity spectra (arb. units) for each channel
            wvl: double array or list of double arrays
                Wavelengths corresponding to the transmissivity spectra. Units
                must be consistent with the specified laser wavelength.
                Must be either a single array with the same number of elements
                as each array in `transm` or a list of arrays of equal lengths 
                to the corresponding arrays in `transm`.
            dtransm: list of double arrays (optional)
                Uncertainty in the transmissivity spectra for each channel.
                Will assumed to be zero if not provided. 
        """

        # Check inputs
        n_channels = len(transm)

        transm_list = []
        for i in range(n_channels):
            transm_list.append(np.copy(transm[i]))

        wvl_list = []
        if type(wvl) == list:
            if len(wvl) != 1 and len(wvl) != n_channels:
                raise ValueError('If `wvl` is supplied as a multi-element ' \
                                 + 'list, the numer of elements must agree \n' \
                                 + 'with the length of `transm`')
            elif len(wvl) == 1:
                wvl_arr = np.array(wvl[0])
                for i in range(n_channels):
                    if len(wvl_arr) != len(transm[i]):
                        raise ValueError('Dimensions of array in `wvl` are ' \
                                  + 'not consistent with those in `transm`')
                    wvl_list.append(np.copy(wvl_arr))
            else:
                for i in range(n_channels):
                    if len(wvl[i]) != len(transm[i]):
                        raise ValueError('Dimensions of array in `wvl` are ' \
                                  + 'not consistent with those in `transm`')
                    wvl_list.append(np.array(wvl[i]))
                    
        elif type(wvl) == np.ndarray:
            for i in range(n_channels):
                if len(wvl) != len(transm[i]):
                    raise ValueError('Dimensions of array in `wvl` are ' \
                                     + 'not consistent with those in `transm`')
                wvl_list.append(np.copy(wvl))

        dtransm_list = []
        if not dtransm:
            for i in range(n_channels):
                dtransm_list.append(np.zeros(np.shape(transm[i])))
        else:
            if len(dtransm) != n_channels:
                raise ValueError('Lengths of `transm` and `dtransm` disagree')
            for i in range(n_channels):
                if len(dtransm[i]) != len(transm[i]):
                    raise ValueError('Dimensions of arrays in `dtransm` are ' \
                                     + 'not consistent with those in `transm`')
                dtransm_list.append(np.copy(dtransm[i]))
            

        # Update number of channels and arrays for wavelength and transmissivity
        self.n_channels = n_channels
        self.wvl = wvl_list
        self.transm = transm_list
        self.dtransm = dtransm_list
 
        # (Re)calculate expected signal spectra
        self.calc_spectra()

    def set_wvl_laser(self, wvl_laser):
        """
        Set the value of the wavelength of the Thomson scattering laser.
        Units must be consistent with those supplied with the transmissivity 
        spectra.
        """

        if (not np.isscalar(wvl_laser)) or (wvl_laser <= 0):
            raise ValueError('`wvl_laser` must be a positive scalar')

        self.wvl_laser = wvl_laser

        # (Re)calculate expected signal spectra if transm data are available
        self.calc_spectra()

    def set_scat_ang(self, scat_ang):
        """
        Set the angle (radians), relative to the laser beam line, of the 
        scattered light received by the polychromator.
        """

        if not np.isscalar(scat_ang):
            raise ValueError('`scat_ang` must be a scalar')

        # (Re)calculate expected signal spectra if transm data are available
        self.calc_spectra()

    def set_spec_te(self, spec_te):
        """
        Set the array of temperatures (keV) at which the spectra of expected
        signals should be evaluated.
        """

        if (not type(spec_te) == np.ndarray) or (any(spec_te <= 0)):
            raise ValueError('`spec_te` must be a Numpy array of positive ' \
                             + 'values')

        self.spec_te = spec_te

        # (Re)calculate expected signal spectra if transm data are available
        self.calc_spectra()

    def calc_spectra(self):
        """
        Calculate the expected signal spectra of each channel by integrating 
        the transmissivity spectra with the scattering spectrum (using the
        Selden formula). 

        The calculation will only be performed if all necessary data are 
        present; i.e., scattering angle, laser
        """

        # Make sure all necessary data are present
        if not (self.wvl_laser and self.scat_ang and self.wvl and self.transm):
            self.spec_calculated = False
            return   

        self.spec  = np.zeros((len(self.spec_te), self.n_channels))
        self.dspec = np.zeros((len(self.spec_te), self.n_channels))

        for i in range(self.n_channels):
            te_ev = 1.e3 * np.reshape(self.spec_te, (-1,1)) 
            wvl_mat, te_ev_mat = np.meshgrid(self.wvl[i], te_ev)
            transm_mat, _ = np.meshgrid(self.transm[i], te_ev)
            dtransm_mat, _ = np.meshgrid(self.dtransm[i], te_ev)
            x_mat = wvl_mat/self.wvl_laser - 1.0

            photons = scattered_photon_spectrum(x_mat, te_ev_mat, self.scat_ang)
            integrand = photons*transm_mat
            integrand_err = photons*dtransm_mat
            self.spec[:,i] = np.trapz(integrand, x=wvl_mat, axis=1)
            self.dspec[:,i] = integral_error(wvl_mat, integrand, integrand_err)

        self.spec_calculated = True
            

    def transmissivity(self, channel):
        """
        Return the transmissivity spectrum for a selected channel.

        Parameters
        ----------
            channel: integer
                Index number of the channel (counting from zero)

        Returns
        -------
            transm: double array
                Transmissivity spectrum 
            dtransm: double array
                Uncertainty of each value in the transmissivity spectrum
            wvl: double array
                Wavelength of each element of `transm` as specified by the user
        """

        # Verify that data are present
        if (not self.transm) or (not self.wvl):
            raise ValueError('Transmissivity spectra have not been added')

        return self.transm[channel], self.dtransm[channel], self.wvl[channel]

    def expected_relative_signal(self, channel):
        """
        Return the expected relative signal strength vs. temperature for a 
        specified channel

        Parameters
        ----------
            channel: integer
                Index number of the channel (counting from zero)

        Returns
        -------
            exp_rel_sig: double array
                Expected signal strength
            exp_rel_sig_err: double array
                Uncertainty in `exp_rel_sig`
            spec_te: double array
                Temperatures (keV) corresponding to each element of `exp_sig`
        """

        # Make sure expected signals are available
        if not self.spec_calculated:
            raise ValueError('More data are needed to calculate spectra, '\
                             'possibly wvl_laser, scat_ang, or transmissivity')

        return np.copy(self.spec[:,channel]), np.copy(self.dspec[:,channel]), \
               np.copy(self.spec_te)

    def expected_relative_signal_spectra(self):
        """
        Return the expected relative signals from all channels along an array 
        of electron temperatures

        Returns
        -------
            exp_rel_sig: 2D double array
                Array in which the element in row i, column j is the expected
                relative signal from channel j when the electron
                temperature is equal to `spec_te[i]`
            exp_rel_sig_err: 2D double array
                Uncertainty in each element of `exp_rel_sig`
            spec_te: double array
                Temperatures (keV) for each column in exp_sig
        """

        # Make sure expected signals are available
        if not self.spec_calculated:
            raise ValueError('More data are needed to calculate spectra, '\
                             'possibly wvl_laser, scat_ang, or transmissivity')

        return np.copy(self.spec), np.copy(self.dspec), np.copy(self.spec_te)

    def expected_relative_signal_interp(self, channel, te):
        """
        Estimate the expected relative signal from a given channel by 
        interpolation.

        Parameters
        ----------
            channel: integer
                Index number of the channel (counting from zero)
            te: double or double array
                Temperature (keV) at which to estimate the signal level

        Returns
        -------
            exp_sig: double or double array
                Expected signal levels for `te` if a scalar, or at each 
                temperature in `te` if an array
            exp_sig_err: double or double array
                Uncertainty in `exp_sig`
        """

        # Make sure expected signals are available
        if not self.spec_calculated:
            raise ValueError('More data are needed to calculate spectra, '\
                             'possibly wvl_laser, scat_ang, or transmissivity')

        return np.interp(te, self.spec_te, self.spec[:,channel]), \
               np.interp(te, self.spec_te, self.dspec[:,channel])

    def synthetic_signals(self, te, am, noise_std=0., signal_dt=4.e-9, \
                          signal_duration=2.e-6, pulse_width=2.e-8, \
                          pulse_center=[]):
        """
        Generates synthetic signals to simulate the output of the polychromator
        in response to a scattering signal with Gaussian form for a given
        electron temperature and amplitude multiplier.

        NOTE: the amplitude multiplier is equivalent to the scaling parameter
        for spectral fitting IF the signals are evaluated using the 'peak'
        method.

        Parameters
        ----------
            te: double
                Electron temperature (keV)
            am: double
                Amplitude multiplier (arb. units)
            noise_std: double (optional)
                Standard deviation of normally-distributed random noise 
                contribution to be added to each element of the output signal.
                Units are those of the output signal: 
                    [`s`]*[expected relative signal]
            signal_dt: double (optional)
                Time interval between data points in each signal (seconds)
            signal_duration: double (optional)
                Time duration of the simulated signal, to the nearest multiple
                of `signal_dt` (seconds)
            pulse_width: double (optional)
                Standard deviation (seconds) of thesimulated Gaussian signal
                pulses
            pulse_center: double (optional)
                Time point during the simulated signal at which the pulse 
                reaches its maximum. If no input is provided, the signal will
                be centered at 0.5*`signal_duration`.

        Returns
        -------
            signals: list of 1D Numpy arrays
                List of simulated signals from each of the polychromator
                channels
            times: 1D Numpy array
                Time associated with each element in the arrays in `signals`
        """

        # Make sure expected signals are available
        if not self.spec_calculated:
            raise ValueError('More data are needed to calculate spectra, '\
                             'possibly wvl_laser, scat_ang, or transmissivity')

        times = np.arange(0., signal_duration, signal_dt)
        n_times = len(times)
        signals = []

        if not pulse_center:
            t0 = 0.5 * signal_duration
        else:
            t0 = pulse_center

        for i in range(self.n_channels):
            relSig, _ = self.expected_relative_signal_interp(i, te)
            ampl = am * relSig
            signals.append(ampl * np.exp(-((times-t0)**2/(2*pulse_width**2))) \
                           + noise_std * np.random.randn(n_times))

        return signals, times

def selden(x, te, theta):
    """
    Calculates the Selden curve for the scattered power per unit wavelength
    due to Thomson scattering at a given angle and electron temperature for an 
    array of wavelengths.

    Source: A. C. Selden, Phys. Lett. 79A, 405 (1980)

    Parameters
    ----------
        x: scalar or Numpy array
            fractional difference between the scattered wavelength(s) at which 
            to evaluate the Selden function and the wavelength of the incident
            light: x = (lambda_scattered - lambda_incident)/lambda_incident
        te: float or double
            Electron temperature in eV
        theta: float or double
            scattering angle in radians

    Returns
    -------
        s: scalar or Numpy array
            Selden function evaluated at all wavelengths in x
    """

    alpha = 0.5 * me * c**2 / (qe * te)
    A = (1 + x)**3 * np.sqrt(2. * (1. - np.cos(theta)) * (1. + x) + x**2)
    B = np.sqrt(1. + x**2/(2. * (1. - np.cos(theta)) * (1. + x))) - 1.
    C = np.sqrt(alpha/np.pi) * (1. - 15./16./alpha + 345./512./alpha**2)
    
    return C / A * np.exp(-2. * alpha * B)

def scattered_photon_spectrum(x, te, theta):
    """
    Returns a value proportional to the number of scattered photons per unit 
    wavelength due to Thomson scattering at a given angle and electron 
    temperature for an array of wavelengths.

    This is attained by multiplying the Selden curve (see function `selden()`)
    for scattered power by the wavelength (scaled by the incident laser
    wavelength).

    Parameters
    ----------
        x: scalar or Numpy array
            fractional difference between the scattered wavelength(s) at which 
            to evaluate the Selden function and the wavelength of the incident
            light: x = (lambda_scattered - lambda_incident)/lambda_incident
        te: float or double
            Electron temperature in eV
        theta: float or double
            scattering angle in radians

    Returns
    -------
        p: scalar or Numpy array
            Selden function evaluated at all wavelengths in x
    """

    return (1. + x) * selden(x, te, theta)

def integral_error(x, y, dy):
    """
    Estimate the error of an integral (performed through trapezoidal 
    integration) in which the integrand has random uncertainties. This
    implementation is intended to correspond to integrals taken along multiple
    rows of a 2D array (i.e., each row is integrated separately)

    Parameters
    ----------
        x: 2D Numpy double array
            Values of the independent variable associated with each discrete
            value of the integrand; should increase along each row
        y: 2D Numpy double array
            Values of the integrand along each row
        dy: 2D Numpy double array
            Uncertainty in `y` (values of the integrand)

    Returns
    -------
        dInt: 1D Numpy double array
            Uncertainty in the integrals; contains one value for each row of
            the input arrays
    """

    # Note that the trapezoidal integral is a linear combination of each
    # element of the integrand; evaluate its partial derivatives with respect
    # to each component of the integrand vector
    dIdy = np.zeros(np.shape(y))
    dIdy[:,0] = 0.5 * (x[:,1] - x[:,0])
    dIdy[:,-1] = 0.5 * (x[:,-1] - x[:,-2])
    dIdy[:,1:-1] = x[:,2:] - x[:,:-2]
 
    return np.sqrt(np.sum(dIdy**2 * dy**2, axis=1))

