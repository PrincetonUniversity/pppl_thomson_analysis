"""
pulse_eval.py

Procedures for calculating an intensity metric for a time-resolved laser
pulse signal output by a polychromator. The pulses are assumed to be
negative, such that the "peak" value is the minimum value.
"""

import numpy as np

# Default values for signal fitting parameters
WINDOW_HALF_WIDTH = 20  
APPROX_2ND_PEAK_OFFS = -7
SNR_MIN = 2.0
N_BACKGND_AVG = 20
N_CLIP = 6
DELTA_T = 0.25
POLARITY = 1
FIT_TOL = 0.001

NP_SINGLE_GAUSS = 4
NP_DOUBLE_GAUSS = 7

def processPulse(data, method='fit', isDouble=False, **kwargs):
    '''
    Calculates a scalar strength parameter from a raw, time-resolved 
    polychromator signal according to a user-specified method. Currently 
    supports peak detection, direct (trapezoidal) integration, and Gaussian 
    curve-fitting.

    Parameters
    ----------
        data: int or double array
            Raw signal from one polychromator channel from one laser pulse
        method: string
            'peak': determine the peak value relative to the background
            'integrate': determine the integral of the signal through numerical
                integration
            'fit': determine the integral of the signal through Gaussian 
                curve-fitting
        isDouble: Boolean
            If true, will be fit as a double Gaussian; otherwise, will be fit
            as a single gaussian. Only relevant if `method=='fit'`.

        kwargs:
            window_half_width: integer
                Half the number of elements in the data array that must fully
                contain the laser pulse signal. 
            approx_2nd_peak_offs: integer
                Expected/approximate number of elements by which a second peak
                (if it exists) would follow the first peak
            snr_min: double
                Minimum signal-to-noise ratio for an intensity metric to be
                determined. Otherwise, `signal` will be returned as nan.
            n_backgnd_avg: integer
                Number of elements of `data` to be used to compute the mean
                background level (evaluated outside of the pulse window)
            n_clip: integer
                Number of elemenents at the end of `data` to be ignored when
                searching for the peak value
            delta_t: double
                Time interval between subsequent points in `data` to use
                for integrals
            fit_tol: double
                Relative tolerance to determine convergence of fits to the
                curve (default is 0.001). Only relevant if `method=='fit'`.

    Returns
    -------
        signal: double
            Intensity metric for the time-resolved input signal
        err: double
            Standard error in the value of `signal`
        info: dictionary with the following key/value pairs:
            ind_peak: integer
                Index of the location of the (negative) peak value of `data`
            params: double array (`'fit'` mode only)
                Optimized fitting parameters
            dparams: double array (`'fit'` mode only)
                Standard error in optimized fitting parameters 
    '''

    # Name default signal fitting parameters (sfp) and update with kwargs
    sfp = dict(window_half_width = WINDOW_HALF_WIDTH, \
               approx_2nd_peak_offs = APPROX_2ND_PEAK_OFFS, \
               snr_min = SNR_MIN, \
               n_backgnd_avg = N_BACKGND_AVG, \
               n_clip = N_CLIP, \
               delta_t = DELTA_T, \
               polarity = POLARITY, \
               fit_tol = FIT_TOL)
    sfp.update(kwargs)

    peak, x, pulseSig, ampBkg, meanBkg, peakIndex = \
        preProcessPeak(data, **kwargs)

    if np.abs(peak)/ampBkg < sfp['snr_min']:
        pulse_ok = False
    else:
        pulse_ok = True

    if method == 'peak':

        if not pulse_ok:
            return np.nan, ampBkg, dict(indPeak=peakIndex)

        return peak, ampBkg, dict(indPeak=peakIndex)

    elif method == 'integrate':

        if not pulse_ok:
            return np.nan, ampBkg, dict(indPeak=peakIndex)

        # Integrated value
        dx = x[1:]-x[:-1]
        integrand = sfp['polarity']*(pulseSig - meanBkg)
        integral = 0.5*np.sum(dx*(integrand[:-1]+integrand[1:]))
        integr_err = np.sqrt(len(integrand))*sfp['delta_t']*ampBkg

        return integral, integr_err, dict(indPeak=peakIndex)

    elif method == 'fit':

        peak_std_guess = 0.1 * sfp['delta_t']*sfp['window_half_width']

        if not isDouble:

            if not pulse_ok:
                return np.nan, np.nan, dict(indPeak=peakIndex, \
                       params=np.ones((4,1))*np.nan, \
                       dparams=np.ones(4)*np.nan)

            params0 = [sfp['polarity']*float(peak), peak_std_guess, \
                       x[sfp['window_half_width']], meanBkg]

            params, errs = fitPulse(params0, x, pulseSig, singleGaussian, \
                               jacobian_singleGaussian, fit_tol=sfp['fit_tol'])

            integral = np.sqrt(2.*np.pi) * np.abs(params[0]*params[1])
            err = np.sqrt(2.*np.pi) * np.sqrt((params[0]*errs[1])**2 \
                                              + (params[1]*errs[0])**2)

            return integral[0], err[0], \
                   dict(indPeak=peakIndex, params=params, dparams=errs)

        else:

            if not pulse_ok:
                return np.nan, np.nan, dict(indPeak=peakIndex, \
                       params=np.ones((7,1))*np.nan, \
                       dparams=np.ones(7)*np.nan)

            params0 = [sfp['polarity']*float(peak), peak_std_guess, \
                       x[sfp['window_half_width']], \
                       0.5*sfp['polarity']*float(peak), peak_std_guess, \
                       x[sfp['window_half_width'] \
                           + sfp['approx_2nd_peak_offs']], \
                       meanBkg]

            params, errs = fitPulse(params0, x, pulseSig, doubleGaussian, \
                               jacobian_doubleGaussian, fit_tol=sfp['fit_tol'])

            integral = np.sqrt(2.*np.pi) * (abs(params[0]*params[1]) + \
                                            abs(params[3]*params[4]))
            err = np.sqrt(2.*np.pi) * np.sqrt((params[0]*errs[1])**2 + \
                                              (params[1]*errs[0])**2 + \
                                              (params[3]*errs[4])**2 + \
                                              (params[4]*errs[3])**2)

            return integral[0], err[0], \
                   dict(peakInd=peakIndex, params=params, dparams=errs)

def preProcessPeak(data, **kwargs):
    """
    Internal subroutine for assessing key signal characteristics
    """

    # Name default signal fitting parameters (sfp) and update with kwargs
    sfp = dict(window_half_width = WINDOW_HALF_WIDTH, \
               approx_2nd_peak_offs = APPROX_2ND_PEAK_OFFS, \
               snr_min = SNR_MIN, \
               n_backgnd_avg = N_BACKGND_AVG, \
               n_clip = N_CLIP, \
               delta_t = DELTA_T, \
               polarity = POLARITY)
    sfp.update(kwargs)

    # Find the minimum value and its location (index)
    if sfp['polarity'] == -1:
        peakValue = np.min(data[:-sfp['n_clip']])
        peakIndex = np.argmin(data[:-sfp['n_clip']])
    elif sfp['polarity'] == 1:
        peakValue = np.max(data[:-sfp['n_clip']])
        peakIndex = np.argmax(data[:-sfp['n_clip']])
    else:
        raise ValueError('polarity parameter must be 1 or -1')

    # Starting index for calculating the average background
    if peakIndex - sfp['window_half_width'] - sfp['n_backgnd_avg'] > 0:
        newIndex = peakIndex - sfp['window_half_width'] - sfp['n_backgnd_avg']
    elif peakIndex + sfp['window_half_width'] + sfp['n_backgnd_avg'] \
             < len(data) - sfp['n_clip']:
        newIndex = peakIndex + sfp['window_half_width']
    else:
        raise RuntimeError('Sampling window too short for desired ' + \
                           'background averaging window')

    # Determine the index to integrate or fit around (if close to buffer edge)
    if peakIndex < sfp['window_half_width']:
        ctrIndex = sfp['window_half_width']
    elif len(data) - 1 - peakIndex < sfp['window_half_width']:
        ctrIndex = len(data) - 1 - sfp['window_half_width']
    else:
        ctrIndex = peakIndex

    # Calculate average and lowest noise level
    minBkg  = np.min( data[newIndex:(newIndex+sfp['n_backgnd_avg'])])
    maxBkg  = np.max( data[newIndex:(newIndex+sfp['n_backgnd_avg'])])
    ampBkg  = 0.5*(maxBkg - minBkg)
    meanBkg = np.mean(data[newIndex:(newIndex+sfp['n_backgnd_avg'])])

    # Absolute peak relative to mean background
    if sfp['polarity'] == -1:
        peak = meanBkg - peakValue
    elif sfp['polarity'] == 1:
        peak = peakValue - meanBkg

    # Subset of the signal around the peak
    peakStart = ctrIndex - sfp['window_half_width']
    peakStop  = ctrIndex + sfp['window_half_width']
    signalTime = sfp['delta_t']*(np.arange(peakStart,peakStop))
    tPeak = sfp['delta_t']*peakIndex
    pulseSignal = data[peakStart:peakStop]

    return peak, signalTime, pulseSignal, ampBkg, meanBkg, peakIndex

def fitPulse(params0, x, data, model, jacobian, verbose=0, \
             itmax=20, lm_init_deg=-2, lm_max_deg=6, fit_tol=FIT_TOL):
    """
    Uses the Levenberg-Marquardt least-squares fitting algorithm to find
    best-fit parameters for a user-specified model to a dataset.

    Parameters
    ----------
        params0: double array
            Initial guesses for each of the fitting parameters
        x: 1D double array
            Independent variable associated with the dataset
        data: 1D double array
            Dependent variable associated with the dataset
        model: function with arguments (`params`, `x`)
            Must return an array of output values corresponding to each
            value of the independent variable `x`, according to the parameters
            in `params`. For the fit to succeed, output values should vary
            smoothly with each parameter in `params`.
        jacobian: function with arguments (`params`, `x`)
            Must return a 2D array in which the element in (row i, column j) is
            the partial derivative of the `model(params,x[i])` with respect
            to `params[j]`.
        verbose: logical (optional)
            If true, parameter values will be printed to the standard output
            at each iteration of the fit   
        itmax: integer (optional)
            Maximum allowable number of iterations
        lm_init_deg: integer (optional)
            Initial base-10 logarithm of the Levenberg-Marquardt damping factor
            for the inner iteration loop.
        lm_max_deg: integer (optional)
            Maximum allowable base-10 logarithm of the Levenberg-Marquardt
            damping factor for the inner iteration loop.
        fit_tol: double (optional)
            Maximum relative change in chi^2 between successive iterations
            for a fit to be considered converged. Fitting iterations will 
            terminate once this tolerance is reached.

    Returns
    -------
        params: double array
            Optimized values for the fitting parameters
        errs: double array
            Standard errors for each element of `params`
    """

    data_arr = np.reshape(np.array(data), (-1,1))
    x_arr = np.reshape(np.array(x), (-1,1))
    params0_vec = np.reshape(np.array(params0),(-1,1))

    resid = data_arr - model(params0_vec, x_arr)
    params = params0_vec
    chi2 = np.sum(resid**2)
    chi2_prev = chi2

    for i in range(itmax):

        if verbose > 0:
            print('Iter %d: p = [%.3f, %.3f, %.3f, %.3f], chi2 = %.3e' \
                  % (i, params[0], params[1], params[2], params[3], chi2))

        J = jacobian(params, x)
        JtJ = np.dot(J.T, J)
        gradChi2 = np.dot(J.T, resid)

        # Adjust the step toward the descent direction if chi2 does not decrease
        lm_deg = lm_init_deg
        while chi2 >= chi2_prev and lm_deg <= lm_max_deg:

            diagInds = np.diag_indices_from(JtJ)
            if lm_deg > lm_init_deg:
                JtJ[diagInds] = JtJ[diagInds] * (1.0 + 10.**lm_deg)

            JtJinv = np.linalg.inv(JtJ)

            deltap = np.dot(JtJinv, gradChi2)

            test_params = params + deltap
            resid = data_arr - model(test_params, x_arr)
            chi2 = np.sum(resid**2)

            lm_deg = lm_deg + 1

        params = test_params

        if (np.abs(chi2 - chi2_prev)/chi2 < fit_tol) or chi2 < 1.e-23:
            break

        chi2_prev = chi2

    covar = chi2/len(x_arr)
    errs = np.reshape(np.sqrt(covar * np.diag(JtJinv)), (-1,1))

    return params, errs

def singleGaussian(params, x):
    """
    Computes the Gaussian function.

    Parameters
    ----------
        params: double column vector (4x1 array)
            Model parameters defined as follows:
                params[0,0]: amplitude
                params[1,0]: standard deviation
                params[2,0]: center location
                params[3,0]: constant offset
        x: 1D double array
            Values of the independent variable at which to compute the function 

    Returns
    -------
        f: 1D double array
            Values of the singleGaussian computed at each element of `x`
    """


    return params[0,0]*np.exp(-(x-params[2,0])**2/(2.0*params[1,0]**2)) \
           + params[3,0]

def doubleGaussian(params, x):
    """
    Computes a double Gaussian function (the sum of two Gaussian functions).

    Parameters
    ----------
        params: double column vector (7x1 array)
            Model parameters defined as follows:
                params[0,0]: amplitude of the first Gaussian
                params[1,0]: standard deviation of the first Gaussian
                params[2,0]: center location of the first Gaussian   
                params[3,0]: amplitude of the second Gaussian
                params[4,0]: standard deviation of the second Gaussian
                params[5,0]: center location of the second Gaussian   
                params[6,0]: constant offset
        x: 1D double array
            Values of the independent variable at which to compute the function 

    Returns
    -------
        f: 1D double array
            Values of the doubleGaussian computed at each element of `x`
    """

    return params[0,0]*np.exp(-(x-params[2,0])**2/(2.0*params[1,0]**2)) \
           + params[3,0]*np.exp(-(x-params[5,0])**2/(2.0*params[4,0]**2)) \
           + params[6,0]

def jacobian_singleGaussian(params, x):
    """
    Computes the Jacobian matrix for the function singleGaussian(). 

    Parameters
    ----------
        params: double array
            Model parameters as defined in the documentation for 
            singleGaussian()
        x: 1D double array
            Values of the independent variable at which to compute the Jacobian

    Returns
    -------
        jacobian: 2D double array
            Matrix with the form df(x[i])/dp_j; i.e. the element in (row i,
            column j) is the partial derivative of singleGaussian with respect 
            to parameter j, evaluated at the i^th component of x
    """


    jacobian = np.zeros((len(x), NP_SINGLE_GAUSS))

    expTerm = np.exp(-(x-params[2,0])**2/(2.*params[1,0]**2))

    jacobian[:,0] = expTerm
    jacobian[:,1] = params[0,0] * (x-params[2,0])**2 / params[1,0]**3 * expTerm
    jacobian[:,2] = params[0,0] * (x-params[2,0]) / params[1,0]**2 * expTerm
    jacobian[:,3] = 1.0

    return jacobian

def jacobian_doubleGaussian(params, x):
    """
    Computes the Jacobian matrix for the function doubleGaussian(). 

    Parameters
    ----------
        params: double array
            Model parameters as defined in the documentation for 
            singleGaussian()
        x: 1D double array
            Values of the independent variable at which to compute the Jacobian

    Returns
    -------
        jacobian: 2D double array
            Matrix with the form df(x[i])/dp_j; i.e. the element in (row i,
            column j) is the partial derivative of doubleGaussian with respect 
            to parameter j, evaluated at the i^th component of x
    """

    jacobian = np.zeros((len(x), NP_DOUBLE_GAUSS))

    expTerm1 = np.exp(-(x-params[2,0])**2/(2.*params[1,0]**2))
    expTerm2 = np.exp(-(x-params[5,0])**2/(2.*params[4,0]**2))

    jacobian[:,0] = expTerm1
    jacobian[:,1] = params[0,0] * (x-params[2,0])**2 / params[1,0]**3 * expTerm1
    jacobian[:,2] = params[0,0] * (x-params[2,0]) / params[1,0]**2 * expTerm1
    jacobian[:,3] = expTerm2
    jacobian[:,4] = params[0,0] * (x-params[5,0])**2 / params[4,0]**3 * expTerm2
    jacobian[:,5] = params[0,0] * (x-params[5,0]) / params[4,0]**2 * expTerm2
    jacobian[:,6] = 1.0

    return jacobian
