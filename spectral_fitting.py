"""
spectral_fitting.py

Emulates the real-time code's tool for performing a fit of polychromator
data to the values expected from Thomson scattering spectra
"""

import numpy as np


def fitSpectrum(sigOk, signals, dSignal, s0, te0, specTe, spec, dSpec, \
                saturation=np.Inf, teMin=0.1, verbose=False, \
                tol=0.001, itmax=200):
    """
    Performs a fit of polychromator data to polychromator spectral calibration 
    data with a Levenberg-Marquardt least-squares algorithm.

    Parameters
    ----------
        sigOk: logical array
            Masking array in which the n^th element is true if the n^th 
            polychromator channel signal should be included in the fit and false
            if the signal should be ignored for the fit
        signals: double array
            Signals obtained from the polychromator channels
        dSignal: double array
            Standard errors for the values in `signals`
        s0: double
            Initial guess for the scaling factor parameter for the fit
        te0: double
            Initial guess for the electron temperature (keV) parameter for the 
            fit
        specTe: double array
            Values of electron temperature (keV) associated with each row of 
            `spec` and `dSpec`
        spec: 2D double array
            Array in which the element in row i, column j is the expected
            relative signal from polychromator channel j when the electron
            temperature is equal to `specTe[i]`
        dSpec: 2D double array
            Array in which the element in row i, column j is the standard error
            in the corresponding element of `spec`
        saturation: double (optional)
            Maximum acceptable value for any polychromator channel signal;
            signals greater than this will be ignored in fits.
        teMin: double (optional)
            Minimum temperature at which fitting will be attempted
        verbose: logical (optional)
            If true, parameter values will be printed to the standard output
            at each iteration of the fit   
        tol: double (optional)
            Maximum relative change in chi^2 between successive iterations
            for a fit to be considered converged. Fitting iterations will 
            terminate once this tolerance is reached.
        itmax: integer (optional)
            Maximum allowable number of iterations

    Returns
    -------
        s: double
            Best-fit scaling factor
        te: double
            Best-fit electron temperature (keV)
        ds: double
            Standard error in `s`
        dte: double
            Standard error in `te`
        chi2: double
            Reduced chi^2 associated with the fit
        sSeries: double array
            Values of `s` from each iteration step, starting with `s0`
        teSeries: double array
            Values of `te` from each iteration step, starting with `te0`
        chi2Series: double array
            Values of `chi2` from each iteration step, starting with the initial
            guess
    """

    nterms = 2

    # Ensure that input arrays are Numpy arrays
    signals = np.array(signals)
    dSignal = np.array(dSignal)
    sigOk   = np.array(sigOk)
    spec  = np.array(spec)
    dSpec = np.array(dSpec)

    # Reject values above saturation level
    for i in range(len(sigOk)):
        if np.abs(signals[i]) > saturation or np.isnan(signals[i]) \
                                   or np.isnan(dSignal[i]):
            sigOk[i] = False
    okInds = np.where(sigOk)[0]
    nOk = len(okInds)
    nfree = nOk - nterms

    flambda = 0.001

    # Initial model, residuals, weights, and chi2
    model, jacobian = lookup(s0, te0, specTe, spec, sigOk=sigOk)  
    residuals = signals[sigOk] - model
    rel_dSpec = [np.interp(te0, specTe, dSpec[:,i]) \
                      /np.interp(te0, specTe, spec[:,i]) for i in okInds]
    weights = 1./(dSignal[sigOk]**2 + (signals[sigOk]*rel_dSpec)**2)
    chi2 = np.sum(residuals**2 * weights)/nfree

    s  = s0
    te = te0
    s_i = s0
    te_i = te0
    sSeries = [s0]
    teSeries = [te0]
    chi2Series = [chi2]

    if verbose:
        print('%5d: s = %11.3e; te = %6.3f; chi2 = %11.3e' % (0, s, te, chi2))

    for nIter in range(itmax):

        if chi2 < np.sum(signals[sigOk])/1.e7/nfree:
            break

        # grad-chi2 vector and curvature matrix
        JtWR  = np.dot(jacobian.T, weights*residuals)
        JtWJ  = np.dot(jacobian.T, np.dot(np.diag(weights), jacobian))

        # Normalize JtWJ
        normalization = np.array([[JtWJ[0,0], np.sqrt(JtWJ[0,0]*JtWJ[1,1])], \
                                  [np.sqrt(JtWJ[0,0]*JtWJ[1,1]), JtWJ[1,1]]])
        JtWJn = JtWJ/normalization
        
        # Adjust the step toward the descent direction if chi2 does not decrease
        chi2_i = chi2
        lambdaCount = 0
        done = False
        while (chi2_i >= chi2 or te_i < teMin or s_i < 0) and not done:

            lambdaCount = lambdaCount + 1

            Curv = JtWJn + flambda*np.eye(nterms)
            detCurv = Curv[0,0]*Curv[1,1] - Curv[0,1]*Curv[1,0]
            iCurv = (1./detCurv) \
                *np.array([[Curv[1,1], -Curv[0,1]], [-Curv[1,0], Curv[0,0]]])
            #iCurv = np.linalg.inv(Curv)

            params = np.dot(iCurv/normalization, JtWR)
            s_i = s + params[0]
            te_i = te + params[1]

            model_i, jacobian_i = lookup(s_i, te_i, specTe, spec, sigOk=sigOk)
            residuals_i = signals[sigOk] - model_i
            chi2_i = np.sum(residuals_i**2 * weights)/nfree

            if verbose:
                print(('%5d: s = %11.3e; te = %6.3f; chi2 = %12.6e; ' \
                          + 'chi2_0 = %12.6e') % \
                      (nIter+1, s_i, te_i, chi2_i, chi2))

            # In certain cases, reject changes made in this iteration
            if (chi2_i < 1.e-30 or chi2_i > 1.e30) or \
               (lambdaCount > 30): # and chi2_i > chi2):

                done = True

            flambda = flambda*10.

        flambda = flambda*0.01
        s  = s_i
        te = te_i
        sSeries.append(s)
        teSeries.append(te)

        # New fit metrics
        model, jacobian = lookup(s, te, specTe, spec, sigOk=sigOk)  
        residuals = signals[sigOk] - model
        rel_dSpec = [np.interp(te, specTe, dSpec[:,i]) \
                          /np.interp(te, specTe, spec[:,i]) for i in okInds]
        weights = 1./(dSignal[okInds]**2 + (signals[okInds]*rel_dSpec)**2)
        chi2_prev = chi2
        chi2 = np.sum(residuals**2 * weights)/nfree
        chi2Series.append(chi2)

        # Stop condition
        if np.abs(chi2_prev - chi2)/chi2 <= tol:
            break;

    JtWJdet = JtWJ[0,0]*JtWJ[1,1] - JtWJ[0,1]*JtWJ[1,0]
    JtWJ_inv = 1./JtWJdet \
        *np.array([[JtWJ[1,1], -JtWJ[0,1]], [-JtWJ[1,0], JtWJ[0,0]]])
    ds  = np.sqrt(JtWJ_inv[0,0])
    dte = np.sqrt(JtWJ_inv[1,1])

    return s, te, ds, dte, chi2, sSeries, teSeries, chi2Series


def lookup(s, te, specTe, spec, sigOk=np.array([True, True, True, True, True])):
    """
    Calculates expected polychromator signals and model Jacobian for a given 
    temperature and scaling factor by interpolating the calibration spectra

    Parameters
    ----------
        s: double
            Scaling factor
        te: double
            Electron temperature (keV)
        specTe: double array
            Values of electron temperature (keV) associated with each row of 
            `spec` and `dSpec`
        spec: 2D double array
            Array in which the element in row i, column j is the expected
            relative signal from polychromator channel j when the electron
            temperature is equal to `specTe[i]`
        sigOk: logical array
            Masking array in which the n^th element is true if the n^th 
            polychromator channel signal should be included in the fit and false
            if the signal should be ignored for the fit

    Returns
    -------
        model: double array
            Modeled polychromator signals, only channels that are marked as
            True in `sigOk`
        jacobian: 2D double array
            Jacobian associated with the model values: first column contains
            the derivatives of the values in `model` with respect to `s`;
            second column contains their derivatives with  respect to `te`.
    """

    s  = np.abs(s)
    te = np.abs(te)

    hte = 0.01 # interval for numerical differentiation
    if te < hte:
        hte = -hte

    okInds = np.where(sigOk)[0]
    nOk    = len(okInds)
    nTe, nPoly = np.shape(spec)
    specDist = np.zeros(nOk)
    dSpecDist = np.zeros(nOk)
    jacobian = np.zeros((nOk,2))
    for i in range(nOk):
        specDist[i] = np.interp(te, specTe, spec[:,okInds[i]])
        dSpecDist[i] = (specDist[i] \
                           - np.interp(te-hte, specTe, spec[:,okInds[i]]))/hte

    jacobian[:,0] = specDist         # dSpectrum/ds
    jacobian[:,1] = s*dSpecDist      # dSpectrum/dte

    model = s*specDist

    return model, jacobian


