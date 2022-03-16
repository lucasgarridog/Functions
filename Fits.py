import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def Linear(x, M, N):
    # Linear function
    # M: slope
    # N: intercept
    return M * x + N

def Gaussian(x, A, MU, SIGMA):
    # Gaussian function
    # A: amplitude
    # MU: mean
    # SIGMA: standard deviation
    return A * np.exp(-(x-MU)**2 / (2*SIGMA**2))

def Linear_fit(x,y):
    init = [(y[2]-y[1]) / (x[2]-x[1]), 0]                                 # initial guess of the parameters
    fit = scipy.optimize.curve_fit(Linear, x, y, init)                    # fit
    opt_param = fit[0]                                                    # optimal parameters that fit the data
    opt_param_cov = fit[1]                                                # matrix of covariance
    opt_param_error = np.sqrt(np.diag(opt_param_cov))                     # errors are in the diagonal of the matrix
    M = opt_param[0]                                                      # slope of the fit
    N = opt_param[1]                                                      # intercept of the fit
    delta_M = opt_param_error[0]                                          # error of the slope
    delta_N = opt_param_error[1]                                          # error of the intercept
    values = {'slope': M, '\Delta(slope)': delta_M,                       # the fit is returned in a dictionary
              'intercept': N, '\Delta(intercept)': delta_N
              }
    return values

def Gaussian_fit(x,y):
    init = [max(y), x[0], (x[1] - x[0]) * 5]                                # initial guess of the parameters
    fit = scipy.optimize.curve_fit(Gaussian, x, y, init)                    # fit
    opt_param = fit[0]                                                      # optimal parameters that fit the data
    opt_param_cov = fit[1]                                                  # matrix of covariance
    opt_param_error = np.sqrt(np.diag(opt_param_cov))                       # errors are in the diagonal of the matrix
    A = opt_param[0]                                                        # amplitude of the fit
    MU = opt_param[1]                                                       # mean of the fit
    SIGMA = opt_param[2]                                                    # std dev of the fit
    delta_A = opt_param_error[0]                                            # error of the amplitude
    delta_MU = opt_param_error[1]                                           # error of the mean
    delta_SIGMA = opt_param_error[2]                                        # error of the std dev
    FWHM = 2 * SIGMA * np.sqrt(2 * np.log(2))                               # Full Width at Half Maximum
    delta_FWHM = 2 * np.sqrt(2 * np.log(2)) * delta_SIGMA                   # error of FWHM
    R = 100 * FWHM / MU                                                     # Resolution
    delta_R = R * np.sqrt((delta_FWHM / FWHM) ** 2 + (delta_MU / MU) ** 2)  # error of R
    values = {'amplitude': A, '\Delta(amplitude)': delta_A,                 # the fit is returned in a dictionary
              'mean': MU, '\Delta(mean)': delta_MU,
              'sigma': SIGMA, '\Delta(sigma)': delta_SIGMA,
              'FWHM': FWHM, '\Delta(FWHM)': delta_SIGMA,
              'R[%]': R, '\Delta(R[%])': delta_R
              }
    return values

def Gaussian_fit_spectra(x,y,plot=True):
    peaks, _ = find_peaks(y, prominence=300)                 # find the gaussian peaks
    amplitudes = np.zeros((len(peaks),2))                    # amplitudes with their errors will be saved here
    means = np.zeros((len(peaks),2))                         # means with their errors will be saved here
    sigmas = np.zeros((len(peaks),2))                        # std devs with their errors will be saved here
    fwhms = np.zeros((len(peaks),2))                         # fwhms with their errors will be saved here
    resolutions = np.zeros((len(peaks),2))                   # resolutions with their errors will be saved here
    for i in range(len(peaks)):                                          # fit each peak
        init = [y[peaks[i]], x[peaks[i]], (x[1]-x[0])*5]                 # initial guess of the parameters
        k = peaks[i]
        while y[k] > int(y[peaks[i]]/ 2):                                # iterate until it reaches the half height
            k = k+1

        dist = k - peaks[i]                                              # distance between the peak and its half height
        fit_x = x[peaks[i] - dist:peaks[i] + dist]                       # x fit range
        fit_y = y[peaks[i] - dist:peaks[i] + dist]                       # y fit range
        fit = scipy.optimize.curve_fit(Gaussian, fit_x, fit_y, init)     # fit
        opt_param = fit[0]                                               # optimal parameters that fit the data
        opt_param_cov = fit[1]                                           # matrix of covariance
        opt_param_error = np.sqrt(np.diag(opt_param_cov))                # errors are in the diagonal of the matrix
        amplitudes[i,0] = opt_param[0]                                   # amplitude of the fit
        means[i,0] = opt_param[1]                                        # mean of the fit
        sigmas[i,0] = opt_param[2]                                       # std dev of the fit
        amplitudes[i,1] = opt_param_error[0]                             # error of the amplitude
        means[i,1] = opt_param_error[1]                                  # error of the mean
        sigmas[i,1] = opt_param_error[2]                                 # error of the std dev
        fwhms[i,0] = 2*sigmas[i,0]*np.sqrt(2*np.log(2))                  # Full Width at Half Maximum
        fwhms[i,1] = 2*np.sqrt(2*np.log(2)) * sigmas[i,1]                # error of FWHM
        resolutions[i,0] = 100*fwhms[i,0] / means[i,0]                   # Resolution
        resolutions[i,1] = resolutions[i,0] *\
                           np.sqrt( (fwhms[i,1] / fwhms[i,0])**2 + (means[i,1] / means[i,0])**2 )      # error of R

        if plot == True:
            plt.figure(1)                                         # Plotting
            plt.step(x, y, color="tab:blue", zorder=0)
            plt.fill_between(x, y, step="pre", color="tab:blue", zorder=0)
            new_x = x[peaks[i] - 4*dist:peaks[i] + 4*dist]
            plt.plot(new_x, Gaussian(new_x, amplitudes[i,0], means[i,0], sigmas[i,0]), color="tab:red", zorder=1)
            plt.ylim(0)

    plt.show()
    values = {'amplitude': amplitudes,              # the fit is returned in a dictionary
              'mean': means,
              'sigma':sigmas,
              'FWHM': fwhms,
              'R[%]': resolutions
              }
    return values