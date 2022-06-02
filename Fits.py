import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths

# This file contains several fit routines and functions

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

def Linear_fit(x,y, plot=False):
    x = np.array(x)
    y = np.array(y)
    init = [(y[2]-y[1]) / (x[2]-x[1]), 0]                                 # initial guess of the parameters
    fit = scipy.optimize.curve_fit(Linear, x, y, init)                    # fit
    opt_param = fit[0]                                                    # optimal parameters that fit the data
    opt_param_cov = fit[1]                                                # matrix of covariance
    opt_param_error = np.sqrt(np.diag(opt_param_cov))                     # errors are in the diagonal of the matrix
    M = opt_param[0]                                                      # slope of the fit
    N = opt_param[1]                                                      # intercept of the fit
    delta_M = opt_param_error[0]                                          # error of the slope
    delta_N = opt_param_error[1]                                          # error of the intercept
    residuals = y - Linear(x, M, N)                                       # calculate residuals
    ss_res = np.sum(residuals ** 2)                                       # sum of the squared residuals
    ss_tot = np.sum((y - np.mean(y)) ** 2)                                # total sum of squares
    r_squared = 1 - (ss_res / ss_tot)                                     # r^2 (goodness of fit)
    values = {'slope': M, 'delta_slope': delta_M,                         # the fit is returned in a dictionary
              'intercept': N, 'delta_intercept': delta_N,
              'r_squared': r_squared
              }
    if plot:                                                              # plots the regression line and data
        new_x = np.linspace(x[0], x[-1])                                  # x for the regression line
        fig = plt.figure(1)
        plt.plot(x, y, "kx", zorder=1)
        plt.plot(new_x, Linear(new_x, values.get("slope"), values.get("intercept")), "r--", zorder=0)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Calibration")
        M = values.get("slope")
        M_err = values.get("delta_slope")
        N = values.get("intercept")
        N_err = values.get("delta_intercept")
        r_2 = values.get("r_squared")
        text = "m = " + "%.3f" % M + " $\pm$ " + "%.3f" % M_err + "\n" + "n = " + "%.3f" % N + " $\pm$ " + "%.3f" % N_err + "\n" + "$r^2$ = " + "%.4f" % r_2
        plt.text(0.15, 0.74, text, transform=fig.transFigure, bbox=dict(facecolor="white"))
        plt.show()
    return values

def Gaussian_fit(x,y):
    init = [max(y), x[0], (x[1] - x[0]) * 5]                                # initial guess of the parameters
    fit = scipy.optimize.curve_fit(Gaussian, x, y, init)                    # fit
    opt_param = fit[0]                                                      # optimal parameters that fit the data
    opt_param_cov = fit[1]                                                  # matrix of covariance
    opt_param_error = np.sqrt(np.diag(opt_param_cov))                       # errors are in the diagonal of the matrix
    A = opt_param[0]                                                        # amplitude of the fit
    MU = opt_param[1]                                                       # mean of the fit
    SIGMA = abs(opt_param[2])                                               # std dev of the fit (abs value)
    delta_A = opt_param_error[0]                                            # error of the amplitude
    delta_MU = opt_param_error[1]                                           # error of the mean
    delta_SIGMA = opt_param_error[2]                                        # error of the std dev
    FWHM = 2 * SIGMA * np.sqrt(2 * np.log(2))                               # Full Width at Half Maximum
    delta_FWHM = 2 * np.sqrt(2 * np.log(2)) * delta_SIGMA                   # error of FWHM
    R = 100 * FWHM / MU                                                     # Resolution
    delta_R = R * np.sqrt((delta_FWHM / FWHM) ** 2 + (delta_MU / MU) ** 2)  # error of R
    values = {'amplitude': A, 'delta_amplitude': delta_A,                 # the fit is returned in a dictionary
              'mean': MU, 'delta_mean': delta_MU,
              'sigma': SIGMA, 'delta_sigma': delta_SIGMA,
              'FWHM': FWHM, 'delta_FWHM': delta_SIGMA,
              'R[%]': R, 'delta_R[%]': delta_R
              }
    return values