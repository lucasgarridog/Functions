from Fits import *
from matplotlib import rc
rc("text", usetex=True)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['figure.figsize'] = (9, 6)

class Spectrum:
    def __init__(self, xvals, yvals, is_calibrated=False):
        self.is_calibrated = is_calibrated
        if is_calibrated:
            self.Evals = xvals                    # x values of the spectrum (E)
        self.xvals = xvals                        # x values of the spectrum (ch)
        self.yvals = yvals                        # y values of the spectrum (height of each bin)

    def peaks(self, prominence=200, distance=10, width=1):
        """ Finds peaks in the spectrum \n
            returns: peaks indexes in x """
        peaks, _ = find_peaks(self.yvals, prominence=prominence, distance=distance, width=width) # change prominence if the fit fails
        return peaks

    def calibrate(self, energies, plot=False, unit=0, prominence=300):
        """ Calibrates the spectrum using the given energies and the peaks found \n
            unit = 0 -> keV (default) \n
            unit = 1 -> MeV \n
            plot = True -> plots linear regression \n
            returns: calibrated x, dictionary with the fit parameters & information """
        if self.is_calibrated:
            print("Warning: This spectrum is already calibrated")
        energy_unit = ["keV", "MeV"]                                # available units for the x-axis
        peaks = self.peaks(prominence=prominence)                   # finds peaks in the spectrum
        peaks_channel = [self.xvals[i] for i in peaks]              # channel where the peaks are
        fit = Linear_fit(peaks_channel, energies)                   # performs linear calibration
        calibrated_x = [fit.get("slope") * i + fit.get("intercept") for i in self.xvals]   # calibrate x axis
        self.Evals = calibrated_x                                   # saves calibrated values
        self.is_calibrated = True                                   # changes the calibration attribute
        if plot:                                                    # plots the regression line and data
            x = np.linspace(peaks_channel[0], peaks_channel[-1])    # x for the regression line
            fig = plt.figure(1)
            plt.plot(peaks_channel, energies, "kx", zorder=1)
            plt.plot(x, Linear(x, fit.get("slope"), fit.get("intercept")), "r--", zorder=0)
            plt.xlabel("Channel")
            plt.ylabel("Energy ("+energy_unit[unit]+")")
            plt.title("Calibration")
            M = fit.get("slope")
            M_err = fit.get("delta_slope")
            N = fit.get("intercept")
            N_err = fit.get("delta_intercept")
            r_2 = fit.get("r_squared")
            text = "m = " + "%.2f" % M + " $\pm$ " + "%.2f" % M_err + "\n" + "n = " + "%.0f" % N + " $\pm$ " + "%.0f" % N_err + "\n" + "$r^2$ = " + "%.6f" % r_2
            plt.text(0.15, 0.74, text, transform=fig.transFigure, bbox=dict(facecolor="white"))
            plt.show()

        return calibrated_x, fit

    def fit_peaks(self,prominence=200, distance=10, width=1):
        """ Finds peaks in the spectrum and fits them \n
            returns: dictionary with fit parameters & information """
        if self.is_calibrated:
            x = self.Evals
        else:
            x = self.xvals
        peaks = self.peaks(prominence=prominence,distance=distance,width=width)   # finds peaks in the spectrum
        widths = peak_widths(self.yvals, peaks, rel_height=0.5)        # calculates the width of each peak (at 50% height)
        amplitudes = np.zeros((len(peaks), 2))                         # amplitudes with their errors will be saved here
        means = np.zeros((len(peaks), 2))                              # means with their errors will be saved here
        sigmas = np.zeros((len(peaks), 2))                             # std devs with their errors will be saved here
        fwhms = np.zeros((len(peaks), 2))                              # fwhms with their errors will be saved here
        resolutions = np.zeros((len(peaks), 2))                        # resolutions with their errors will be saved here
        for i in range(len(peaks)):                                    # fit each peak
            fit_x = x[peaks[i] - int(2*widths[0][i]):peaks[i] + int(2*widths[0][i])]              # x fit range
            fit_y = self.yvals[peaks[i] - int(2*widths[0][i]):peaks[i] + int(2*widths[0][i])]     # y fit range
            init = [self.yvals[peaks[i]], x[peaks[i]], fit_x[-1]- fit_x[0]]                       # initial guess of the parameters
            fit = scipy.optimize.curve_fit(Gaussian, fit_x, fit_y, init)                      # fit
            opt_param = fit[0]                                       # optimal parameters that fit the data
            opt_param_cov = fit[1]                                   # matrix of covariance
            opt_param_error = np.sqrt(np.diag(opt_param_cov))        # errors are in the diagonal of the matrix
            amplitudes[i, 0] = opt_param[0]                          # amplitude of the fit
            means[i, 0] = opt_param[1]                               # mean of the fit
            sigmas[i, 0] = abs(opt_param[2])                         # std dev of the fit
            amplitudes[i, 1] = opt_param_error[0]                    # error of the amplitude
            means[i, 1] = opt_param_error[1]                         # error of the mean
            sigmas[i, 1] = opt_param_error[2]                        # error of the std dev
            fwhms[i, 0] = 2 * sigmas[i, 0] * np.sqrt(2 * np.log(2))  # Full Width at Half Maximum
            fwhms[i, 1] = 2 * np.sqrt(2 * np.log(2)) * sigmas[i, 1]  # error of FWHM
            resolutions[i, 0] = 100 * fwhms[i, 0] / means[i, 0]      # Resolution
            resolutions[i, 1] = resolutions[i, 0] * \
                                np.sqrt((fwhms[i, 1] / fwhms[i, 0]) ** 2 + (means[i, 1] / means[i, 0]) ** 2)  # error of R

        values = {'amplitude': amplitudes,                           # the fit is returned in a dictionary
                      'mean': means,
                      'sigma': sigmas,
                      'FWHM': fwhms,
                      'R[%]': resolutions
                      }
        return values

    def plot(self, xlim=None, unit=0):
        """ Plots the spectrum \n
            unit = 0 -> keV (default) \n
            unit = 1 -> MeV \n """
        if self.is_calibrated:
            x = self.Evals
        else:
            x = self.xvals
        energy_unit = ["keV", "MeV"]
        # plt.figure()
        plt.step(x, self.yvals, color="tab:blue")
        plt.fill_between(x, self.yvals, color="tab:blue", step="pre")
        plt.ylim(0)
        plt.xlim(x[0], x[-1])
        plt.ylabel("Counts")
        plt.title("Spectrum")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        if self.is_calibrated:
            plt.xlabel("Energy ("+energy_unit[unit]+")")
            plt.xlim(xlim)
        else:
            plt.xlabel("ADC Channel")
            plt.xlim(xlim)
        # plt.show()
        return

    def plot_fit(self, prominence=200, distance=10, width=1, xlim=None, unit=0):
        """ Plots the spectrum and the gaussian fits to its peaks \n
            unit = 0 -> keV (default) \n
            unit = 1 -> MeV \n """
        if self.is_calibrated:
            x = self.Evals
        else:
            x = self.xvals
        peaks = self.peaks(prominence=prominence, distance=distance, width=width)     # finds peaks
        widths = peak_widths(self.yvals, peaks, rel_height=0.5)  # calculates the width of each peak (at 50% height)
        energy_unit = ["keV", "MeV"]
        plt.figure(1)
        plt.step(x, self.yvals, color="tab:blue", where="mid")
        plt.fill_between(x, self.yvals, color="tab:blue", step="mid")
        plt.ylim(0)
        plt.ylabel("Counts")
        plt.title("Spectrum + Fits")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        if self.is_calibrated:
            plt.xlabel("Energy (" + energy_unit[unit] + ")")
            plt.xlim(xlim)
        else:
            plt.xlabel("ADC Channel")
            plt.xlim(xlim)
        fit = self.fit_peaks(prominence=prominence, distance=distance, width=width)      # fit
        for i in range(len(peaks)):
            fit_x = x[peaks[i] - int(2*widths[0][i]):peaks[i] + int(2*widths[0][i])]     # x fit range
            plt.plot(x[peaks[i]], self.yvals[peaks[i]], "kx")
            plt.plot(np.linspace(fit_x[0],fit_x[-1], num=200), Gaussian(np.linspace(fit_x[0],fit_x[-1], num=200), fit.get("amplitude")[i, 0], fit.get("mean")[i, 0], fit.get("sigma")[i, 0]),color="tab:red")

        plt.ylim(0)
        # plt.show()
        return fit