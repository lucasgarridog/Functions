from Fits import *
from matplotlib import rc
rc("text", usetex=True)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['figure.figsize'] = (8, 5)

class Spectrum:
    def __init__(self, xvals, yvals):
        self.xvals = xvals
        self.yvals = yvals
        self.channelvals = xvals
        self.iscalibrated = False

    def peaks(self):
        peaks_index, _ = find_peaks(self.yvals, prominence=300)  # find the gaussian peaks
        return peaks_index

    def calibrate(self, energies, plot=False, unit=0):
        """ Unit = 0 -> keV (default) \n
            Unit = 1 -> MeV \n
            plot = True -> plots linear regression """
        energy_unit = ["keV", "MeV"]
        peaks = self.peaks()
        peaks_channel = [self.xvals[i] for i in peaks]  # channel where the peaks are
        fit = Linear_fit(peaks_channel, energies)  # performs linear calibration
        calibrated_x = [fit.get("slope") * i + fit.get("intercept") for i in self.xvals]  # calibrate x axis
        self.xvals = calibrated_x
        self.iscalibrated = True
        if plot:
            plt.figure(1)
            plt.plot(peaks_channel, energies, "kx")
            plt.xlabel("Channel")
            plt.ylabel("Energy ("+energy_unit[unit]+")")
            plt.show()

        return calibrated_x, fit

    def fit_peaks(self):
        peaks = self.peaks()
        widths = peak_widths(self.yvals, peaks, rel_height=0.8)  # calculates the width of each peak (at 20% height)
        amplitudes = np.zeros((len(peaks), 2))  # amplitudes with their errors will be saved here
        means = np.zeros((len(peaks), 2))  # means with their errors will be saved here
        sigmas = np.zeros((len(peaks), 2))  # std devs with their errors will be saved here
        fwhms = np.zeros((len(peaks), 2))  # fwhms with their errors will be saved here
        resolutions = np.zeros((len(peaks), 2))  # resolutions with their errors will be saved here
        for i in range(len(peaks)):  # fit each peak
            init = [self.yvals[peaks[i]], self.xvals[peaks[i]], (self.xvals[1] - self.xvals[0]) * 5]  # initial guess of the parameters
            fit_x = self.xvals[peaks[i] - int(widths[0][i]):peaks[i] + int(widths[0][i])]  # x fit range
            fit_y = self.yvals[peaks[i] - int(widths[0][i]):peaks[i] + int(widths[0][i])]  # y fit range
            fit = scipy.optimize.curve_fit(Gaussian, fit_x, fit_y, init)  # fit
            opt_param = fit[0]  # optimal parameters that fit the data
            opt_param_cov = fit[1]  # matrix of covariance
            opt_param_error = np.sqrt(np.diag(opt_param_cov))  # errors are in the diagonal of the matrix
            amplitudes[i, 0] = opt_param[0]  # amplitude of the fit
            means[i, 0] = opt_param[1]  # mean of the fit
            sigmas[i, 0] = abs(opt_param[2])  # std dev of the fit
            amplitudes[i, 1] = opt_param_error[0]  # error of the amplitude
            means[i, 1] = opt_param_error[1]  # error of the mean
            sigmas[i, 1] = opt_param_error[2]  # error of the std dev
            fwhms[i, 0] = 2 * sigmas[i, 0] * np.sqrt(2 * np.log(2))  # Full Width at Half Maximum
            fwhms[i, 1] = 2 * np.sqrt(2 * np.log(2)) * sigmas[i, 1]  # error of FWHM
            resolutions[i, 0] = 100 * fwhms[i, 0] / means[i, 0]  # Resolution
            resolutions[i, 1] = resolutions[i, 0] * \
                                np.sqrt((fwhms[i, 1] / fwhms[i, 0]) ** 2 + (means[i, 1] / means[i, 0]) ** 2)  # error of R

        values = {'amplitude': amplitudes,  # the fit is returned in a dictionary
                      'mean': means,
                      'sigma': sigmas,
                      'FWHM': fwhms,
                      'R[%]': resolutions
                      }
        return values

    def plot(self, unit=0):
        energy_unit = ["keV", "MeV"]
        plt.figure(1)
        plt.step(self.xvals, self.yvals, color="tab:blue")
        plt.fill_between(self.xvals, self.yvals, color="tab:blue", step="pre")
        plt.ylim(0)
        plt.ylabel("Counts")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        if self.iscalibrated:
            plt.xlabel("Energy ("+energy_unit[unit]+")")
        else:
            plt.xlabel("ADC Channel")
        plt.show()
        return

    def plot_fit(self,unit=0):
        peaks = self.peaks()
        widths = peak_widths(self.yvals, peaks, rel_height=0.8)  # calculates the width of each peak (at 20% height)
        energy_unit = ["keV", "MeV"]
        plt.figure(1)
        plt.step(self.xvals, self.yvals, color="tab:blue")
        plt.fill_between(self.xvals, self.yvals, color="tab:blue", step="pre")
        plt.ylim(0)
        plt.ylabel("Counts")
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        if self.iscalibrated:
            plt.xlabel("Energy (" + energy_unit[unit] + ")")
        else:
            plt.xlabel("ADC Channel")
        fit = self.fit_peaks()
        for i in range(len(peaks)):
            fit_x = self.xvals[peaks[i] - int(widths[0][i]):peaks[i] + int(widths[0][i])]  # x fit range
            plt.plot(self.xvals[peaks[i]], self.yvals[peaks[i]], "kx")
            plt.plot(fit_x, Gaussian(fit_x, fit.get("amplitude")[i, 0], fit.get("mean")[i, 0], fit.get("sigma")[i, 0]),color="tab:red")

        plt.ylim(0)
        plt.show()
        return
