from Fits import *

# This function calibrates the imput spectra with the given energies

def SpectraCalibrator(x, y, energies, plot=False):
    peaks, _ = find_peaks(y, prominence=300)                                    # find the gaussian peaks
    if len(energies) != len(peaks):                                             # terminate the program if number of energies
        print("Error. \# of peaks does not match \# of energies.")              # given does not match number of peaks
        return

    peaks_channel = [x[i] for i in peaks]                                       # channel where the peaks are
    fit = Linear_fit(peaks_channel, energies)                                   # performs linear calibration
    calibrated_x = [fit.get("slope") * i + fit.get("intercept") for i in x]     # calibrate x axis
    if plot == True:
        plt.figure(1)                                                           # plotting
        plt.step(calibrated_x, y, color="tab:blue", zorder=0)
        plt.fill_between(calibrated_x, y, step="pre", color="tab:blue", zorder=0)
        plt.ylim(0)

    return fit, calibrated_x                                                    # fit parameters are returned in a dictionary