import uproot
import matplotlib.pyplot as plt
import numpy as np

# This function loads the ROOT files created by COMPASS, both the event and histo files

def CompassLoader(filename, events=False, plot=False):
    file = uproot.open(filename)                          # open ROOT file
    if events:
        data = file["Data_F"]["Energy"].array()           # for event files, saves the channel of each event
        if plot:                                          # so it needs to be plotted with plt.hist
            plt.figure(1)
            plt.hist(data, bins=200)                      # TO BE CONTINUED...
            plt.show()
        return data                                       # this contains all the channels, be careful!

    channels = np.arange(16)
    channels_array = []
    for channel in channels:
        channels_array.append(file["Energy"]["_F_EnergyCH" + str(channel) + "@V1725S_646"])

    data = np.zeros((2,len(channels_array[0].values),16))                  # 2 rows (x,y), 16 channels
    for channel in channels:                                               # data is returned in an array
        data[1,:,channel] = np.array(channels_array[channel].values)
        data[0,:,channel] = np.arange(len(data[1,:,channel]))

        if plot:
            plt.figure(1)  # plotting
            plt.step(data[0,:,channel], data[1,:,channel], color="tab:blue", zorder=0)
            plt.fill_between(data[0,:,channel], data[1,:,channel], step="pre", color="tab:blue", zorder=0)
            plt.ylim(0)
            plt.title("Channel " + str(channel))
            plt.xlabel("ADC Channel")
            plt.ylabel("Counts")
            plt.show()

    return data