import uproot
import matplotlib.pyplot as plt
import numpy as np

# CompassLoader loads the ROOT files created by COMPASS, both the event and histo files
# MvmeLoader loads the ROOT files created by COMPASS, both the event and histo files

def CompassLoader(filename, PSD=False, plot=False):
    file = uproot.open(filename)                          # open ROOT file
    # if events:
    #     E_data = file["Data_F"]["Energy"].array()           # for event files, saves the channel of each event
    #     if plot:                                          # so it needs to be plotted with plt.hist
    #         plt.figure(1)
    #         plt.hist(E_data, bins=200)                      # TO BE CONTINUED...
    #         plt.show()
    #     return E_data                                       # this contains all the channels, be careful!
    channels = np.arange(16)
    E_array = []
    for channel in channels:
        E_array.append(file["Energy"]["_F_EnergyCH" + str(channel) + "@V1725S_646"])

    E_data = np.zeros((2,len(E_array[0].values),16))                  # 2 rows (x,y), 16 channels
    for channel in channels:                                               # data is returned in an array
        E_data[1,:,channel] = np.array(E_array[channel].values)
        E_data[0,:,channel] = np.arange(len(E_data[1,:,channel]))

    PSD_data=None
    if PSD:
        PSD_array = []
        for channel in channels:
            PSD_array.append(file["PSD_E"]["_F_PSDvsECH" + str(channel) + "@V1725S_646"])
        PSD_data = []
        for channel in channels:
            data, bins = PSD_array[channel].numpy()
            data = data.T
            data = np.ma.masked_where(data == 0, data)                     # this converts 0 values to none
            PSD_data.append(data)

        if plot:
            plt.figure(channel)  # plotting
            plt.step(E_data[0,:,channel], E_data[1,:,channel], color="tab:blue", zorder=0)
            plt.fill_between(E_data[0,:,channel], E_data[1,:,channel], step="pre", color="tab:blue", zorder=0)
            plt.ylim(0)
            plt.title("Channel " + str(channel))
            plt.xlabel("ADC Channel")
            plt.ylabel("Counts")
    plt.show()

    return E_data, PSD_data

def MvmeLoader(filename):
    file = uproot.open(filename)
    channels = np.arange(16)
    E_array = []
    for channel in channels:
        E_array.append(file["MDPP16"]["MDPP16_" + str(channel)])
    E_data = np.zeros((2,len(E_array[0].values),16))
    for channel in channels:  # data is returned in an array
        E_data[1, :, channel] = np.array(E_array[channel].values)
        E_data[0, :, channel] = np.arange(len(E_data[1, :, channel]))

    return E_data
