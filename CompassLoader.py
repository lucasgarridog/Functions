import uproot
import numpy as np

def CompassLoader(filename, events=False):
    file = uproot.open(filename)                          # open ROOT file
    if events == True:
        data = file["Data_F"]["Energy"].array()           # for event files, saves the channel of each event
        return data                                       # so it needs to be plotted with plt.hist

    data = file["Energy"]["_F_EnergyCH0@V1725S_646"]      # for histo files, its only reads channel 0
    y = np.array(data.values)                             # edit for more channels
    x = np.array(range(1,len(y)+1))
    arr = np.zeros((2,len(y)))
    arr[0,:] = x                                          # data is returned in an array
    arr[1,:] = y
    return arr