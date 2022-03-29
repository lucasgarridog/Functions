import uproot
import numpy as np

# This function loads the ROOT files created by COMPASS, both the event and histo files

def CompassLoader(filename, events=False):
    file = uproot.open(filename)                          # open ROOT file
    if events == True:
        data = file["Data_F"]["Energy"].array()           # for event files, saves the channel of each event
        return data                                       # so it needs to be plotted with plt.hist
                                                          # this contains all the channels, be careful!

    ch0 = file["Energy"]["_F_EnergyCH0@V1725S_646"]
    ch1 = file["Energy"]["_F_EnergyCH1@V1725S_646"]
    ch2 = file["Energy"]["_F_EnergyCH2@V1725S_646"]
    ch3 = file["Energy"]["_F_EnergyCH3@V1725S_646"]
    ch4 = file["Energy"]["_F_EnergyCH4@V1725S_646"]
    ch5 = file["Energy"]["_F_EnergyCH5@V1725S_646"]
    ch6 = file["Energy"]["_F_EnergyCH6@V1725S_646"]
    ch7 = file["Energy"]["_F_EnergyCH7@V1725S_646"]
    ch8 = file["Energy"]["_F_EnergyCH8@V1725S_646"]
    ch9 = file["Energy"]["_F_EnergyCH9@V1725S_646"]
    ch10 = file["Energy"]["_F_EnergyCH10@V1725S_646"]
    ch11 = file["Energy"]["_F_EnergyCH11@V1725S_646"]
    ch12 = file["Energy"]["_F_EnergyCH12@V1725S_646"]
    ch13 = file["Energy"]["_F_EnergyCH13@V1725S_646"]
    ch14 = file["Energy"]["_F_EnergyCH14@V1725S_646"]
    ch15 = file["Energy"]["_F_EnergyCH15@V1725S_646"]
    arr = np.zeros((2,len(ch0.values),16))                         # 2 rows (x,y), 16 channels
    arr[1,:,0] = np.array(ch0.values)                              # data is returned in an array
    arr[1,:,1] = np.array(ch1.values)
    arr[1,:,2] = np.array(ch2.values)
    arr[1,:,3] = np.array(ch3.values)
    arr[1,:,4] = np.array(ch4.values)
    arr[1,:,5] = np.array(ch5.values)
    arr[1,:,6] = np.array(ch6.values)
    arr[1,:,7] = np.array(ch7.values)
    arr[1,:,8] = np.array(ch8.values)
    arr[1,:,9] = np.array(ch9.values)
    arr[1,:,10] = np.array(ch10.values)
    arr[1,:,11] = np.array(ch11.values)
    arr[1,:,12] = np.array(ch12.values)
    arr[1,:,13] = np.array(ch13.values)
    arr[1,:,14] = np.array(ch14.values)
    arr[1,:,15] = np.array(ch15.values)
    arr[0,:,0] = np.array(range(1,len(arr[1,:,0])+1))
    arr[0,:,1] = np.array(range(1,len(arr[1,:,1])+1))
    arr[0,:,2] = np.array(range(1,len(arr[1,:,2])+1))
    arr[0,:,3] = np.array(range(1,len(arr[1,:,3])+1))
    arr[0,:,4] = np.array(range(1,len(arr[1,:,4])+1))
    arr[0,:,5] = np.array(range(1,len(arr[1,:,5])+1))
    arr[0,:,6] = np.array(range(1,len(arr[1,:,6])+1))
    arr[0,:,7] = np.array(range(1,len(arr[1,:,7])+1))
    arr[0,:,8] = np.array(range(1,len(arr[1,:,8])+1))
    arr[0,:,9] = np.array(range(1,len(arr[1,:,9])+1))
    arr[0,:,10] = np.array(range(1,len(arr[1,:,10])+1))
    arr[0,:,11] = np.array(range(1,len(arr[1,:,11])+1))
    arr[0,:,12] = np.array(range(1,len(arr[1,:,12])+1))
    arr[0,:,13] = np.array(range(1,len(arr[1,:,13])+1))
    arr[0,:,14] = np.array(range(1,len(arr[1,:,14])+1))
    arr[0,:,15] = np.array(range(1,len(arr[1,:,15])+1))
    return arr