import csv
import numpy as np


def OsciReader(fname):
    sample_period = 0.0
    raw_samples = []

    with open(fname, 'r') as csvfile:
        c = csv.reader(csvfile)

        # Sample period is in cell B2 (1,1)

        for row_num, row in enumerate(c):
            if row_num == 1: # get the sample period
                sample_period = float(row[1])
                break

        # Sample data starts after the last header line
        # containing the firmware version.
        in_header = True
        for row in c:
            if in_header:
                if row[0] == 'Firmware Version':
                    in_header = False
            else:
                raw_samples.append(float(row[4]))

    y = np.array(raw_samples)
    x = np.arange(0, len(y)*sample_period, sample_period)
    x = x*1e9                                                 # x in nanoseconds

    return x, y, sample_period