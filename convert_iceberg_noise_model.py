import numpy as np
import matplotlib.pyplot as plt

file = np.load("ICEBERG_Noise_Model_For_Wirecell.npz")
nch, nfreq = file['FFT'].shape
nsmps = int(file['N'])

# # plot frequency spectrum for a channel
# fig, ax = plt.subplots()
# ax.plot(file['freq']/1e6, file['FFT'][200])
# 
# ax.set(xlabel='Frequency (MHz)', ylabel='FFT Amplitude (mV / MHz)')
# ax.grid()
# 
# plt.show()

import json
import bz2
 
# example of writing a compressed json
# data = [ {"groupID":4, "channels":[4]}, {"groupID":5, "channels":[5]}]
# with bz2.BZ2File("sample.json.bz2", "wb") as f:
#     f.write(json.dumps(data))

# groupID to channel map (one channel per group)
data = []
for ich in range(nch):
    data.append({"groupID": ich, "channels":[ich]})
with bz2.open("iceberg_group_to_channel_map_incoh.json.bz2", "wt", compresslevel=1) as f:
    f.write(json.dumps(data))

# noise model per group
data = []
for ich in range(nch): # one channel per group
    freqs = file['freq']*1e-9
    amps = file['FFT'][ich]*1e-9
    data.append({\
        "period": 500.0,\
        "nsamples": nsmps,\
        "const": 0.0,\
        "groupID": ich,\
        "freqs":freqs.tolist(),\
        "amps": amps.tolist()})
with bz2.open("iceberg_noise_model_incoh.json.bz2", "wt") as f:
    f.write(json.dumps(data))