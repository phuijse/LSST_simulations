import numpy as np
import pickle
from os.path import join
from sys import argv
import time
import P4J

if len(argv) < 2:
    raise ValueError("Please run as python get_periods_AOV.py OBJ, where OBJ = RRab, RRc, CEPH or EB")
varstar = argv[1]

LSST_path = '/home/phuijse/Data/LSST/'
with open(join(LSST_path, 'templates', varstar+'_LSST.pkl'), "rb") as f:
    lc_data, lc_info, lc_per = pickle.load(f, encoding='latin1')
obj_id = sorted(list(lc_data.keys()))

perAOV = P4J.periodogram(method='MHAOV', debug=False)

"""
Multiband AoV implementation following 
https://github.com/Mondrik/Multiband_AoV_Demo/blob/master/Multiband_AoV.py

If using this code please cite:
N. Mondrik, et al. "A multiband generalization of the multiharmonic analysis of variance period 
estimation algorithm and the effect of inter-band observing cadence on period recovery rate." 
The Astrophysical Journal Letters vol. 811, n. 2, pp. 34,  2015
"""
Npoints = [12, 24, 36, 48]
filters = [b'u', b'g', b'r', b'i', b'z']
Nharmonics = [1, 2, 3]
res = np.zeros(shape=(len(obj_id), 5, len(Npoints), len(filters)+1, len(Nharmonics)))
start_time = time.time()
rng = np.random.RandomState(0)

for idx in range(len(obj_id)):
    print(idx)
    # Get a template
    data = lc_data[obj_id[idx]]
    for idx_r in range(res.shape[1]):
        # Draw Gaussian noise N(0, 1)
        white_noise = rng.randn(48, 5)
        for idx_n, N in enumerate(Npoints):
            for idx_p, Nh in enumerate(Nharmonics): 
                per_sum = 0.0
                d1 = float(2.*Nh)
                d2 = float(N - 2*Nh - 1)
                wvariance_acum = 0.0
                for idx_f, filt in enumerate(filters):
                    mjd = data[filt]['mjd'][:N]
                    mag = data[filt]['mag'][:N]                    
                    err = data[filt]['error'][:N]
                    # Add N(0, err**2) to the template mags
                    noisy_mag = mag + err*white_noise[:N, idx_f]                    
                    waverage = np.sum(noisy_mag/err**2)/np.sum(1./err**2)
                    wvariance = np.sum((noisy_mag - waverage)**2/err**2)
                    wvariance_acum += wvariance
                    perAOV.set_data(mjd, noisy_mag, err, whitten=True, Nharmonics=Nh)
                    perAOV.frequency_grid_evaluation(fmin=0.0, fmax=4.0, fresolution=1e-4)
                    freq, per = perAOV.get_periodogram()
                    per_sum += d1*per*wvariance/(d2 + d1*per)
                    res[idx, idx_r, idx_n, idx_f, idx_p] = 1.0/freq[np.argmax(per)]
                # Weighted average of the single band AoV periodograms
                per_sum = d2*per_sum/(d1*(wvariance_acum - per_sum))
                res[idx, idx_r, idx_n, idx_f+1, idx_p] = 1.0/freq[np.argmax(per_sum)]

pickle.dump(res, open("res_"+varstar+"_AOV.pkl", "wb"), protocol=2)
print("Elapsed time: %f [m]" %((time.time()-start_time)/60.0))
