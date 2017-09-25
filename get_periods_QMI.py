import numpy as np
import pickle
from os.path import join
from sys import argv
import time
import P4J

if len(argv) < 4:
    raise ValueError("Please run as python get_periods_QMI.py OBJ hm hp, where OBJ = RRab, RRc, CEPH or EB, and hm/hp are the kernel sizes")
varstar = argv[1]

with open(join('templates', varstar+'_LSST.pkl'), "rb") as f:
    lc_data, lc_info, lc_per = pickle.load(f, encoding='latin1')
obj_id = sorted(list(lc_data.keys()))

perEU = P4J.periodogram(method='QMIEU', debug=False)
perCS = P4J.periodogram(method='QMICS', debug=False)

hmm, hp = 1., 1. 
hmm, hp = float(argv[2]), float(argv[3])

Npoints = [12, 24, 36, 48]
filters = [b'u', b'g', b'r', b'i', b'z']
Nrealizations = 5
methods = [perEU, perCS]
res = np.zeros(shape=(len(obj_id), Nrealizations, len(Npoints), len(filters)+1, len(methods)))
start_time = time.time()
rng = np.random.RandomState(0)

for idx in range(len(obj_id)):
    print(idx)
    # Get a template
    data = lc_data[obj_id[idx]]
    for idx_r in range(Nrealizations):
        # Draw Gaussian noise N(0, 1)
        white_noise = rng.randn(Npoints[-1], len(filters))
        for idx_n, N in enumerate(Npoints):
            for idx_p, method in enumerate(methods):
                per_sum = 0.0
                for idx_f, filt in enumerate(filters):
                    mjd = data[filt]['mjd'][:N]
                    mag = data[filt]['mag'][:N]                    
                    err = data[filt]['error'][:N]
                    # Add N(0, err**2) to the template mags
                    noisy_mag = mag + err*white_noise[:N, idx_f]                    
                    method.set_data(mjd, noisy_mag, err, whitten=False, h_KDE_M=hmm, h_KDE_P=hp)
                    method.frequency_grid_evaluation(fmin=0.0, fmax=4.0, fresolution=1e-4)
                    #method.finetune_best_frequencies(fresolution=2e-5)
                    freq, per = method.get_periodogram()
                    per_sum += per
                    res[idx, idx_r, idx_n, idx_f, idx_p] = 1.0/freq[np.argmax(per)]          
                res[idx, idx_r, idx_n, idx_f+1, idx_p] = 1.0/freq[np.argmax(per_sum)]

pickle.dump(res, open("res_QMI_"+varstar+"_"+str(hmm)+"_"+str(hp)+".pkl", "wb"), protocol=2)
print("Elapsed time: %f [m]" %((time.time()-start_time)/60.0))
