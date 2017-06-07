from __future__ import with_statement, division
import os
import numpy as np
import json
import re
from sys import argv
import pickle

if len(argv) < 2:
    raise ValueError("Please run as python create_obj.py OBJ, where OBJ = RRab, RRc, CEPH or EB")

OBJ = argv[1]
if OBJ == 'CEPH':
    vmethod = 'applyCepheid'
    vpath = 'cepheid_lc'
elif OBJ == 'RRab':
    vmethod = 'applyRRly'
    vpath = 'rrly_lc/RRab'
elif OBJ == 'RRc':
    vmethod = 'applyRRly'
    vpath = 'rrly_lc/RRc'
elif OBJ == 'EB':
    vmethod = 'applyEb'
    vpath = 'eb_lc'
else:
    raise ValueError("Wrong object identifier, try RRab, RRc, CEPH or EB")

# Set the location of minion_1016 database
opsimdb = os.path.join("/home/phuijse/Data", "LSST")
opsimdb = os.path.join(opsimdb, "minion_1016_sqlite.db")

from lsst.sims.catUtils.utils import StellarLightCurveGenerator
from lsst.sims.catalogs.db import fileDBObject
from lsst.utils import getPackageDir
library_dir = getPackageDir('sims_sed_library')

# How many light curves to create:
n_stars = 1000
rng = np.random.RandomState(0)

# Stars will be placed randomly in RA \in [65,75], DEC \in [-30, -20]
RA_list = rng.random_sample(n_stars)*10.0 + 65.0
DEC_list = rng.random_sample(n_stars)*10.0 - 30.0

# Get Spectral Energy Distributions (SEDs) and normalizing magnitudes from main sequence stars
possible_seds = sorted(os.listdir(os.path.join(library_dir, 'starSED', 'kurucz')))
seds = []
for i in rng.randint(0, len(possible_seds), n_stars):
    seds.append(possible_seds[i])
norm_mag = rng.random_sample(n_stars)*5.0 + 20.0

# Random initial MJD 
t0_list = 59580.0 - rng.random_sample(n_stars)*1000.0

# Get Cepheid templates and periods
possible_lc = sorted(os.listdir(os.path.join(library_dir, vpath)))
print("%d templates found" %(len(possible_lc)))
possible_periods = np.zeros(shape=(len(possible_lc),))
re_match = re.compile('\d+\.\d+')
for i, lc in enumerate(possible_lc):
    with open(os.path.join(library_dir, vpath, lc), 'r') as fid:
        per_line = fid.readlines()[1]
        possible_periods[i] = float(re_match.findall(per_line)[0])
        # print(possible_periods[i])
lc_files = []
periods = np.zeros(shape=(n_stars,))
for j, i in enumerate(rng.randint(0, len(possible_lc), n_stars)):
    lc_files.append(os.path.join(vpath, possible_lc[i]))
    periods[j] = possible_periods[i]

# Prepare param dictionary
param_list = []
for i in range(n_stars):
    if OBJ == 'CEPH' or OBJ == 'EB':
        var_dict = {'varMethodName': vmethod, 
                'pars':{'t0':t0_list[i], 'lcfile': lc_files[i], 'period':periods[i]}}
    elif OBJ == 'RRab' or OBJ == 'RRc':  # pars field names are different
         var_dict = {'varMethodName': vmethod, 
                 'pars':{'tStartMjd':t0_list[i], 'filename': lc_files[i], 'period': periods[i]}}
    param_list.append(json.dumps(var_dict))

# Create catalog file
catalog_file = 'tmp_cat.txt'
with open(catalog_file, 'w') as fid:
    fid.write("# lore ipsum\n")
    for i in range(n_stars):
        fid.write("%d; %f; %f; %s; %f; %s; 0.1; 0.0\n" % (i+1, RA_list[i], DEC_list[i],
            seds[i], norm_mag[i], param_list[i]))

# Create object DB
class genObjectDB(fileDBObject):
    objid = 'dummyVariableStar'
    skipRegistration = True
    raColName = 'raDeg'
    decColName = 'decDeg'
    objectTypeId = 0
    columns = [('raJ2000', 'raDeg*PI()/180.0', np.float),
            ('decJ2000', 'decDeg*PI()/180.0', np.float),
            ('ebv', 'galacticAv/3.1', np.float)]

periodic_variable_db = genObjectDB(catalog_file, idColKey='id', delimiter=';',
        dtype = np.dtype([('id', int), ('raDeg', float), ('decDeg', float),
            ('sedFilename', str, 300), ('magNorm', float), ('varParamStr', str, 300),
            ('galacticAv', float), ('parallax', float)]))

# Query DB to get the light curves
RA = (64.0, 76.0)
DEC = (-31.0, -19.0)
gen = StellarLightCurveGenerator(periodic_variable_db, opsimdb)
pointings = gen.get_pointings(RA, DEC, bandpass=('u', 'g', 'r', 'i', 'z'))
lc_dict, lc_info = gen.light_curves_from_pointings(pointings)

# Get periods and write all in a pkl file
lc_period = {}
for i, id_lc in enumerate(lc_info.keys()):
    lc_period[id_lc] = float(json.loads(lc_info[id_lc])["pars"]["period"])

pickle.dump([lc_dict, lc_info, lc_period], open(OBJ+"_LSST.pkl", "wb"), protocol=2)

