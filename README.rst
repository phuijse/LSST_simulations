LSST simulations
================

This repository holds codes associated to:

Pablo Huijse et al., "Robust period estimation using mutual information for multi-band light curves in the synoptic survey era", to be published in the ApJ Supplement Series "Special Issue on Stellar Astronomy Big Data", 2017, https://arxiv.org/abs/1709.03541

The pickled data in the template folder contain light curves generated using the OpSim and CatSim tools from the LSST simulations framework. See:

https://github.com/phuijse/LSST_simulations/blob/master/Example.ipynb 

for a example on how to draw an ab-type RR Lyrae light curve using these templates.

The get_periods scripts were used in the paper mentioned above to compare period estimation accuracy between different methods. To run these you will need to install the https://github.com/phuijse/P4J/ library

Generating you own LSST-like templates using OpSim and CatSim
-------------------------------------------------------------

The template files were built using the create_obj script. If you want to run the script you will need to install the LSST stack. To install the LSST stack first setup the core packages by following steps 1 and 2 (Installing from Source) at 

https://confluence.lsstcorp.org/display/SIM/Catalogs+and+MAF

After that you will need to download the **Reference simulated survey** from 

https://www.lsst.org/scientists/simulations/opsim/opsim-survey-data 

Download and extract the minion_2016 SQLite DB and point the opsimdb variable in create_obj to it.

To run the script::

    cd /my_lsst_stack_dir
    source loadLSST.bash
    setup sims_catUtils -t sims
    create_obj RRL

For more details please see 

http://confluence.lsstcorp.org/display/SIM/Catalog+Simulations+Documentation. 




