# Mparloc: Bayesian Method for Earthquake Location Using Multi-Parameter Data
(c) 2020-2021 Alessandro Caruso <alessandro.caruso@unina.it>,

[DOI](https://doi.org/10.1029/2020JB020359)

## Introduction
I implemant a Bayesian method for the earthquake location which combines 
the information derived from the differential P-wave and S-wave arrival times, 
amplitude ratios and back-azimuths measured at a minimum of two stations. 
A complete posterior pdf of location is provided as a probability map on the 
computational grids and uncertainty of location is estimate as the range where 
the probability is over a prefixed threshold. 
In order to improve the speed performance, a pseudo-global search (like octree)
algorithm is implemented to explore in a fast way a even very dense computational grids.


## Config and Installing

Mparloc is written in pure python3 code tested with version 3.7.3, with several external library dependency:

    obspy, numpy, nllgrid, math, pyproj, taup, matplotlib
	
Mparloc library consists of two main functions implemented in the source file Mparloc.py:
 
    grid_creatore(Name_network,Name_event) # Tool for creating grids starting from NnLinLoc grids.
    Launcher(Name_network,Name_event) # Function that initialize and launch the Location algorithm.

`Name_network` is a string object that contain the name of the folder where are:

The `.buf` and `.hdr` files of time grids created with Nonlinloc tools for the 
specific network and located in a dedicated sub-folder named `grids`.

The `sac` files of the accelerometric data contained in a sub-folder whose name is 
contained in the string `Name_event` with the station fields and transduction 
coefficient correctly set in the header.

The text file `list_stations_4loc.dat` with the station name list used from 
location algorithm, such as:

    AND3
    AVG3
    BEL3
    BENI
    BSC3
    CGG3 ....

The configuration of the module is done by compiling two files located in the `config_files` folder:

    run_loc.conf # Main configuration file
    weights_table.dat # Correspondence table between SNR(Signal to Noise Ratio) and the error in second on seismic pick phase
	
Example of `weights_table` file:

    # SNRmin SNRmax Wp Ws
    10 0.25 0.55  # In the first line there is only SNRmax
    10 20 0.15 0.25
    20 50 0.1 0.15
    50 0.015 0.015 # In the last line there is only SNRmin

Example of `run_loc` file:

    #Setting parameters:
    grid_precision: f
    #Confidence level for error estimate (in percentage)
    Confidence_probability_threshold_level: 0.5
    #Time window in seconds for SNR estimate for P phase
    snr_wind_p: 0.5
    #Time window in seconds for SNR estimate for S phase
    snr_wind_s: 0.5
    #Minimum number of P phases for eqk location
    no_min_p: 2
    #Minimum number of S phases for eqk location
    no_min_s: 2
    #Use back-azimuth data for eqk location
    back_az: True
    #back-azimuth standard deviation in degree
    back_az_error: 30
    #Use Pv differential amplitude data for eqk location
    diff_Amp: True
    #differential amplitude b coefficient:
    diff_Amp_b: -1.78
    #differential amplitude standard deviation
    diff_Amp_error: 0.5
    #Use S phases for eqk location
    s_phases: True


After the installation of all dependence you can test the library run in the main folder the command:

    python3 Example.py (Ubuntu environment)
	python Example.py (Windows environment)
	
if the code work well the output files are provided 
at path: @localpath/Mparloc/EXAMPLE/sac/

## Important Note

This is a research implementation of the code, with no actual
real-time capabilities.
Feel free to pick the relevant parts of this (GPL licensed) code
and build your own real-time module!


## References

> Zollo, A. , Caruso, A. , De Landro , G. , Colombelli, S. , Elia, L. 2020.
> A Bayesian Method for Real-time Earthquake Location Using Multi-Parameter Data.
> JOURNAL OF GEOPHYSICAL RESEARCH: SOLID EARTH doi: [https://doi.org/10.1029/2020JB020359]



### DISCLAIMER
THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY.