# Mparloc: Bayesian Method for Earthquake Location Using Multi-Parameter Data
(c) 2020-2021 Alessandro Caruso <alessandro.caruso@unina.it>,


Post-doctoral researcher at RISSCLab - Department of Physics (University of Naples Federico II)

[DOI](https://doi.org/10.1029/2020JB020359)

## Introduction
I implement a Bayesian method for the earthquake location, which combines 
the information derived from the differential P-wave and S-wave arrival times, 
amplitude ratios and back-azimuths measured at a minimum of two stations. 
A complete posterior pdf of location is provided as a probability map
and uncertainty of location is estimate as the range where 
the probability is over a prefixed threshold. 
In order to improve the speed performance, a pseudo-global search algorithm is implemented
to explore in a fast way an even very dense computational grids.


## Config, Installing and running

Mparloc is written in pure python3 code tested with version 3.7.3 and several external library dependency installed:

    obspy, numpy, nllgrid, math, pyproj, taup, matplotlib
	
Mparloc library consists of two main functions implemented in the source file `Mparloc.py` :
 
    grid_creatore(Name_network,Name_event) # Tool for creating grids starting from NnLinLoc grids.
    Launcher(Name_network,Name_event) # Function that initialize and launch the location algorithm.

`Name_network` and `Name_event` are string object that contain respectively the name of the network folder 
and the name of records folder. The network folder contains: a sub-folder `grids`, a records folder and 
a text file named `list_stations_4loc.dat` containing the station name list used by Mparloc.
The `grids` folder contains `.buf` and `.hdr` files of time grids created with tools of Nonlinloc,
while in the records folder there are the accelerometric `sac` files with their station 
fields and transduction coefficient correctly fixed in header.
Following an example of `list_stations_4loc.dat`:

    AND3
    AVG3
    BEL3
    BENI
    BSC3
    CGG3 ....

The configuration of module is done by compiling two files located in the `config_files` folder:

    run_loc.conf # Main configuration file
    weights_table.dat # Correspondence table between SNR(Signal to Noise Ratio) and the error in second on seismic pick phase
	
Example of `weights_table` file:

    # SNRmin SNRmax Wp(error on P) Ws(error on S)
    10 0.25 0.55  # In the first line there is only SNRmax
    10 20 0.15 0.25
    20 50 0.1 0.15
    50 0.015 0.015 # In the last line there is only SNRmin

Example of `run_loc` file:

    #Setting parameters:
	
    #Parameter to set grid computation precision [f for single-precision, d for double-precision]
    grid_precision: f
    #Confidence level for error estimate [float number between 0 and 1]
    Confidence_probability_threshold_level: 0.5
    #Time window in seconds for SNR estimate for P phase [float positive number]
    snr_wind_p: 0.5
    #Time window in seconds for SNR estimate for S phase [float positive number]
    snr_wind_s: 0.5
    #Minimum number of P phases for eqk location [integer positive number >2]
    no_min_p: 2
    #Minimum number of S phases for eqk location [integer positive number >1]
    no_min_s: 1
    #Use back-azimuth data for eqk location [False or True]
    back_az: True
    #back-azimuth standard deviation in degree [integer positive number]
    back_az_error: 30
    #Use Pv differential amplitude data for eqk location [False or True]
    diff_Amp: True
    #Phase Amplitude attenuation coefficient [float negative number]
    diff_Amp_b: -1.78
    #differential amplitude standard deviation [float positive number]
    diff_Amp_error: 0.5
    #Use S phases for eqk location [False or True]
    s_phases: True

After configuring the software and correctly formatting the input files, it is possible to execute
the routine to perform the construction of the grids. For example, naming the network folder as `ISNet`
and the event folder as `12458a`:

    #!/usr/bin/python3
    #IMPORT OF LIBRERY
    from MParLoc import grid_creator
	# RUN THE GRID CREATOR CODE THAT READ THE NONLINLOC GRIDS 
	#AND CREATE DIFFERENTIAL AMPLITUDE AND BACK-AZIMUTH GRIDS
    grid_creator('ISNet','12458a')

Finally you can run the location algorithm:

    #!/usr/bin/python3
    #IMPORT OF LIBRERY
    from MParLoc import Launcher
	# RUN THE LOCATION ALGORITHM
    Launcher('ISNet','12458a')

After running, the code provide in the folder of event `12458a` an hypo71 format file `MParLoc_event.h71` with all
available information about the event location, the location Map `MParLoc_event_FIGURE.pdf` which gives the shape
of the spatial trend of the probability and a textual log `MParLoc_event.log` of main operations performed by the algorithm.



## Testing


After the installation of all dependence, you can test the library run in the main folder the command:

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
