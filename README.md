# Mparloc: Bayesian Method for Earthquake Location Using Multi-Parameter Data
(c) 2020-2021 Alessandro Caruso <alessandro.caruso@unina.it>,

[DOI](https://doi.org/10.1002/essoar.10503363.1)

## Introduction
I implemant a Bayesian method for the earthquake location which combines 
the information derived from the differential P-wave and S-wave arrival times, 
amplitude ratios and back-azimuths measured at a minimum of two stations. 
A complete posterior pdf of location is provided as a probability map on the 
computational grids and uncertainty of location is estimate as the range where 
the probability is over a prefixed threshold. 
In order to improve the speed performance, a pseudo-global search algorithm 
is implemented to explore in a fast way a even very dense computational grids.


## Config and Installing

Mparloc is written in pure python3 code tested with version 3.7.3, with several external library dependency:

    obspy, numpy, nllgrid, math, pyproj, taup, matplotlib

After the installation of all dependence you can test the library run in the main folder the command:

    python3 Example.py (Ubuntu environment) or python Example.py (Windows environment)

## Testing




## Important Note

This is a research implementation of the code, with no actual
real-time capabilities.
Feel free to pick the relevant parts of this (GPL licensed) code
and build your own real-time module!


## References

> Zollo, A. , Caruso, A. , De Landro , G. , Colombelli, S. , Elia, L. 2020.
> A Bayesian Method for Real-time Earthquake Location Using Multi-Parameter Data.
> Journal of geophysical research - solid earth doi: [https://doi.org/10.1002/essoar.10503363.1]



### DISCLAIMER
THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY.