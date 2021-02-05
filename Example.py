#!/usr/bin/python3
#IMPORT OF LIBRERY
from MParLoc import Launcher 
from MParLoc import grid_creator
##############################
# RUN THE GRID CREATOR CODE THAT READ THE NONLINLOC GRIDS 
#AND CREATE A DIFFERENTIAL AMPLITUDE AND BACK-AZIMUTH GRIDS
grid_creator('EXAMPLE','sac') 
#RUN MPARLOC LOCATION
Launcher('EXAMPLE','sac')
