This folder contains the raw data from the online Moser database for 
retau 180, 395, 590 (saved as .txt files).
1. Run the ConvertPlotMoserDNS_june.m to read in the raw (half channel) data from .txt files
and convert to full channel data set (through mirroring). 

2. Run the interpolateWNgrid_june.m to interpolate the data from chebychev WN grid to 
hyperbolic tangent stretced grid using in ARNL LESGO simulations. 
This saves the input data files used for the resolvent codes from the Moser mean profiles.
