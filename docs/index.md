#Simple Circulation Based Downscaling code

This code is designed to provide a simple statistical downscaling method relying on 
regressions and statistical transformations from various inputs (e.g. precipitation, 
humidity, wind, PCA, etc.) to various outputs (e.g. precipitation, temperature, etc.)

##Useful commands
Use the following to generate a list of e.g. GEFS precipitation files for input. 
    ls -1 gefs/2010/*/apcp_sfc_*_mean.nc | sed 's/*//g;s/$/"/g;s/^/"/g'>gefs_pr_file.txt


##Requirements
T.B.D

##Developing
T.B.D

##Reference
T.B.D.
