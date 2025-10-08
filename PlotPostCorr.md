**plot_postcorr**

The plot_postcorr.py script reads  the postcorrelation beam datafile created by svfits and plots it, along with the data selection done by svfits, as recorded in the svfits.log file created while running svfits.   The dynamic spectrum is plotted along with the computed start and end time of the burst in each channel. Only the data lying between these two cuts is converted to FITS in the default mode. The program also creates a de-dispersed profile, and plots it along with the burst time reported by the trigger software (i.e. the time and DM used by svfits to determine the location of the burst signal in the time-frequency plane). 

The following options are available:

-d (--dm) The dispersion measure to use when de-dispersing the signal. By default use the one stored in the header of the post correlation beam file.

-f (--freq) The start and stop frequencies (GHz) of the observation. By default use the values stored in the header of the post correlation beam file.

-i (--integ) The integration time of the visibilities (sec). By default use the values stored in the header of the post correlation beam file.

-l (--log_file) The svfits log file from which to read the locations of the burst signal. By default use svfits.log

-m (--mjd) The MJD of the burst. By default use the value stored in the header of the post correlation beam file.

-p (--post_corr) The name of the post correlation beam file. By default use svfits_postcorr.dat


