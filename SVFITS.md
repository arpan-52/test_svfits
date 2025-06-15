# SVFITS

## Overview

svfits is a C utility  for conversion of the SPOTLIGHT raw 1.3ms visibility dump files into a multi-source random group UVFITS file. This can then be imported into standard interferometric imaging packages. The default mode is to identify the time-frequency subset of the data that contains the burst signal and convert only that data into a UVFITS file. This results in  the output file being significantly smaller than the input file

## Background

The [SPOTLIGHT](https://spotlight.ncra.tifr.res.in/) system at the GMRT can be triggered to dump the visibilities for candidate transients. The visibilities have an integration time of 1.31072 ms and 4096 spectral channels and are contained in 16 files which are time multiplexed.  The total time interval of the data that is dumped depends on the candidate's dispersion measure, but the maximum duration is expected to be about 30s. The visibility files contain a timestamp, but not metadata. The metadata has to be read in from a separate file. 

The svfits utility reads the 16 time multiplexed data files as well as the associated metadata file and produces a random group UVFITS file which can be used by the imaging and localization pipeline.

The default mode is to write out only those visibility timestamps and frequency channels that contain the burst signal. The time of arrival and DM of the burst is used to identify the time-frequency region that contains the burst signal. These visibilities are then converted from half-float to float, and written out as a single channel UVFITS file. The (u,v,w) coordinates are corrected to a standard reference frequency, and are also rotated to J2000. There are options for using the off source data near the location of the burst signal to do an "amplitude bandpass" (i.e. flatten the band response) as well as subtract off the "baseline" (as well as RFI). Use of these two options generally significantly improves the noise properties.

Options to write out all of the channels for the records containing the burst signal, as well as to convert all of the dumped visibilities (i.e. irrespective of whether they contain the burst signal or not) to UVFITS are also provided. The last option would be useful for example in making an image containing bright background sources that could be used for astrometry.

Parts of the code that do the actual conversion to random group UVFITS are modified versions of the code used for gvfits.

## Command Line options

***-a*** File containing the metadata regarding the SPOTLIGHT correlator configuration parameters, sampler connections, antenna coordinates etc. (The antsamp.hdr file) It is important to confirm that the parameters in this file conform to those actually used in the observation. At the moment this confirmation has to be done manually, which is a stop gap arrangement. In future one would expect these parameters to be taken from the "corr"  data structures used by the SPOTLIGHT online software, and which should be written to disk along with the visibility files.

***-b*** File containing the binary headers dumped by the DAS chain (i.e. the corr and the scan structures). The correlator set up and and observational parameters read from this file over ride the built in defaults, and are in turn over ridden by user supplied values provided by the svfits_par.txt file. The final intention is to have a minimum set of paramters read from the svfits_par.txt file (ideally only the burst related paramters and the parameters that control conversion.) **As of now however, it appears that the SPOTLIGHT binary headers are incomplete - although the -b option is available it is not advised to use it until the issue of unfilled parameters in the binary headers is resolved.**

***-u*** file containing the user parameters, (default svfits_par.txt), these include the observational parameters such as the phase centre coordinates, the frequency settings etc., as well as the candidate burst parameters returned by the trigger software, as well as various processing options chosen by the user. In future one would expect the observational parameters to be taken from the "scan" structure used by the SPOTLIGHT online software, which should be written to disk along with the visibility files.

## User Supplied Parameters

The general format of the user parameter file (i.e. the one specified using the -u command line option; the default is svfits_par.txt) is a set of KEYWORD VALUE pairs, which each KEYWORD and corresponding VALUE appearing on a separate line. Lines starting with an asterisk '*' indicate a comment, and are not parsed. Similarly any text following a bang '!' after the VALUE are regarded as a comment and are not parsed. For boolean parameters (e.g. **ALL_CHAN, ALL_DATA**,  0 means to not set the corresponding option, while 1 means it is set). The currently available parameters are listed below.

- **NFILE** The number of 1.3ms raw visibility files to process. All time slice and other calculations are made assuming that there are a total of 16 files. This entry only specifies how many of those files are to be read for conversion to FITS.

- **PATH**  The directory containing the raw visibility files. It is assumed that all files are in the same directory.

- **INPUT**  The names of the raw visibility files. The total number of files listed here needs to match the number given in NFILE. Currently all files have identical timestamps, so **the order here is important**. It is assumed that the first file in this list contains the first time slice, the second file contains the second timeslice and so on. In future it is expected that each file will contain enough information to determine which timeslice it contains.

- **HAVE_IDX** the first set of visibility dump files did not include the file index in the meta data. That meant that the file index (i.e. which in the set of 16 slices the data in this file corresponds to) had to be supplied by hand. Later (i.e. from June 2025) versions of the raw visibility files include the index in the meta-data. HAVE_IDX=0 is meant for processing historical data where the meta-data does not contain the index. The default is to assume that the file meta-data contains the index.

- **N_SLICE** the number of slices of data to process from each file. This is relevant only when the **ALL_DATA** flag is set. The default is to process all slices.

- **FITS** The name of the output FITS file. By default it is "TEST.FITS". **Existing files will be overwritten!**

- **ANTMASK**  Mask specifying which antenna's data to use. Baselines containing antennas which are not part of ANTMASK will be dropped. There is only one ANTMASK which is applied to both polarizations. If Bit n of the ANTMASK is 1 then antenna n (in the antenna list contained in the correlator parameters, currently the antsamp.hdr file) will be **included**.

- **ALL_CHAN** copy all channels for records containing the burst signal to the output FITS file. In this case the output FITS file will have 4096 spectral channels. By default the output FITS file has only one channel. 

- **ALL_DATA** convert all of the data into the UVFITS format, irrespective of the burst parameters. This could be useful to make a continuum map to identify background point sources that could be used for astrometry. To reduce the data size, all of the records in a given slice are averaged together in time. Averaging in frequency is also possible by setting NCHAV to > 1.

- **NCHAV** The number of channels to average together when creating the output UVFITS file. This applies only when the **ALL_DATA** flag is set.

- IATUTC offset in seconds between IAT and UTC. The system has a built in default which can be over ridden here.

- **EPOCH**  No longer used. The output visibilities are forced to J2000

- **OBS_MJD** the start MJD of the observation. 

- **FREQ_SET**  The frequency setting, i.e. a colon separated frequency of the first, last channel (Hz) and the number of channels. For example 5.5e6:7.5e8:4096 for Band 4. In future this should be taken from the corr and scan structures.

- **COORD_TYPE** No longer used

- **RECENTRE** change the phase centre to the BURST_{RA,DEC}.

- **RA_APP** The apparent right ascension of the phase centre. In future this should be taken from the scan structure

- **DEC_APP** The apparent declination of the phase centre. In future this should be taken from the scan structure.

- **RA_MEAN** The mean right ascension of the phase centre. In future this should be taken from the scan structure

- **DEC_MEAN** The mean declination of the phase centre. In future this should be taken from the scan structure.   

- **STOKES_TYPE** Label to be given to the stokes parameters, i.e. "XY" or "RL"

- **BURST_NAME** Identifier for the burst

- **BURST_MJD** The start time of the burst, as reported by the trigger software

- **BURST_TIME** The start time (sec) of the burst. This is in general computed from the MJD, (if **UPDATE_BURST** is set), however it can also be given by hand here. In general this is useful only for debugging.

- **BURST_DT** the error in the burst time. Currently unused.

- **BURST_INTWD** the intrinsic  width of the burst (sec) as reported by the trigger software

- **BURST_WIDTH**  The total duration of the burst signal (including the DM sweep). In general this is computed from the other parameters when UPDATE_BURST is set. In general useful only for debugging.

- **BURST_DM** The dispersion measure of the burst as reported by the trigger software.

- **BURST_DDM** The error in the DM, currently unused.

- **BURST_FREQ** The frequency that the BURST_MJD (or BURST_TIME) refers to. (Hz)

- **BURST_BM_ID**  The id (i.e. number) of the beam in which the burst was detected

- **BURST_RA** The right ascension of the beam in which the burst was detected. One could rotate the visibilities to this coordinate before imaging to minimize the w-term error.

- **BURST_DEC** The declination of the beam in which the burst was detected. One could rotate the visibilities to this coordinate before imaging to minimize the w-term error.

- **UPDATE_BURST** update the remaining burst parameters using the BURST_MJD, BURST_INTWD and BURSt_DM

- **DO_FLAG** Do MAD based flagging of the output visibilities. For burst data the MAD is computed only using the visibilities which contain the burst signal. When averaging over channels (i.e. when **ALL_DATA** is set) the MAD is computed separately for each baseline and time slice. Unused otherwise (i.e. unused when **ALL_CHAN** is set)

- **THRESH** The threshold for MAD based flagging. All visibilities whose deviation from the media is greater than THRESH times the MAD are flagged.

- **DO_BAND** use the off source (i.e. visibilities without the burst signal) data to compute the "amplitude bandpass", in practice this is just the mean of the off source visibility amplitudes (computed separately per baseline per timeslice). Unused when **ALL_CHAN** and **ALL_DATA** are set. The bandpass is normalized, so the overall scaling is not affected.

- **DO_BASE** subtract off the off burst. This does a reasonable job of removing RFI. In general the noise properties of the data with DO_BASE and DO_BAND applied are significantly better than the original.

- **POST_CORR** compute the post-correlation beam. Currently this is computed at the phase centre , but in future the option of rotating to the beam co-ordinates before making the beam will be added. This makes for a useful diagnostic plot (see plot_postcorr.py)

- **DROP_CSQ** (drop the visibilities between central square antennas, (actually antennas with numbers less than 11) when computing the post-correlation beam.

- **NUM_THREADS** number of threads to use when parallelizing the more compute intensive parts of the code.
