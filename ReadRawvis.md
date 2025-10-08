**read_rawvis**



The *read_rawvis.py* python script reads  **SPOTLIGHT** 1.3ms raw visibility files, extracts (and averages, if needed) selected data and make a plots. This program is mainly meant for debugging use, and may not cater to all possible situations that may arise during **SPOTLIGHT** observations. **SPOTLIGHT** raw visibilities are time multiplexed over 16 files, with the time multiplexing being done over a set of 50 records. Each 1.3ms record cosists  of 4096 channels for 1056 baselines. Visibilities are complex, half float, and the script assumes that the sampler to antenna connection is the default one, with the last two data streams corresponding to dummy antennas. 



The main modes are :

- print the polarization averaged data for some antenna pair (*--ants*)

- print the polarization averaged data for all selfs (*--self*)

- print the polarization averaged data for all cross (default)



The switches that control the action of the script are detailed below.

- *-a (--ants)* a0 a1 . Plot only the amplitude data for the baseline between antennas a0 and a1 (both stokes parameters are averaged together). By default it will plot the average of all cross baselines.

- *-b (--bandpass)*  Apply an amplitude bandpass correction. The "bandpass"" is computed as the median amplitude over all records in the timeslice.

- *-B (--baseline)* subtract the median amplitude (over all records in the time slice).

- *-C (--collate)* collate all the data over the specified files and timeslices into a single plot. By default a separate plot is made for each time slice. This is useful, when for e.g. a given burst dynamic spectrum is spread over multiple slices.

- *-d (--drop_csq)* drop the central square baselines from the set of baselines being averaged together. The "central square baselines" are defined as those between antennas with antenna number < 11. In any case a set of the shortest baselines are always dropped from the average (even if --drop_csq is not set). This is ignored if a specific baseline is requested via the --ants option.

- *-f (--files*) space separated list of 1.3ms raw visibility files to process.

- *-I (--interactive*) plot is displayed on an interactive matplotlib window. By default the plot is saved in a pdf file.

- *-n (--nan_stats)* plot the nan statistics.  nan's arise typically when the visibility values exceed what can be contained in a half float. (Mostly this appears to happen only in the self correlations)

- *-o (--overplot) *svfits.log overplot the channel selection made in svfits when identifying the channels containing the burst signal. The channel selection is read from the file name given as the  argument, which should be  the log file created by svfits (default name svfits.log).

- *-p (--power)* add the visibility powers instead of making a phase coherent sum of the visibilities. This may be of some use when looking for bright off centre bursts.

- *-s (--slice)* if collation is not being done, then this is a set of two numbers giving the start and stop timeslice number (the same for all files). If the plots are being collated then it is an array of slice numbers to use, one for each input file.

- *-S (--self) * plot the sum of all the self correlations.






