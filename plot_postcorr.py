# script to plot the postcorrelation beam made by svfits, overlay the selection
# done by svfits etc. Seems to be working for both Band 3 (frequency decreases
# with channel number) and Band 4 (frequency increases with channel number).
#
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticky
from   scipy import stats
from   scipy.ndimage import shift
from   astropy.time import Time
import copy

import argparse

K=4.15e-3 # constant in dispersion calculation

def dedisperse(dyn_spc,dt,fs,fe,DM):
    """ Does incoherent dedispersion."""
    nrec,channels=dyn_spc.shape;
    fm=0.5*(fs+fe)
    df=(fe-fs)/channels
    prof=np.zeros(nrec)
    for c in range(channels):
        f=fs+c*df
        t_dsp=K*DM*(1.0/(fm*fm)-1.0/(f*f))
        rec_shift=t_dsp/dt
        arr=copy.copy(dyn_spc[:,c])
        med=np.median(arr)
        prof=prof+shift(arr,rec_shift,order=1,mode='constant', cval=med)
    return prof/channels

if __name__=='__main__':

    floatsize=4
    timesize=8
    channels=4096
    recl=channels*floatsize
    bufsize=timesize+recl

    fname="svfits_postcorr.dat"
    chan_sel_file="svfits.log"

    freq_start=0.550 #GHz
    freq_end=0.750   #GHz
    integ=1.310720e-3 # sec
    DM=28
    mjd=60676.59087528
    
    # parse command line arguments
    parser=argparse.ArgumentParser()
    parser.add_argument("-d","--dm",help="dispersion measuure",type=float)
    parser.add_argument("-f","--freq", help="start and stop frequency (GHz)",
                        nargs=2, type=float)
    parser.add_argument("-i","--integ",help="integration time (sec)",type=float)
    parser.add_argument("-l","--log_file",help="svfits log file")
    parser.add_argument("-m","--mjd",help="modified julian date of burst",
                        type=float)
    parser.add_argument("-p","--post_corr",help="post_corr data file")

    args=parser.parse_args()
    if args.post_corr is not None:
        fname=args.post_corr
        if args.log_file is not None:
            chan_sel_file=args.log_file

    # read in the post-correlation beam
    hdr_size=7*8 # size of the postcorrelation data header
    try:
        fp=open(fname,"rb")
    except:
        print("Unable to open ",fname)
        exit(1)
    fp.seek(0,2)
    file_size=fp.tell()-hdr_size
    fp.seek(0,0)
    n_rec=file_size//recl
    #read in the header of the file
    rbuf=fp.read(hdr_size) # size of the header
    mjd,int_wd,DM,freq_start,freq_end,integ=np.frombuffer(rbuf,dtype=np.float64,
                                                          count=6)
    freq_start=freq_start/1.0e9 #GHz
    freq_end=freq_end/1.0e9  #GHz
    channels=np.frombuffer(rbuf[48:],dtype=np.int32,count=1)[0]

    #overwrite with user supplied values (if present)
    if args.dm is not None:
        DM=args.dm
    if args.freq is not None:
        freq_start,freq_end=args.freq
    if args.integ is not None:
        integ=args.integ
    if args.mjd is not None:
        mjd=args.mjd

    # read in the data from the file
    time_stamp=np.zeros(n_rec)
    dyn_spc=np.zeros((n_rec,channels))
    for r in range(n_rec):
        rbuf=fp.read(bufsize)
        time_stamp[r]=np.frombuffer(rbuf,dtype=np.float64,count=1)[0]
        dyn_spc[r]=np.frombuffer(rbuf[8:],dtype=np.float32)

    # read in the data selection done by svfits
    chan0=np.zeros(n_rec)
    chan1=np.zeros(n_rec)
    time_stamp1=np.zeros(n_rec)
    try:
        file=open(chan_sel_file,mode="r")
    except:
        print("Unable to open ",chan_sel_file)
        exit(1)
    r=0
    got_start=False
    while(True):
        line=file.readline()
        if not len(line):
            break # reached EoF
        word_list=line.split()
        if word_list[0] != "COPY:File":
            continue
        if r== n_rec:
            print("Warning:Trucated selection overplot");
            break;
        got_start=True
        time_stamp1[r]=float(word_list[7])
        chan0[r]=int(word_list[9])
        chan1[r]=int(word_list[11])
        r=r+1
    n_rec1=r        

    # ensure that the arrays are in time order
    ts_idx=time_stamp.argsort()
    ts=time_stamp[ts_idx]
    ds=dyn_spc[ts_idx,:]
    ts1_idx=time_stamp1[:n_rec1].argsort()
    ts1=time_stamp1[ts1_idx]
    c0=chan0[ts1_idx]
    c1=chan1[ts1_idx]
    
    # plot the post-correlation beam and data selection
    fig,ax=plt.subplots(1,1)
    med=np.median(ds)
    mad=stats.median_abs_deviation(ds,axis=None)
    v_min=med-3*mad
    v_max=med+3*mad
    y_min=ts[0]
    y_max=ts[n_rec-1]
    x_min=0
    x_max=channels
    im=ax.imshow(ds,origin='lower',vmin=v_min,vmax=v_max,
                 cmap=plt.get_cmap('viridis'),extent=(x_min,x_max,y_min,y_max))
    fig.colorbar(im)
    extent=im.get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2])))
    ax.scatter(c0,ts1,marker='+',s=5,color='blue')
    ax.scatter(c1,ts1,marker='+',s=5,color='red')
    title="Area Identified as having Burst Signal"
    ax.set_title(title)
    ax.set_xlabel("Channels")
    ax.set_ylabel("IST Time (sec)")
    plt.show()
    plt.close()

    # create and plot the de-deispersed profile
    prof=dedisperse(ds,integ,freq_start,freq_end,DM)
    med=np.median(prof)
    mad=np.median(np.abs(prof-med))
    snr=prof.max()/mad
    fig,ax=plt.subplots(1,1)
    ax.plot(ts,prof)
    jd=mjd+2400000.5
    burst_time=Time(jd,format='jd')
    ist_sec=(burst_time.ymdhms)['hour']*3600+(burst_time.ymdhms)['minute']*60
    ist_sec=ist_sec+(burst_time.ymdhms)['second']+5*3600.0+30*60.0
    if ist_sec > 86400.0:
        ist_sec=ist_sec-86400
    fm=0.5*(freq_start+freq_end)
    if(freq_start>freq_end):
        fh=freq_start
    else:
        fh=freq_end
    ist_sec=ist_sec + K*DM*(1.0/(fm*fm)-1.0/(fh*fh))
    plt.axvline(x=ist_sec,color='red')
    title="Incoherent de-dispersion 1.3ms Visibility\n" 
    title=title+str(burst_time.iso)+" ("+str(fh)+" GHz) [UTC]\n"
    title=title+f"Med={med:12.6e} Mad={mad:12.6e} SNR={snr:10.4e}"
    ax.set_title(title)
    ax.set_xlabel("IST Time ("+str(fm)+" GHz)")
    ax.set_ylabel("Power")
    plt.show()
    plt.close()
    


