# Function to read a SPOTLIGHT 1.3ms raw visibility file, extract (and average
# if needed) selected data and make a plots. This program is mainly meant for
# short term, debugging use, and probably does not cater for all possibilities
# that may arise in the SPOTLIGHT program.
#
# Required arguments are the list of filenames with the raw data (--files) as
# well as the list of indices of these files (--index). The timestamp read from
# the file is addjusted for the file order using the supplied index of the file.
#
# The main modes are :
#   print the polarization averaged data for some antenna pair (--ants)
#   print the polarization averaged data for all selfs (--self)
#   print the polarization averaged data for all cross (default)
#
# In addition, the user has a choice to do a coherent or incoherent sum of
# the visiblities [--coherent], collate multiple slices [--slice] from multiple
# visibility files [--collate] (useful for example when a particular burst
# covers multiple slices in different files), apply an amplitude bandpass
# [--bandpass], subtract the median visibility over the slice [--baseline]
# drop the central square baselines [--drop_csq] while doing the sum, as
# well as compute the statistics of the number of NaNs/Infinities per
# channel and baseline, overplot the data selection done by svfits, to check
# if the burst signal has been correctly identified[--overplot_sel]
#
# jnc apr 2025
#
import sys
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticky
from   datetime import datetime,timedelta
import argparse
import copy

MaxAnt=32
MaxBase=(MaxAnt*(MaxAnt+1))//2
Channels=4096
Stokes=2
FloatSize=4
TimeSize=16
Recl=MaxBase*Channels*Stokes*FloatSize
NFiles=16
RecPerSlice=50
Integ=1.30172e-3 # lta integration time (sec)
SliceDuration=RecPerSlice*Integ
SliceInterval=NFiles*SliceDuration
AntMask=1073741823 #30 antennas
SliceDate=""
ShortBase=[34,37,38,40,67,68,70,95,152,154,180]#(C03-C04),(C01,C02,C05,C06,C09)
BaseFlags=np.zeros(MaxBase*Stokes)
ChanFlags=np.zeros(Channels)

def get_data(fp,idx,rec,self,coherent,sel_ant0,sel_ant1,drop_csq,nan_stats):
    """Returns [average] of selected data in a 1.3ms SPOTLIGHT raw visibility file."""

    amp=np.zeros(Channels,dtype=np.float32)
    re=np.zeros(Channels,dtype=np.float32)
    im=np.zeros(Channels,dtype=np.float32)
    global SliceDate,ChanFlags,BaseFlags
    
    # read the timestamp (if appropriate) and the visibility data
    time_off=0
    bufsize=Recl
    n_slice=rec//RecPerSlice
    off=rec*Recl+n_slice*TimeSize
    fp.seek(off,0)
    if rec % RecPerSlice ==0:
        time_off=TimeSize
        bufsize=Recl+TimeSize
    try:
        rbuf1=fp.read(bufsize)
    except:
        print("EoF reached")
        return [0],[0]
    if time_off > 0: # we have a timestamp
        tv=np.frombuffer(rbuf1,dtype=np.int64,count=2) #timeval is signed long
        rec_time=tv[0]+tv[1]/1.0e6+idx*SliceDuration #offset for this file
        SliceDate=str(datetime.fromtimestamp(rec_time))
    rbuf=rbuf1[time_off:]
    vis=np.frombuffer(rbuf,dtype=np.float16) #re,im are half float
    if len(vis) < 2*MaxBase*Stokes*Channels:
        print("EoF reached")
        return [0],[0]

    # read in the selfs and average over all antennas and all pols
    for ant0 in range(0,MaxAnt):
        if not (1<<ant0) & AntMask:
            continue
        if sel_ant0 >=0 and ant0 != sel_ant0:
            continue
        for ant1 in range(ant0,MaxAnt):
            if not (1<<ant1) & AntMask:
                continue
            if sel_ant1 >= 0 and ant1 !=sel_ant1:
                continue
            if self and (ant0 != ant1):
                continue;
            if not self and (ant0 == ant1):
                continue
            ibase=0
            found=False
            for a0 in range(0,MaxAnt):
                for a1 in range(a0,MaxAnt):
                    if a0==ant0 and a1==ant1:
                        found=True; break
                    ibase=ibase+1
                if found:
                    break
            if not found:
                print("Could not find baseline ",ant0,ant1)
                exit(1)
            #drop short baselines unless we have actually asked for them
            if sel_ant0 <0 and sel_ant1 <0 and not self:
                if ibase in ShortBase:
                    continue
                if drop_csq and (ant0 < 11 or ant1 < 11):
                    continue
            off=ibase*Channels*2
            for pol in [0,1]:
                if pol ==1:
                    off=off+2*MaxBase*Channels
                vis1=vis[off:]
                if nan_stats: # per baseline & per channel NaN/Infinity
                    count=np.count_nonzero(np.isnan(vis1[:2*Channels]))
                    count=count+np.count_nonzero(np.isinf(vis1[:2*Channels]))
                    BaseFlags[ibase+pol*MaxBase]=BaseFlags[ibase+pol*MaxBase]+count
                    ChanFlags=ChanFlags+np.isnan(vis1[0:2*Channels:2])
                    ChanFlags=ChanFlags+np.isnan(vis1[1:2*Channels:2])
                    ChanFlags=ChanFlags+np.isinf(vis1[0:2*Channels:2])
                    ChanFlags=ChanFlags+np.isinf(vis1[1:2*Channels:2])
                re1=np.nan_to_num(vis1[0:2*Channels:2].astype(float),
                                  nan=0.0, posinf=0.0,neginf=0.0)
                im1=np.nan_to_num(vis1[1:2*Channels:2].astype(float),
                                  nan=0.0,posinf=0.0,neginf=0.0)
                if coherent:
                    re=re+re1
                    im=im+im1
                else:
                    amp=amp+np.sqrt(re1*re1+im1*im1)
    if coherent:
        return re,im
    else:
        return amp,0*amp

if __name__=='__main__':

    ant0=-1
    ant1=-1
    cfname=""
    chan_sel_file=""
    parser=argparse.ArgumentParser()
    parser.add_argument("-a","--ants", help="antenna pair in baselines",
                        nargs=2, type=int)
    parser.add_argument("-b","--bandpass", action="store_true",
                        help="do amp bandpass calibration")
    parser.add_argument("-B","--baseline", action="store_true",
                        help="subtract median from visibilities")
    parser.add_argument("-C","--collate", action="store_true",
                        help="collate multiple slices")
    parser.add_argument("-d","--drop_csq", action="store_true",
                        help="drop central square baselines")
    parser.add_argument("-f","--files", help="list of files to process",
                        nargs='+',required=True)
    parser.add_argument("-i","--index", help="index of the files to process",
                        nargs='+',type=int,required=True)
    parser.add_argument("-I","--interactive", action="store_true",
                        help="interactive plotting")
    parser.add_argument("-n","--nan_stats", action="store_true",
                        help="compute and plot NaN statistics")
    parser.add_argument("-o","--overplot", nargs=1,
                        help="overplot channels selection in svfits")
    parser.add_argument("-p","--power", action="store_true",
                        help="add visibility powers [inchoherent addition]")
    parser.add_argument("-s","--slice", help="[start and [stop] slice num",
                        nargs='*',type=int)
    parser.add_argument("-S","--self", action="store_true",
                        help="plot the sum all self visibilities")

    
    args=parser.parse_args()
    bandpass=args.bandpass
    baseline=args.baseline
    coherent=not args.power
    collate=args.collate
    drop_csq=args.drop_csq
    interactive=args.interactive
    nan_stats=args.nan_stats
    self=args.self
    files=args.files
    index=args.index
    if args.ants is not None:
        ant0,ant1=args.ants
        drop_csq=False;self=False;
    if len(files) != len(index):
        print("Number of files does not match number of indices")
        parser.print_help()
        exit(1)
    if args.overplot != None:
        overplot_sel=True
        chan_sel_file=args.overplot
    else:
        overplot_sel=False
    if args.slice != None:
        slice_a=np.array(args.slice,dtype=np.int32)
    else:
        slice_a = np.array([0,-1],dtype=np.int32)
    if collate:
        if args.slice==None:
            print("Slice numbers must be specified in collate mode")
            parser.print_help()
            exit(1)
        if len(slice_a) != len(files):
            print("For collate mode Nfiles must be == Nslices")
            parser.print_help()
            exit(1)

    recnum=np.zeros(len(files)*RecPerSlice)
    chan0=np.zeros(len(files)*RecPerSlice)
    chan1=np.zeros(len(files)*RecPerSlice)
    if overplot_sel:
        with open(chan_sel_file,"r") as file:
            for i,line in enumerate(file):
                if i== len(files)*RecPerSlice:
                    print("Warning:Trucated selection overplot");
                    break;
                r,c0,c1=line.split()
                recnum[i]=int(r)
                chan0[i]=int(c0)
                chan1[i]=int(c1)
            nrec=i
         
    ds_re=np.zeros([RecPerSlice,Channels],dtype=np.float32)
    ds_im=np.zeros([RecPerSlice,Channels],dtype=np.float32)
    dyn_spc1=np.zeros([RecPerSlice,Channels],dtype=np.float32)
    if collate:
        dyn_spc2=np.zeros([RecPerSlice*len(files),Channels],dtype=np.float32)
    for ifile,ftuple in enumerate(zip(files,index)):
        fname,idx=ftuple
        try:
            fp=open(fname,"rb")
        except:
            print("Unable to open ",fname)
        if collate:
            r=int(slice_a[ifile])*RecPerSlice
        else:
            r=slice_a[0]*RecPerSlice
        if collate:
            slice_date0=copy.deepcopy(SliceDate) #start first slice collation
        while(True):
            re,im=get_data(fp,idx,r,self,coherent,ant0,ant1,drop_csq,nan_stats)
            if len(re) < Channels:
                print("Short Read: EoF reached for %s after %d recs" %(fname,r))
                break
            else:
                ds_re[r%RecPerSlice]=re; ds_im[r%RecPerSlice]=im
            r=r+1
            slice=r//RecPerSlice-1
            if (r >0 and r % RecPerSlice) == 0:
                dyn_spc1=np.sqrt(np.square(ds_re)+np.square(ds_im))
                median_a=np.median(dyn_spc1,axis=0)
                med=np.median(median_a)
                if baseline:
                    ds_re=ds_re-np.median(ds_re,axis=0)
                    ds_im=ds_im-np.median(ds_im,axis=0)
                if bandpass:
                    dyn_spc1=np.sqrt(np.square(ds_re)+np.square(ds_im))
                    dyn_spc1=med*(dyn_spc1/median_a[np.newaxis,:])
                if collate:
                    r1=ifile*RecPerSlice
                    r2=(ifile+1)*RecPerSlice
                    dyn_spc2[r1:r2]=dyn_spc1
                    if ifile==len(files)-1:
                        dyn_spc=dyn_spc2
                        slice_date=slice_date0
                        fname=cfname+fname+"Slice"+str(slice_a[ifile])
                    else:
                        cfname=cfname+fname+"Slice"+str(slice_a[ifile])
                        break; # go to next file
                else:
                    fname=fname+"Slice"+str(slice)
                    slice_date=SliceDate
                    dyn_spc=dyn_spc1
                mad_a=stats.median_abs_deviation(dyn_spc,nan_policy='omit')
                mad=np.median(mad_a)
                med=np.median(dyn_spc)
                min=dyn_spc.min()
                max=dyn_spc.max()
                print(f"Min= %12.4e Max= %12.4e Med= %12.4e Mad= %12.4e" %
                      (min,max,med,mad))
                fig,ax=plt.subplots(1,1)
                if not bandpass:
                    im=ax.imshow(dyn_spc,origin='lower',
                                 cmap=plt.get_cmap('viridis'),
                                 norm=matplotlib.colors.LogNorm())
                else:
                    v_min=med-3*mad
                    v_max=med+3*mad
                    im=ax.imshow(dyn_spc,origin='lower',vmin=v_min,vmax=v_max,
                                 cmap=plt.get_cmap('viridis'))
                fig.colorbar(im)
                extent=im.get_extent()
                ax.set_aspect(abs((extent[1]-extent[0])/
                                  (extent[3]-extent[2])))
                ax.set_xlabel("Channels")
                ax.set_ylabel("Records")
                if overplot_sel:
                    ax.scatter(chan0[:nrec],recnum[:nrec],marker='+',
                               color='red')
                    ax.scatter(chan1[:nrec],recnum[:nrec],marker='+',
                               color='red')
                title=f"Dynamic Spectrum {fname}\n Slice {slice} {slice_date}\n Min={min:10.2e} Max={max:10.2e} Med={med:10.2e} Mad={mad:10.2e}"
                ax.set_title(title)
                if self:
                    foot="Total Power for all Antennas"
                else:
                    if ant0 >=0 and ant1 >= 0:
                        foot="Visibility Amp ["+str(ant0)+"-"+str(ant1)+"]"
                    else:
                        if coherent:
                            foot="Coherent Sum of All Visibilities"
                        else:
                            foot="Inoherent Sum of All Visibilities"
                plt.annotate(foot,xy=(0.3,0.0),xycoords='figure fraction')
                if interactive:
                    plt.show()
                plt.savefig(fname+".pdf")
                plt.close()
                if collate:
                    exit(1) # all done
                else:
                    if slice_a[1] >= 0 and slice >= slice_a[1]:
                        break # go to next file
                if nan_stats:
                    x=np.arange(0,2*MaxBase)
                    fig,ax=plt.subplots()
                    ax.scatter(x,BaseFlags,color='blue')
                    ax.set_xlabel("Baselines",color='blue')
                    ax.tick_params(axis='x',labelcolor='blue')
                    
                    ax1=ax.twiny()
                    x1=np.arange(0,Channels)
                    ax1.scatter(x1,ChanFlags,color='red')
                    ax1.set_xlabel("Channels",color='red')
                    ax1.tick_params(axis='x',labelcolor='red')
                    plt.show()
                    plt.close()
                            
