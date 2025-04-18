import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticky
from   scipy import stats
import argparse

floatsize=4
timesize=8
channels=4096
recl=channels*floatsize
bufsize=timesize+recl

fname="svfits_postcorr.dat"
chan_sel_file="svfits.log"

parser=argparse.ArgumentParser()
parser.add_argument("-p","--post_corr",help="post_corr data file")
parser.add_argument("-l","--log_file",help="svfits log file")
args=parser.parse_args()
if args.post_corr is not None:
    fname=args.post_corr
if args.log_file is not None:
    chan_sel_file=args.log_file
    
try:
    fp=open(fname,"rb")
except:
    print("Unable to open ",fname)
    exit(1)
fp.seek(0,2)
file_size=fp.tell()
fp.seek(0,0)
n_rec=file_size//recl
time_stamp=np.zeros(n_rec)
dyn_spc=np.zeros((n_rec,channels))
for r in range(n_rec):
    rbuf=fp.read(bufsize)
    time_stamp[r]=np.frombuffer(rbuf,dtype=np.float64,count=1)[0]
    dyn_spc[r]=np.frombuffer(rbuf[8:],dtype=np.float32)

chan0=np.zeros(n_rec)
chan1=np.zeros(n_rec)
time_stamp1=np.zeros(n_rec)
try:
    file=open(chan_sel_file,mode="r")
except:
    print("Unable to open ",fname)
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
fig,ax=plt.subplots(1,1)
med=np.median(dyn_spc)
mad=stats.median_abs_deviation(dyn_spc,axis=None)
v_min=med-3*mad
v_max=med+3*mad
y_min=time_stamp[0]
y_max=time_stamp[n_rec-1]
x_min=0
x_max=channels
im=ax.imshow(dyn_spc,origin='lower',vmin=v_min,vmax=v_max,cmap=plt.get_cmap('viridis'),extent=(x_min,x_max,y_min,y_max))
fig.colorbar(im)
extent=im.get_extent()
ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2])))
ax.scatter(chan0[:n_rec1],time_stamp1[:n_rec1],marker='+',s=5,color='red')
ax.scatter(chan1[:n_rec1],time_stamp1[:n_rec1],marker='+',s=5,color='red')
ax.set_title("Area Identified as having Burst Signal")
ax.set_xlabel("Channels")
ax.set_ylabel("Records")
plt.show()
plt.close()


    


