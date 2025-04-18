import sys
import numpy as np
from   scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticky
import argparse


flag_file="svfits.log"
log_scale=1
parser=argparse.ArgumentParser()
parser.add_argument("-l","--log_file",help="svfits log file")
parser.add_argument("-L","--linear_scale",action="store_true",
                    help="use linear scale (default log")
args=parser.parse_args()
if args.log_file is not None:
    flag_file=args.log_file
if args.linear_scale is not None:
    log_scale=not args.linear_scale
    

try:
    file=open(flag_file,mode="r")
except:
    print("Unable to open ",flag_file)
    exit(1)

flag_stats=np.zeros((32,32))    
while(True):
    got_start=False
    file_end=False
    while not file_end:
        line=file.readline()
        if len(line)==0:
            print("EoF for "+flag_file)
            exit(0)
        word_list=line.split()
        if word_list[0]!= "FLAG:File":
            if got_start:
                break
            continue
        got_start=True
        file_idx=word_list[1]
        slice=word_list[3]
        a0=int(word_list[5])
        a1=int(word_list[7])
        flag_stats[a0][a1]=int(word_list[9])
        flag_stats[a1][a0]=int(word_list[11])
                

    fig,ax=plt.subplots(1,1)
    if not log_scale:
        med=np.median(flag_stats)
        mad=stats.median_abs_deviation(flag_stats,axis=None)
        v_min=med-30*mad
        v_max=med+30*mad
        im=ax.imshow(flag_stats,origin='lower',vmin=v_min,vmax=v_max,
                     cmap=plt.get_cmap('viridis'))
    else:
        im=ax.imshow(flag_stats,origin='lower',cmap=plt.get_cmap('viridis'),
                     norm=matplotlib.colors.LogNorm())
        title="Flags for File "+file_idx+" Slice "+slice
        title=title+"\nTopLeft(RR) BotRight(LL)"
        ax.set_title(title)
        ax.set_xlabel("Antenna Number")
        ax.set_ylabel("Antenna Number")
        fig.colorbar(im)
        plt.show()
        plt.close()
        
