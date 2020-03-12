#! /usr/bin/env python
import sys
from argparse import ArgumentParser
import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import casacore.tables as pt 
import matplotlib.dates as mdates
from datetime import datetime as dt


parser = ArgumentParser("create dynspectrum plot from raw data");
parser.add_argument('-i','--rawdata', nargs='+',help='list of rawdata files',dest="rawfile",required=True)
parser.add_argument('-k','--skip_samples', help='number of samples to skip',dest='skip_samples',type=int,default=0)
parser.add_argument('-s','--vmin', help='vmin of plot',dest='vmin',type=float,default=0.5)
parser.add_argument('-m','--vmax', help='vmax of plot',dest='vmax',type=float,default=2)
parser.add_argument('-c','--cmap', help='matplotlib colormap',dest='cmap',default="viridis")
parser.add_argument('-n','--show_normalization', help='plot normalization',dest='show_normalization', action='store_true')

#To DO: fix complex_voltages plotting (make sure they start at the same time). Save as png. Do not keep full buffer in memory. Get info from h5files. Maybe keep x-axis fixed. # plotmedian #add axes etc #find parset# figure ut how to stop the plotting

bytes_per_sample=4

def get_metadata_from_h5(h5file):
     #metadata=h5file.attrs[u'NOF_SUB_ARRAY_POINTINGS'] 
     metadata=dict(h5file[h5file.visit(lambda x: x if 'STOKES' in x else None)].attrs)
     metadata['freqs'] = h5file[h5file.visit(lambda x: x if 'COORDINATE_1' in x else None)].attrs[u'AXIS_VALUES_WORLD']
     metadata=dict(metadata,**dict(h5file[h5file.visit(lambda x: x if 'BEAM' in x else None)].attrs))
     metadata["starttime"]= h5file.attrs[u'OBSERVATION_START_MJD']
     metadata["endtime"]= h5file.attrs[u'OBSERVATION_END_MJD']
     return metadata


def plot_real_time2(fig,axarr,rawfile,nch,nSB,freqs,vmin,vmax,readTime=10,averagetime=1,averagefreq=.2,skipSamples=0,starttime=None,endtime=None,sampleSize=1./125.,plottime=1800,cmap='Reds',show_norm=False):
     rawfile.seek(skipSamples*bytes_per_sample*nch*nSB)
     '''averagetime in seconds,averagefreq in MHz'''
     
     dfreq = (freqs[-1]-freqs[0])/freqs.shape[0]
     averagech=int(averagefreq*1e6//dfreq)
     print ("freq",dfreq,averagefreq,averagech,freqs.shape)
     if not starttime is None:
          starttime += (skipSamples*sampleSize)/(24*3600.)
          starttime_dt =dt.strptime(pt.taql('select str(mjdtodate({})) as date'.format(starttime)).getcol('date')[0], '%Y/%m/%d/%H:%M:%S')
          fig.suptitle(starttime_dt.strftime("%m/%d/%Y"))
     data=np.zeros((plottime*averagetime,nSB*nch//averagech),dtype=np.float32)
     time=starttime
     currentsample=0
     readSamples=int(readTime//sampleSize)
     averageSamples=int(averagetime//sampleSize)
     averagetime=averageSamples*sampleSize
     readSamples-=readSamples%averageSamples  #make sure you read an integer number of averageSamples
     ax=axarr
     if show_norm:
          ax=ax[0]
     ax.cla()
     myextent=[0,data.shape[0],freqs[0]*1e-6,freqs[-1]*1e-6]
     myextent[0]=mdates.date2num(starttime_dt)
     myextent[1]=mdates.date2num(dt.strptime(pt.taql('select str(mjdtodate({})) as date'.format(starttime+plottime/(24.*3600.))).getcol('date')[0], '%Y/%m/%d/%H:%M:%S'))  # thanks to TJD
     myimshow = ax.imshow((data).T,origin='lower',interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax,extent=myextent,cmap=cmap)
     ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
     ax.set_ylabel("freq (MHz)")
     if show_norm:
          ax2=axarr[1]
          ax2.cla()
          mymedian = np.ones(nSB*nch//averagech,dtype=np.float)
          ax2.plot(freqs[::averagech]*1e-6,mymedian,'k')
          ax2.axes.set_yscale("log")
          ax2.set_xlabel("freq (MHz)")
     while(time<endtime):
          mybuffer = rawfile.read(readSamples*bytes_per_sample*nch*nSB) # TODO:before reading check that data is there?
          tmpdata = np.frombuffer(mybuffer,dtype=np.float32) #N.B.: np.float = np.float64!!
          if not tmpdata.shape:
               #out of data increase waittime
               waittime=readTime-1  #1 second for processing?
               sys.sleep(10)
               continue
          else:
               waittime=0.01
          nSam = tmpdata.shape[0]//(nch*nSB)
          tmpdata = np.average(tmpdata.reshape(nSam,(nch*nSB))[:nSam-nSam%averageSamples].reshape((-1,averageSamples,nch*nSB)),axis=1)
          if averagech>1:
               tmpdata = np.average(tmpdata[:,:(nch*nSB)-(nch*nSB)%averagech].reshape((tmpdata.shape[0],-1,averagech)),axis=2)
          if currentsample+tmpdata.shape[0]>data.shape[0]:
               if currentsample<data.shape[0]:
                    data[:-tmpdata.shape[0]] = data[currentsample-(data.shape[0]-tmpdata.shape[0]):currentsample]
               else:
                    data[:-tmpdata.shape[0]]=data[tmpdata.shape[0]:]
               data[-tmpdata.shape[0]:]=tmpdata[:]
               mymedian=np.median(data,axis=0)
               plottimestart=dt.strptime(pt.taql('select str(mjdtodate({})) as date'.format(time+(tmpdata.shape[0]*averagetime)/(24.*3600.) - plottime/(24.*3600.))).getcol('date')[0], '%Y/%m/%d/%H:%M:%S')
               myextent[0]=mdates.date2num(plottimestart)
               plottimeend=dt.strptime(pt.taql('select str(mjdtodate({})) as date'.format(time+(tmpdata.shape[0]*averagetime)/(24.*3600.))).getcol('date')[0], '%Y/%m/%d/%H:%M:%S')
               myextent[1]=mdates.date2num(plottimeend)
          else:
               data[currentsample:currentsample+tmpdata.shape[0]]=tmpdata[:]
               mymedian=np.median(data[:currentsample],axis=0)
          currentsample+=tmpdata.shape[0]

          time+=(tmpdata.shape[0]*averagetime)/(24.*3600.)
          myimshow.set_data((data/mymedian).T)
          myimshow.set_extent(myextent)
          if show_norm:
               ax2.cla()
               ax2.plot(freqs[::averagech]*1e-6,mymedian,'k')
               ax2.axes.set_yscale("log")
               ax2.set_xlabel("freq (MHz)")
          plt.pause(waittime)
#imshow.setdata(data/mymedian)
         

def plot_real_time(fig,axarr,rawfile,nch,nSB,freqs,vmin,vmax,maxSamples=1000,skiptime=125,skipch=1,skipSamples=0,starttime=None,endtime=None,sampleSize=1./125.,cmap='Reds',show_norm=False):
     hasdata=False
     rawfile.seek(skipSamples*bytes_per_sample*nch*nSB)
     if not starttime is None:
          starttime += (skipSamples*sampleSize)/(24*3600.)
          starttime_dt =dt.strptime(pt.taql('select str(mjdtodate({})) as date'.format(starttime)).getcol('date')[0], '%Y/%m/%d/%H:%M:%S')
          fig.suptitle(starttime_dt.strftime("%m/%d/%Y"))
     while(True):   
        mybuffer = rawfile.read(maxSamples*bytes_per_sample*nch*nSB)
        tmpdata = np.frombuffer(mybuffer,dtype=np.float32) #np.float = np.float64!!
        nSam = tmpdata.shape[0]//(nch*nSB)
        tmpdata = np.average(tmpdata.reshape(nSam,(nch*nSB))[:nSam-nSam%skiptime].reshape((-1,skiptime,nch*nSB)),axis=1)
        tmpdata = np.average(tmpdata[:,:(nch*nSB)-(nch*nSB)%skipch].reshape((tmpdata.shape[0],-1,skipch)),axis=2)
        if not hasdata:
             #data=tmpdata.reshape(nSam,(nch*nSB))[::skiptime,::skipch]
             data=tmpdata[:]
             hasdata=True
        else:
             data=np.concatenate((data,tmpdata),axis=0)  #faster fill predefined array
        mymedian=np.median(data,axis=0)
        #fig.clf()
        ax=axarr
        if show_norm:
             ax=ax[0]
        ax.cla()
        myextent=[0,data.shape[0],freqs[0]*1e-6,freqs[::skipch][-1]*1e-6]
        if not (starttime is None):
             myextent[0]=mdates.date2num(starttime_dt)
             myextent[1]=mdates.date2num(dt.strptime(pt.taql('select str(mjdtodate({})) as date'.format(starttime+(data.shape[0]*skiptime*sampleSize)/(24.*3600.))).getcol('date')[0], '%Y/%m/%d/%H:%M:%S'))  # thanks to TJD
             ax.imshow((data/mymedian).T,origin='lower',interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax,extent=myextent,cmap=cmap)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_ylabel("freq (MHz)")
        if show_norm:
             ax=axarr[1]
             ax.cla()
             ax.plot(freqs[::skipch]*1e-6,mymedian,'k')
             ax.set_xlabel("freq (MHz)")
        #print ("waiting",0.1*maxSamples*sampleSize,"seconds")
        plt.pause(0.1*maxSamples*sampleSize)
        #savefig("fig.png")


def plot_real_time_complex_voltages(rawfiles,nch,nSB,freqs,vmin,vmax,maxSamples=1000,skiptime=100,skipch=1):
     hasdata=False
     buffers=['','','','']
     prevminsize=0
     while True:
        minsize=-1
        for i in xrange(len(rawfiles)):
            buffers[i] = buffers[i][prevminsize:]+rawfiles[i].read(maxSamples*bytes_per_sample*nch*nSB)
            #print minsize,len(buffers[i])
            if minsize<0 or len(buffers[i])<minsize:
                minsize=len(buffers[i])
        prevminsize=minsize
        tmpdata=np.frombuffer(buffers[0][:minsize],dtype=np.float32)
        nSam=tmpdata.shape[0]/(nch*nSB)
        print "reshaping",tmpdata.shape[0],"to",nSam,nch*nSB,minsize
        tmpdata=tmpdata.reshape(nSam,(nch*nSB))[::skiptime,::skipch]
        tmpdata= tmpdata**2 + \
              np.frombuffer(buffers[1][:minsize],dtype=np.float32).reshape(nSam,(nch*nSB))[::skiptime,::skipch]**2 + \
              np.frombuffer(buffers[2][:minsize],dtype=np.float32).reshape(nSam,(nch*nSB))[::skiptime,::skipch]**2 + \
              np.frombuffer(buffers[3][:minsize],dtype=np.float32).reshape(nSam,(nch*nSB))[::skiptime,::skipch]**2 
        tmpdata=np.sqrt(tmpdata)
        if not hasdata:
             data=tmpdata
             hasdata=True
        else:
             data=np.concatenate((data,tmpdata),axis=0)
        mymedian=np.median(data,axis=0)
        clf()
        subplot(2,1,1)
        imshow((data/mymedian).T,origin='lower',interpolation='nearest',aspect='auto',vmin=vmin,vmax=vmax,extent=(0,data.shape[0],freqs[0]*1e-6,freqs[::skipch][-1]*1e-6))
        ylabel("freq (MHz)")
        subplot(2,1,2)
        plot(freqs[::skipch]*1e-6,mymedian,'k')
        xlabel("freq (MHz)")
        pause(.1)
        #savefig("fig.png")

def main(argv):
    args=parser.parse_args(argv)
    skipSamples=args.skip_samples
    rawfiles=[open(i,'rb') for i in args.rawfile]
    metadata = get_metadata_from_h5(h5py.File(args.rawfile[0].replace('.raw','.h5')))
    #a=[rawfile.seek(-1000*60*16*4,os.SEEK_END) for rawfile in rawfiles]
    #buffers=[rawfile.read(1000*400*1*4) for rawfile in rawfiles]
    if not metadata[u'COMPLEX_VOLTAGE']:
            for i,rawfile in enumerate(rawfiles):
                 fig,axarr=plt.subplots(1+args.show_normalization,1)
                 plot_real_time2(fig,axarr,rawfile,metadata['CHANNELS_PER_SUBBAND'],metadata[u'NOF_SUBBANDS'],metadata['freqs'],args.vmin,args.vmax,skipSamples=skipSamples,starttime=metadata['starttime'],endtime=metadata['endtime'],sampleSize=metadata[u'SAMPLING_TIME'],cmap=args.cmap,show_norm=args.show_normalization)
    else:
            for i in range(len(rawfiles))[::4]:
                rawfile4=rawfiles[i*4:i*4+4]
                fig,axarr=plt.subplots(2,1)
                plot_real_time_complex_voltages(fig,axarr,rawfile4,metadata['CHANNELS_PER_SUBBAND'],metadata[u'NOF_SUBBANDS'],metadata['freqs'],args.vmin,args.vmax)
            

if __name__ == '__main__':
    main(sys.argv[1:])    
