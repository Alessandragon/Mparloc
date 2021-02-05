#!/usr/bin/python3
from obspy import read
from functools import partial
import threading
import numpy as np
import os
from nllgrid import NLLGrid
from obspy.core.trace import Trace
import time
import glob
import math
import pyproj
from obspy.taup import *
from obspy.taup.taup_geo import calc_dist
from obspy.taup import taup_create
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel
from matplotlib.ticker import FormatStrFormatter
from threading import Lock
global lock_log
lock_log = Lock()
global map
import copy 
from array import array



def grid_creator(Network_folder,Event_folder):
 os.chdir(Network_folder+'\grids')
 grd=glob.glob('*buf')
 grd = NLLGrid(grd[0])
 dX=grd.dx#km
 dY=grd.dy #km
 dZ=grd.dz #km
 NodeX=grd.nx#km
 NodeY=grd.ny #km
 NodeZ=grd.nz #km
 Xo=grd.x_orig #km
 Yo=grd.y_orig #km
 Zo=grd.z_orig  #km
 Lat0=grd.orig_lat
 Lon0= grd.orig_lon
 lambert_area = {'proj': 'laea',
 'lat_0':Lat0, 
 'lon_0':Lon0, 
 'x_0':Xo, 
 'y_0':Yo,
 'ellps': 'WGS84',
 'datum': 'WGS84',
 'R':6378137.0}
 map = pyproj.Proj(lambert_area)
 os.chdir('..')
 os.chdir(Event_folder)
 file=sorted(glob.glob('MParLoc_event*'))	  
 if len(file)>0:
   for count in file:
     os.remove(count)
 Lista_sac=sorted(glob.glob('*'))
 Nomi=open('../list_stations_4loc.dat')
 Config=open('../../config_files/run_loc.conf')
 Nomi=Nomi.read()
 Nomi=Nomi.rsplit('\n')
 Config=Config.read()
 grid_precision=str(Config.rsplit('grid_precision: ')[1].rsplit('\n')[0])
 diff_Amp_b=float(Config.rsplit('diff_Amp_b: ')[1].rsplit('\n')[0])
 NAME=[]
 tr=[]
 for i in range(0,len(Lista_sac)):
           if 'Z' in read(Lista_sac[i])[0].stats.sac.kcmpnm or 'U' in read(Lista_sac[i])[0].stats.sac.kcmpnm:
            if read(Lista_sac[i])[0].stats.station in Nomi:
             NAME.append(read(Lista_sac[i])[0].stats.station)
             tr.append(read(Lista_sac[i])[0])	
 xx=[]
 yy=[]
 staz_lon=[]
 staz_lat=[]
 for i in range(0,len(tr)):
  lon = tr[i].stats.sac.stlo
  lat = tr[i].stats.sac.stla
  staz_lon.append(lon)
  staz_lat.append(lat)			 
 xx,yy=map(staz_lon,staz_lat)
 os.chdir('..')
 Baz=[]
 for i in range(0,len(tr)):
  Baz.append([0,0,0])
 ynode, xnode, znode = np.meshgrid(np.linspace(Yo*1000,(Yo+NodeY*dY-dY)*1000,NodeY),np.linspace(Xo*1000,(Xo+NodeX*dX-dX)*1000,NodeX), -np.linspace(Zo*1000,(Zo+NodeZ*dZ-dZ)*1000,NodeZ))
 BAznode=[]
 Rnode=[]
 Rdeltanode=[]
 RdeltanodeErr=[]
 TarrivalnodeErr=[]
 BAznodeErr=[]
 delta_arrival_matrix=[]
 os.chdir('grids')

 for z in range(0,len(NAME)):
  arrival_matrix=None
  arrival_matrix_S=None
  arrival_matrixId=None
  BAznodeId=None
  Rnode=None
 
 
  grid_list_arrivals=sorted(glob.glob('*.P.*'+NAME[z]+'*'+'.hdr'))		 
  if 'f' in   grid_precision:
   arrival_matrix=np.float32(NLLGrid(grid_list_arrivals[0]).array.reshape(((NodeX*NodeY*NodeZ))))
  if 'd' in   grid_precision:
   arrival_matrix=np.float64(NLLGrid(grid_list_arrivals[0]).array.reshape(((NodeX*NodeY*NodeZ))))
  arrival_matrixId=open('arrival_matrix'+NAME[z]+'.bin','w')
  arrival_matrix.tofile(arrival_matrixId)
 
  grid_list_arrivals=sorted(glob.glob('*.S.*'+NAME[z]+'*'+'.hdr'))
  if 'f' in   grid_precision:
   arrival_matrix_S=np.float32(NLLGrid(grid_list_arrivals[0]).array.reshape(((NodeX*NodeY*NodeZ))))
  if 'd' in   grid_precision:
   arrival_matrix_S=np.float64(NLLGrid(grid_list_arrivals[0]).array.reshape(((NodeX*NodeY*NodeZ))))
  arrival_matrixId=open('arrival_matrixS'+NAME[z]+'.bin','w') 
  arrival_matrix_S.tofile(arrival_matrixId)
  
  BAznodeId=open('BAznode'+NAME[z]+'.bin','w')  
  RnodeId=open('Rnode'+NAME[z]+'.bin','w')
  if 'f' in   grid_precision: 
   BAznode=np.float32(np.degrees(np.arctan2(xx[z]-xnode,yy[z]-ynode)+math.pi))
   Rnode=np.float32(diff_Amp_b*np.log10(np.sqrt(((xnode - xx[z])**2 + (yy[z] - ynode)**2)+znode**2)/1000))
  if 'd' in   grid_precision:
   BAznode=np.float64(np.degrees(np.arctan2(xx[z]-xnode,yy[z]-ynode)+math.pi))
   Rnode=np.float64(diff_Amp_b*np.log10(np.sqrt(((xnode - xx[z])**2 + (yy[z] - ynode)**2)+znode**2)/1000))
  BAznode.reshape(((NodeX*NodeY*NodeZ))).tofile(BAznodeId)
  BAznodeId.close()
  Rnode.reshape(((NodeX*NodeY*NodeZ))).tofile(RnodeId)
  RnodeId.close()
 os.chdir('../..')








class Location_MAP():
    def __init__(self, parent=None):
        self.figure = plt.figure()
        self.plot()

    def plot(self):
        self.figure.clear()		
        ab = self.figure.add_subplot(1, 1, 1)
        if len(xev)>0:
         ab.plot(xev,yev,'*',markersize=12,markerfacecolor=(1, 1,1, 1),markeredgecolor='k')
        ab.set_title('Location')
    def plot_location(self):
        self.figure.clear()
        spec=self.figure.add_gridspec(ncols=3, nrows=3, width_ratios=[1,4,1],height_ratios=[1,4,1])
		
		
		
###########################################################################################################################################	
        ab = self.figure.add_subplot(spec[2-1,2-1])
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        index=np.where(np.array(pick)!=0)[0]
        for i in range(0,len(xx)):
          ab.plot([xx[i]],[yy[i]],color='w',marker='^',markersize=6,markeredgecolor='k')
        if len(index)>0:
         for i in index:
          ab.plot([xx[i]],[yy[i]],'k^',markersize=6)
        cset=ab.contourf(win.xnode[:,:,index_LocZ[0]], win.ynode[:,:,index_LocZ[0]], win.Locmatrix[:,:,index_LocZ[0]]/np.max(win.Locmatrix),10,cmap='hot_r',vmin = 0, vmax = 1, alpha=.5, antialiased=True)
        if len(xev)>0:
         ab.plot(xev,yev,'*',markersize=12,markerfacecolor=(1, 1,1, 1),markeredgecolor='k')
        ab.plot(Locx,Locy,'*',markersize=12,markerfacecolor=(0, 0,0, 1),markeredgecolor='k')
        ab.set_xticks([])
        ab.set_yticks([])		
        cdfX1=[0]
        cdfX2=[0]
        cdfY1=[0]
        cdfY2=[0]
        cdfZ1=[0]
        cdfZ2=[0]
        pdfX=[]
        pdfY=[]
        pdfZ=[]
        dX=win.xnode[0,1,0]-win.xnode[0,0,0]
        dy=win.ynode[1,0,0]-win.ynode[0,0,0]
        dz=win.znode[0,0,1]-win.znode[0,0,0]
        for i in range(0,NodeX):
         pdfX.append((win.Locmatrix[i,index_LocY[0],index_LocZ[0]]))
        for i in range(0,NodeY):
         pdfY.append(((win.Locmatrix[index_LocX[0],i,index_LocZ[0]])))
        for i in range(0,NodeZ):
         pdfZ.append((win.Locmatrix[index_LocX[0],index_LocY[0],i]))
        Sum=pdfX[index_LocX[0]]
        Index_errorX=[0,len(pdfX)-1]	 
        for s in range(0,len(pdfX)):
         if index_LocX[0]-s-1 >= 0 and pdfX[index_LocX[0]-s] >= (Sum*win.th_pro):
          Index_errorX[0]=index_LocX[0]-s-1
         if index_LocX[0]+s+1 < len(pdfX) and pdfX[index_LocX[0]+s] >= (Sum*win.th_pro):
          Index_errorX[1]=index_LocX[0]+s+1
        Sum=pdfY[index_LocY[0]]
        Index_errorY=[0,len(pdfY)-1]		 
        for s in range(0,len(pdfY)):
         if index_LocY[0]-s-1 >= 0 and pdfY[index_LocY[0]-s] >= (Sum*win.th_pro):
          Index_errorY[0]=index_LocY[0]-s-1
         if index_LocY[0]+s+1 < len(pdfY) and pdfY[index_LocY[0]+s] >= (Sum*win.th_pro):
          Index_errorY[1]=index_LocY[0]+s+1
        Sum=pdfZ[index_LocZ[0]]
        Index_errorZ=[0,len(pdfZ)-1]		 
        for s in range(0,len(pdfZ)):
         if index_LocZ[0]-s-1 >= 0 and pdfZ[index_LocZ[0]-s] >= (Sum*win.th_pro):
          Index_errorZ[0]=index_LocZ[0]-s-1
         if index_LocZ[0]+s+1 < len(pdfZ) and pdfZ[index_LocZ[0]+s] >= (Sum*win.th_pro):
          Index_errorZ[1]=index_LocZ[0]+s+1
        ab.plot(win.xnode[Index_errorX[0],:,0],win.ynode[Index_errorX[0],:,0],'k--')
        ab.plot(win.xnode[Index_errorX[1],:,0],win.ynode[Index_errorX[1],:,0],'k--')		
        ab.plot(win.xnode[:,Index_errorY[0],0],win.ynode[:,Index_errorY[0],0],'k--')
        ab.plot(win.xnode[:,Index_errorY[1],0],win.ynode[:,Index_errorY[1],0],'k--')	
		
###########################################################################################################################################

        ac = self.figure.add_subplot(spec[2-1,3-1])
        ac.yaxis.set_ticks_position("right")  

        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        if len(index)>0:
         for i in index:
          ac.plot([0],[yy[i]]-Locy,'k^')
	
        ac.plot(-win.znode[index_LocX[0],:,Index_errorZ[0]],win.ynode[index_LocX[0],:,Index_errorZ[0]]-Locy,'k--')
        ac.plot(-win.znode[index_LocX[0],:,Index_errorZ[1]],win.ynode[index_LocX[0],:,Index_errorZ[1]]-Locy,'k--')
        if len(xev)>0:
         ac.plot([zev[0]],yev-Locy,'*',markersize=8,markerfacecolor=(1, 1,1, 1),markeredgecolor='k')
        ac.plot(-Locz,Locy-Locy,'*',markersize=8,markerfacecolor=(0, 0,0, 1),markeredgecolor='k')
        cset=ac.contourf(-win.znode[index_LocX[0],:,:],win.ynode[index_LocX[0],:,:]-Locy,win.Locmatrix[index_LocX[0],:,:]/np.max(win.Locmatrix),10, cmap='hot_r',vmin = 0, vmax = 1, alpha=.5, antialiased=True)
        ac.set_xlabel('Z(m)',fontsize=8)
        ac.set_ylabel('N-S(m)',fontsize=8)

        lock_log.acquire()	
        os.chdir(Path_grid)
        os.chdir(Path_OUT)
        win.file_log=open('MParLoc_event.log', 'a')
        Loninf,LatInf=map(win.xnode[Index_errorX[0],0,0],win.ynode[0,Index_errorY[0],0],inverse=True)
        Lonsup,Latsup=map(win.xnode[Index_errorX[1],0,0],win.ynode[0,Index_errorY[1],0],inverse=True)
        win.file_log.write(str(data[0].stats.endtime+0.1)+': NEW LOCATION LAT: '+str(Latev)+'['+str(LatInf)+' , '+str(Latsup)+']'+'\n') 
        win.file_log.write(str(data[0].stats.endtime+0.1)+': NEW LOCATION LON: '+str(Lonev)+'['+str(Loninf)+' , '+str(Lonsup)+']'+'\n')  
        win.file_log.write(str(data[0].stats.endtime+0.1)+': NEW LOCATION DEP(km): '+str(-win.znode[index_LocY[0],index_LocX[0],index_LocZ[0]]/1000)+'['+str(-win.znode[index_LocX[0],index_LocY[0],Index_errorZ[0]]/1000)+' , '+str(-win.znode[index_LocX[0],index_LocY[0],Index_errorZ[1]]/1000)+']'+'\n')  
        win.file_log.write('Error_E-W: '+str(np.abs(win.xnode[Index_errorX[0],0,0]-win.xnode[Index_errorX[1],0,0])/2000)+' Km\n') 
        win.file_log.write('Error_N-S: '+str(np.abs(win.ynode[0,Index_errorY[0],0]-win.ynode[0,Index_errorY[1],0])/2000)+' Km\n') 
        win.file_log.write('Error_Z: '+str(np.abs(win.znode[0,0,Index_errorZ[0]]-win.znode[0,0,Index_errorZ[1]])/2000)+' Km\n') 
        Er_epi=np.sqrt(((win.xnode[Index_errorX[0],0,0]-win.xnode[Index_errorX[1],0,0])/2000)**2+((win.ynode[0,Index_errorY[0],0]-win.ynode[0,Index_errorY[1],0])/2000)**2)
        Er_Z=np.abs(win.znode[0,0,Index_errorZ[0]]-win.znode[0,0,Index_errorZ[1]])/2000
        win.file_log.close()
        os.chdir(Path_grid)
        lock_log.release()	
		

		
######################################################################################################
        ad = self.figure.add_subplot(spec[3-1,2-1])
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        if len(index)>0:
         for i in index:
          ad.plot([xx[i]]-Locx,[0],'k^')
        ad.plot(win.xnode[:,index_LocY[0],:]-Locx,win.znode[:,index_LocY[0],Index_errorZ[0]],'k--')
        ad.plot(win.xnode[:,index_LocY[0],:]-Locx,win.znode[:,index_LocY[0],Index_errorZ[1]],'k--')	
        cset=ad.contourf(win.xnode[:,index_LocY[0],:]-Locx, win.znode[:,index_LocY[0],:], win.Locmatrix[:,index_LocY[0],:]/np.max(win.Locmatrix),10, cmap='hot_r',vmin = 0, vmax = 1, alpha=.5, antialiased=True)
        if len(xev)>0:
         ad.plot(xev-Locx,zev,'*',markersize=8,markerfacecolor=(1, 1,1, 1),markeredgecolor='k')
        ad.plot(Locx-Locx,Locz,'k*',markersize=8)
        ad.set_ylabel('Z(m)',fontsize=8)
        ad.set_xlabel('E-W(m)',fontsize=8)
		
##########################################################################################################################
		
        ag = self.figure.add_subplot(spec[2-1,1-1])
        ag.plot(pdfY/np.max(pdfY),win.ynode[Index_errorX[0],:,index_LocZ[0]]-Locy,'k')

        ag.plot(np.linspace(0,1,len(win.xnode[Index_errorX[0],:,index_LocZ[0]]-Locx)),win.xnode[Index_errorX[0],:,index_LocZ[0]]-Locx,'k--')
        ag.plot(np.linspace(0,1,len(win.xnode[Index_errorX[1],:,index_LocZ[0]]-Locx)),win.xnode[Index_errorX[1],:,index_LocZ[0]]-Locx,'k--')
        ag.set_ylim([np.min(win.ynode[Index_errorX[0],:,index_LocZ[0]]-Locy), np.max(win.ynode[Index_errorX[0],:,index_LocZ[0]]-Locy)])
        ag.set_ylabel('N-S(m)',fontsize=8)
        ag.xaxis.set_label_position("top")
        ag.xaxis.set_ticks_position("top")
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        ag.set_xlabel('Norm. Prob.',fontsize=8)
        ag.xaxis.set_major_formatter(FormatStrFormatter('%.g'))
		
###########################################################################################################################
        ah = self.figure.add_subplot(spec[1-1,2-1])
        ah.plot(win.xnode[:,index_LocY[0],index_LocZ[0]]-Locx,pdfX/np.max(pdfX),'k')
        ah.plot(win.ynode[:,Index_errorY[0],index_LocZ[0]]-Locy,np.linspace(0,1,len(win.ynode[:,Index_errorY[0],index_LocZ[0]])),'k--')
        ah.plot(win.ynode[:,Index_errorY[1],index_LocZ[0]]-Locy,np.linspace(0,1,len(win.ynode[:,Index_errorY[1],index_LocZ[0]])),'k--')
        ah.set_xlim([np.min(win.xnode[:,index_LocY[0],index_LocZ[0]]-Locx), np.max(win.xnode[:,index_LocY[0],index_LocZ[0]]-Locx)])
        ah.xaxis.set_label_position("top")
        ah.xaxis.set_ticks_position("top")
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        ah.set_xlabel('E-W(m)',fontsize=8)
        ah.set_ylabel('Norm. Prob.',fontsize=8)	
        ah.yaxis.set_major_formatter(FormatStrFormatter('%.g'))
		
		
################################################################################################################################
        ai = self.figure.add_subplot(spec[3-1,1-1])
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)
        ai.plot(pdfZ/np.max(pdfZ),win.znode[index_LocX[0],index_LocY[0],:],'k')
        ai.plot(np.linspace(0,1,len(win.znode[:,index_LocY[0],Index_errorZ[0]])),win.znode[:,index_LocY[0],Index_errorZ[0]],'k--')
        ai.plot(np.linspace(0,1,len(win.znode[:,index_LocY[0],Index_errorZ[0]])),win.znode[:,index_LocY[0],Index_errorZ[1]],'k--')
        ai.set_ylabel('Z(m)',fontsize=8)
        ai.set_xlabel('Norm. Prob.',fontsize=6)	
        ai.xaxis.set_major_formatter(FormatStrFormatter('%.g'))
        if True:
         os.chdir(Path_OUT)
         self.figure.savefig('MParLoc_event_FIGURE.pdf',format='pdf')
         self.file_log=open('MParLoc_event.h71', 'w')
         self.file_log.write('  DATE    ORIGIN    LAT      LONG      DEPTH    MAG NO DM GAP M  RMS  ERH  ERZ Q SQD  ADJ IN NR  AVR  AAR NM AVXM SDXM NF AVFM SDFM I \n')
         Num_obs=len(np.where(np.array(pick)!=0)[0])+len(np.where(np.array(pickS)!=0)[0])
         Nomi=[]
         dist=[]
         picchiS=[]
         picchiP=[]
         Azs=[]
         AINs=[]
         Indexs=[]
         for i in range(0,len(NAME)):
          if pick[i]!=0 or pickS[i]!=0:
           Nomi.append(NAME[i])
           Indexs.append(i)
           info=gps2dist_azimuth(Latev,Lonev,staz_lat[i],staz_lon[i])
           Azs.append(np.round((info[2])))
           win.file_log=open('MParLoc_event.log', 'a')
           if Azs[-1]+180>360:
            win.file_log.write('BAZ STAZ: '+str(NAME[i])+': '+str(Azs[-1]+180-360)+' BAZ REALTIME: '+str(Baz[i][0])+'\n') 
           else:
            win.file_log.write('BAZ STAZ: '+str(NAME[i])+': '+str(Azs[-1]+180)+' BAZ REALTIME: '+str(Baz[i][0])+'\n') 	   
           win.file_log.close()		
           AINs.append((AIN[i]))
           dist.append(np.sqrt((info[0]/1000)**2+(-win.znode[index_LocX[0],index_LocY[0],index_LocZ[0]]/1000)**2))
         GAP=int(360-(np.max(Azs)-np.min(Azs)))	 
         Azsorted=np.sort(Azs.copy())	
         gap=0
         for i in range(0,len(Azsorted)-1):
          if gap!=0:
           if gap<Azsorted[i+1]-Azsorted[i]:
            gap=Azsorted[i+1]-Azsorted[i]
          else:
           gap=Azsorted[i+1]-Azsorted[i]
         if gap>GAP:		   
          GAP=gap
         LON=str(Lonev[0]).rsplit('.')[0]+' '+str(round(60*60*(Lonev[0]-float(str(Lonev[0]).rsplit('.')[0]))/100,2))
         LAT=str(Latev[0]).rsplit('.')[0]+' '+str(round(60*60*(Latev[0]-float(str(Latev[0]).rsplit('.')[0]))/100,2))
         self.file_log.write(' '+str(To)[2:10].replace('-','')+' '+str(To)[11:16].replace(':','')+' '+str(To)[17:22]+'  '+LAT+'  '+LON+'  '+str(round(-win.znode[index_LocX[0],index_LocY[0],index_LocZ[0]]/1000,2))+'   0.00  ')
         self.file_log.write(str(Num_obs)+' 0 '+str(GAP)+' 0 '+' '+str(np.round((np.sum(np.abs(RES_P))+np.sum(np.abs(RES_S)))/Num_obs,2))+' '+str(round(Er_epi,2))+' '+str(round(Er_Z,2))+' 0 0 0 0.00  0  0 0.00 0.00  0  0.0  0.0  0  0.0  0.0 0\n\n')		 
         self.file_log.write('  STN  DIST AZM AIN PRMK HRMN P-SEC TPOBS TPCAL DLY/H1 P-RES P-WT AMX PRX CALX K XMAG RMK FMP FMAG SRMK S-SEC TSOBS S-RES  S-WT    DT\n')
         Index_sorted=np.argsort(dist)
         for i in Index_sorted:
          if pick[Indexs[i]]!=0:
           if dist[i]>=100:
            self.file_log.write(' '+Nomi[i]+' '+str('%.1f' % (dist[i])))
           if dist[i]>=10 and dist[i]<100:
            self.file_log.write(' '+Nomi[i]+'  '+str('%.1f' % (dist[i])))
           if dist[i]<10:
            self.file_log.write(' '+Nomi[i]+'   '+str('%.1f' % (dist[i])))
           if Azs[i]<10:
            self.file_log.write('   '+str(int(Azs[i])))  
           if Azs[i]>=10 and Azs[i]<100:
            self.file_log.write('  '+str(int(Azs[i]))) 
           if Azs[i]>=100:
            self.file_log.write(' '+str(int(Azs[i])))
           if AINs[i]<10:
             self.file_log.write('   '+str(int(AINs[i]))+'  P ? ')  
           if AINs[i]>=10 and AINs[i]<100:
             self.file_log.write('  '+str(int(AINs[i]))+'  P ? ') 
           if AINs[i]>=100:
            self.file_log.write(' '+str(int(AINs[i]))+'  P ? ')
           self.file_log.write(str(pick[Indexs[i]])[11:16].replace(':',''))
           if float(str(pick[Indexs[i]])[17:22])<10:
            self.file_log.write('  '+str('%.2f' % (float(str(pick[Indexs[i]])[17:22]))))
           if float(str(pick[Indexs[i]])[17:22])>=10:	
            self.file_log.write(' '+str('%.2f' % (float(str(pick[Indexs[i]])[17:22]))))		
           if float(pick[Indexs[i]]-To)<10:
            self.file_log.write('  '+str('%.2f' % float((pick[Indexs[i]])-To)))
           if float(pick[Indexs[i]]-To)>=10:	
            self.file_log.write(' '+str('%.2f' % float((pick[Indexs[i]])-To)))
           if float(P_teo[Indexs[i]]-To)<10:
            self.file_log.write('  '+str('%.2f' % float((P_teo[Indexs[i]])-To))+'  0.0')
           if float(P_teo[Indexs[i]]-To)>=10:	
            self.file_log.write(' '+str('%.2f' % float((P_teo[Indexs[i]])-To))+'  0.0')
           if float(RES_P[Indexs[i]])<10 and float(RES_P[Indexs[i]])>=0:
            self.file_log.write('  '+str('%.2f' % float(RES_P[Indexs[i]])))
           if float(RES_P[Indexs[i]])>=10 or float(RES_P[Indexs[i]])<0:
            self.file_log.write(' '+str('%.2f' % float((RES_P[Indexs[i]]))))
           snr=(float(SNR[Indexs[i]]))
           if (1/(2*snr))<10:
            self.file_log.write('  '+str('%.2f' % (1/(2*snr)))+' ')
           if (1/(2*snr))>=10:	
            self.file_log.write(' '+str('%.2f' % (1/(2*snr)))+' ')
           self.file_log.write(' 0.0 0.0 0.00 0 0.00 000 0.0 0.00 0.00  0.00  0.00 00.00   0.0   0.0\n')
          if pickS[Indexs[i]]!=0:	
           if dist[i]>=100:		  
            self.file_log.write(' '+Nomi[i]+' '+str('%.1f' % (dist[i])))
           if dist[i]>=10 and dist[i]<100:
            self.file_log.write(' '+Nomi[i]+'  '+str('%.1f' % (dist[i])))
           if dist[i]<10:
            self.file_log.write(' '+Nomi[i]+'   '+str('%.1f' % (dist[i])))	
           if Azs[i]<10:
            self.file_log.write('   '+str(int(Azs[i])))  
           if Azs[i]>=10 and Azs[i]<100:
            self.file_log.write('  '+str(int(Azs[i]))) 
           if Azs[i]>=100:
            self.file_log.write(' '+str(int(Azs[i])))
           if AINs[i]<10:
             self.file_log.write('   '+str(int(AINs[i]))+'  S ? ')  
           if AINs[i]>=10 and AINs[i]<100:
             self.file_log.write('  '+str(int(AINs[i]))+'  S ? ') 
           if AINs[i]>=100:
             self.file_log.write(' '+str(int(AINs[i]))+'  S ? ')
           self.file_log.write(str(pickS[Indexs[i]])[11:16].replace(':',''))
           if float(str(pickS[Indexs[i]])[17:22])<10:
            self.file_log.write('  '+str('%.2f' % (float(str(pickS[Indexs[i]])[17:22]))))
           if float(str(pickS[Indexs[i]])[17:22])>=10:	
            self.file_log.write(' '+str('%.2f' % (float(str(pickS[Indexs[i]])[17:22]))))
           if float(pickS[Indexs[i]]-To)<10:
            self.file_log.write('  '+str('%.2f' % float((pickS[Indexs[i]])-To)))
           if float(pickS[Indexs[i]]-To)>=10:	
            self.file_log.write(' '+str('%.2f' % float((pickS[Indexs[i]])-To)))
           if float(S_teo[Indexs[i]]-To)<10:
            self.file_log.write('  '+str('%.2f' % float((S_teo[Indexs[i]])-To))+'  0.0')
           if float(S_teo[Indexs[i]]-To)>=10:	
            self.file_log.write(' '+str('%.2f' % float((S_teo[Indexs[i]])-To))+'  0.0')
           if float(RES_S[Indexs[i]])<10 and float(RES_S[Indexs[i]])>=0:
            self.file_log.write('  '+str('%.2f' % float(RES_S[Indexs[i]])))
           if float(RES_S[Indexs[i]])>=10 or float(RES_S[Indexs[i]])<0:	
            self.file_log.write(' '+str('%.2f' % float((RES_S[Indexs[i]]))))
           snr=(float(SNR_S[Indexs[i]]))
           if (1/(2*snr))<10:
            self.file_log.write('  '+str('%.2f' % (1/(2*snr)))+' ')
           if (1/(2*snr))>=10:	
            self.file_log.write(' '+str('%.2f' % (1/(2*snr)))+' ')
           self.file_log.write(' 0.0 0.0 0.00 0 0.00 000 0.0 0.00 0.00  0.00  0.00 00.00   0.0   0.0\n')
         self.file_log.close()		 
         os._exit(0)
 
 
 
def BAZ(Trace_orig,accE,accN,staz_index):

#################### DATA processing   ###############################
 AccZ=(((((Trace_orig.slice(pick[staz_index]-3,pick[staz_index]+1)))))).detrend('constant').filter("bandpass",freqmin=1,freqmax=4)
 AccE=(((((accE.slice(pick[staz_index]-3,pick[staz_index]+1)))))).detrend('constant').filter("bandpass",freqmin=1,freqmax=4)
 AccN=(((((accN.slice(pick[staz_index]-3,pick[staz_index]+1)))))).detrend('constant').filter("bandpass",freqmin=1,freqmax=4)
 
 ################### Pre-signal noise-level estimation################
 
 preSignalZ=np.max(np.abs(AccZ.slice(pick[staz_index]-3,pick[staz_index])))
 preSignalN=(np.max((np.abs(AccN.slice(pick[staz_index]-3,pick[staz_index]).data))))
 preSignalE=(np.max((np.abs(AccE.slice(pick[staz_index]-3,pick[staz_index]).data))))
 
 AccZ=AccZ.slice(pick[staz_index],pick[staz_index]+1)
 AccE=AccE.slice(pick[staz_index],pick[staz_index]+1)
 AccN=AccN.slice(pick[staz_index],pick[staz_index]+1)
 
 ################### Use of data until the Maximum of Amplitude ##################

 #AccE.data=AccE.data[0:np.where(np.abs(AccZ.data)==np.max(np.abs(AccZ.data)))[0][0]]
 #AccN.data=AccN.data[0:np.where(np.abs(AccZ.data)==np.max(np.abs(AccZ.data)))[0][0]]
 
 BAZ=np.zeros(shape=(len(AccZ)))-1
 ###################  Threshold for the data usability ##########################
 SNR_max=3
 ##################  Estimation of BAZ for all samples   #########################
 
 ######### iteration == 0  #####################################################
 RZE=0
 RZN=0
 ################################################################### 
 for i in range(0,min([len(AccZ),len(AccE),len(AccN)])): 
 ######### iteration > 1  #####################################################
  if (np.abs(AccN[i])/preSignalN)>SNR_max or (np.abs(AccE[i])/preSignalE)>SNR_max:
   RZE=(0.99*RZE+AccE[i]*AccZ[i])
   RZN=(0.99*RZN+AccN[i]*AccZ[i])
   BAZp=math.degrees(math.atan2(RZE,RZN)+math.pi)
   if BAZp>360:
    BAZp=BAZp-math.floor(BAZp/360)*360
   BAZ[i]=BAZp

####################  Final Computation of median value ####################
 if len(BAZ[np.where(BAZ!=-1)[0]])>0:
  BAZ=BAZ[np.where(BAZ!=-1)[0]]  
  x,y=np.histogram(BAZ, bins=[0,30,60,90,120,150,180,210,240,270,300,330,360])
  BazMax=y[np.where((np.array(x))==np.max(x))[0][0]]+60
  BazMin=y[np.where((np.array(x))==np.max(x))[0][0]]-60
  
  if BazMin<0:
   BazMin=BazMin+360
   BazMax=BazMax+360
   BAZ[np.where(BAZ<180)[0]]=BAZ[np.where(BAZ<180)[0]]+360
  Baz_sum=[]
  for count in BAZ:
   if count >=BazMin and count <=BazMax:
    Baz_sum.append(count)
  if len(Baz_sum)>0:
   BAZ=np.mean((Baz_sum))
   BAZ_error=np.std((Baz_sum))
  else:
   BAZ=0
   BAZ_error=0
 else:
  BAZ=0
  BAZ_error=0
 if BAZ>360:
  BAZ=BAZ-360
 return(BAZ,BAZ_error)
 



class MainCode():
         def simulation(self,start_time,staz_index):
          global data 
          global dataY
          global dataX 
          global AccZ
          global AccN
          global AccE
          global pick
          global Baz
          global change_factor
          global new_sum
          self.frist_pick=0
          change_factor=False
          new_sum=False
          index=0
          while True:
           pakage=self.tr[staz_index].copy()
           pakageX=self.trX[staz_index].copy()
           pakageY=self.trY[staz_index].copy()
           if len(pakage.data)>0:
            if len(data[staz_index])==0:
             data[staz_index]=pakage.copy()
             dataX[staz_index]=pakageX.copy()
             dataY[staz_index]=pakageY.copy() 
            else:
             data[staz_index]=self.tr[staz_index].copy()
             dataX[staz_index]=self.trX[staz_index].copy()
             dataY[staz_index]=self.trY[staz_index].copy()
           time.sleep(0.1)
         def __init__(self, parent=None):
          global model
          global Model_path_high
          self.sensor=[]
          self.logger=[]
          os.chdir(Network_path)
          Model_path_low=os.path.abspath('grids')
          Model_path_high=os.path.abspath('grids')
          taup_create.build_taup_model(Model_path_low+r'/velocity_model.tvel',Model_path_low )
          model = TauPyModel(Model_path_low+r'/velocity_model.npz')
          global map
          global ListaZ
          global ListaE
          global ListaN
          global Path_OUT
          global NodeX
          global NodeY
          global NodeZ
          global RES_P
          global RES_S
          global RMS
          global Path_grid
          global grd
          config_file=open('../config_files/run_loc.conf')
          config_file=config_file.read()		  
          self.th_pro=float(config_file.rsplit('Confidence_probability_threshold_level:')[1].rsplit('\n')[0])
          grd=glob.glob(Model_path_low+r'/*buf')
          grd = NLLGrid(grd[0])
          Path_grid=Model_path_low
          self.grid_precision=str(config_file.rsplit('grid_precision: ')[1].rsplit('\n')[0])
          if 'f' in self.grid_precision:
           self.grid_precision='f'
          if 'd' in self.grid_precision:
           self.grid_precision='d'
          self.min_numP=int(config_file.rsplit('no_min_p: ')[1].rsplit('\n')[0])	
          self.min_numS=int(config_file.rsplit('no_min_s: ')[1].rsplit('\n')[0])			
          self.snr_wind_S=float(config_file.rsplit('snr_wind_s: ')[1].rsplit('\n')[0])		 
          self.snr_wind_p=float(config_file.rsplit('snr_wind_p: ')[1].rsplit('\n')[0])	
          self.diff_Amp=str(config_file.rsplit('diff_Amp: ')[1].rsplit('\n')[0])
          self.diff_Amp_b=float(config_file.rsplit('diff_Amp_b: ')[1].rsplit('\n')[0])
          self.diff_Amp_error=float(config_file.rsplit('diff_Amp_error: ')[1].rsplit('\n')[0])
          self.back_az_error=float(config_file.rsplit('back_az_error: ')[1].rsplit('\n')[0])
          if 'False' in self.diff_Amp:
            self.diff_Amp=False	
          else:		
            self.diff_Amp=True		  
          self.back_az=str(config_file.rsplit('back_az: ')[1].rsplit('\n')[0])
          if 'False' in self.back_az:
            self.back_az=False	
          else:		
            self.back_az=True			  
          self.s_phases=str(config_file.rsplit('s_phases: ')[1].rsplit('\n')[0])	
          if 'False' in self.s_phases:
            self.s_phases=False	
          else:		
            self.s_phases=True		  
          Tabella_Pesi=open('../config_files/weights_table.dat')
          Tabella_Pesi=Tabella_Pesi.read()  
          Tabella_Pesi=Tabella_Pesi.rsplit('\n')
          self.Tabella_P=[[2,0]]
          self.Tabella_S=[[2,0]]
          for count in Tabella_Pesi:
           if len(count.rsplit(' '))==4:
            self.Tabella_P.append([float(count.rsplit(' ')[1]),float(count.rsplit(' ')[2])])
            self.Tabella_S.append([float(count.rsplit(' ')[1]),float(count.rsplit(' ')[3])])
           if len(count.rsplit(' '))==3:
            self.Tabella_P.append([float(count.rsplit(' ')[0]),float(count.rsplit(' ')[1])])
            self.Tabella_S.append([float(count.rsplit(' ')[0]),float(count.rsplit(' ')[2])])
          os.chdir(Event_path)
          file=sorted(glob.glob('MParLoc_event*'))	  
          if len(file)>0:
           for count in file:
            os.remove(count)
          Lista_sac=sorted(glob.glob('*'))	  
          self.tr=[]
          self.trY=[]
          self.trX=[]
          global NAME
          NAME=[]
          Nomi=open('../list_stations_4loc.dat')
          Nomi=Nomi.read()
          Nomi=Nomi.rsplit('\n')
          self.index=0
          for i in range(0,len(Lista_sac)):
           if 'Z' in read(Lista_sac[i])[0].stats.sac.kcmpnm or 'U' in read(Lista_sac[i])[0].stats.sac.kcmpnm:
            if read(Lista_sac[i])[0].stats.station in Nomi:
             NAME.append(read(Lista_sac[i])[0].stats.station)
             self.tr.append(read(Lista_sac[i])[0])
             for j in range(0,len(Lista_sac)):
              if read(Lista_sac[j])[0].stats.station == NAME[-1]:
               if 'N' in read(Lista_sac[j])[0].stats.sac.kcmpnm:
                self.trY.append(read(Lista_sac[j])[0])
               if 'E' in read(Lista_sac[j])[0].stats.sac.kcmpnm:			
                self.trX.append(read(Lista_sac[j])[0])			 
          global RU_sup	
          global xx
          global yy
          global xev
          global yev
          global zev
          global trigger
          global P_teo
          global S_teo
          global Real_BAz
          global latev
          global lonev
          global evdp
          global staz_lon
          global staz_lat
          global iter	
          global Xo	
          global Yo		  
          global Zo		  		  
          iter=[]
          P_teo=[]
          S_teo=[]
          staz_lon=[]
          staz_lat=[]
          Real_BAz=[]
          xev=[]
          yev=[]
          zev=[]
          dX=grd.dx#km
          dY=grd.dy #km
          dZ=grd.dz #km
          NodeX=grd.nx#km
          NodeY=grd.ny #km
          NodeZ=grd.nz #km
          Xo=grd.x_orig #km
          Yo=grd.y_orig #km
          Zo=grd.z_orig  #km
          Lat0=grd.orig_lat
          Lon0= grd.orig_lon
          lambert_area = {'proj': 'laea',
          'lat_0':Lat0, 
          'lon_0':Lon0, 
          'x_0':Xo, 
          'y_0':Yo,
          'ellps': 'WGS84',
          'datum': 'WGS84',
          'R':6378137.0}
          map = pyproj.Proj(lambert_area)
          for i in range(0,len(self.tr)):
           iter.append(0)
           staz_lon.append(self.tr[i].stats.sac.stlo)
           staz_lat.append(self.tr[i].stats.sac.stla)
          xx,yy=map(staz_lon,staz_lat)			  		  
          global Baz
          Baz=[]
          for i in range(0,len(self.tr)):
           Baz.append([-9999,-9999,-9999])
          global pick
          global SNR
          global pickS
          pick=[]
          SNR=[]
          pickS=[]
          for i in range(0,len(self.tr)):
           self.logger.append(self.tr[i].stats.sac.user2)
           self.sensor.append(self.tr[i].stats.sac.user3)
           if self.s_phases:
            if hasattr(self.trY[i].stats.sac, 't0') or hasattr(self.trX[i].stats.sac, 't0'):
             if hasattr(self.trX[i].stats.sac, 't0'):
              pickS.append(self.trX[i].stats.starttime-self.trX[i].stats.sac.b+self.trX[i].stats.sac.t0)
             else:
              pickS.append(self.trY[i].stats.starttime-self.trY[i].stats.sac.b+self.trY[i].stats.sac.t0)
            else:
              pickS.append(0)
           else:
             pickS.append(0)
           if hasattr(self.tr[i].stats.sac, 'a'):
            if self.tr[i].stats.sac.a!=-1:
             pick.append(self.tr[i].stats.starttime-self.tr[i].stats.sac.b+self.tr[i].stats.sac.a)
            else:
             pick.append(0)
           else:
            pick.append(0)
          global data 
          global dataY
          global dataX
          if len(np.where(np.array(pick)!=0)[0])<self.min_numP and len(np.where(np.array(pickS)!=0)[0])<self.min_numS:
           os._exit(0)
          SNR=(np.ones(shape=(len(NAME))))
          data=[]
          dataX=[]
          dataY=[]
          global AccZ
          global AccN
          global AccE
          RES_P=[]
          RES_S=[]
          AccE=[]
          AccN=[]
          AccZ=[]  
          self.Pv=[]	  
          for i in range(0,len(NAME)):
           RES_P.append(0)
           RES_S.append(0)
           self.Pv.append(0)
           data.append(Trace())
           dataX.append(Trace())
           dataY.append(Trace())
           AccE.append(Trace())
           AccN.append(Trace())
           AccZ.append(Trace())
          grd=glob.glob(Model_path_high+r'/*buf')
          grd = NLLGrid(grd[0])
          NodeXh=grd.nx#km
          NodeYh=grd.ny #km
          NodeZh=grd.nz #km
          self.ynode, self.xnode, self.znode = np.meshgrid(np.linspace(Yo*1000,(Yo+NodeY*dY-dY)*1000,NodeY),np.linspace(Xo*1000,(Xo+NodeX*dX-dX)*1000,NodeX), -np.linspace(Zo*1000,(Zo+NodeZ*dZ-dZ)*1000,NodeZ))
          self.BAznode=[]
          grd=glob.glob(Model_path_low+r'/*buf')
          grd = NLLGrid(grd[0])
          self.Rnode=[]
          self.TarrivalnodeErr=[]
          if self.back_az:
           self.BAzsumErr=(np.zeros(shape=(NodeX,NodeY,NodeZ)))
          if self.diff_Amp:
           self.RsumErr=(np.zeros(shape=(NodeX,NodeY,NodeZ)))
          self.TsumErr=(np.ones(shape=(NodeX,NodeY,NodeZ)))
          self.TsumErr_S=(np.ones(shape=(NodeX,NodeY,NodeZ)))
          self.TsumErr_S_P=(np.ones(shape=(NodeX,NodeY,NodeZ)))
          self.pick=0
          for i in range(0,len(self.tr)):
           if i==0:
            start_time=self.tr[i].stats.starttime
           else:
            if start_time-self.tr[i].stats.starttime>0:
             start_time=self.tr[i].stats.starttime
          thread_simulation=[]
          Path_OUT=os.path.abspath('')
          os.chdir(r'../..')           		  
          for i in range(0,len(self.tr)):
           thread_simulation.append(threading.Thread(target=partial(self.simulation,start_time,i)))
           thread_simulation[i].start()
          thread_processing=threading.Thread(target=partial(self.Final_loc))
          thread_processing.start()
          thread_sum=threading.Thread(target=self.Sum_matrix)	 
          thread_sum.start()		  
          self.mappa=Location_MAP()         		  

         def Sum_matrix(self):
          global new_sum
          global Locx
          global Locy
          global Locz
          global MinX
          global MinY
          global MinZ
          global index_LocZ
          global index_LocX
          global index_LocY
          global change_factor
          global P_teo
          global S_teo
          global Scart_P
          global RES_P
          global RES_S		  
          global SNR_P
          global SNR_SX
          global SNR_SY   
          global To
          global AIN
          global Path_grid
          global Xindexs
          global Yindexs
          global Zindexs
          global Lonev
          global Latev
          global grd		  	 
          global NodeX
          global NodeY
          global NodeZ
          second_loc=1
          while True:
           if new_sum==True:
            if len(np.where(np.array(pick)!=0)[0])==0:
             print('No picks avaliable for low SNR!')
             os._exit(0)
            lock_log.acquire()
            os.chdir(Path_OUT)
            self.file_log=open('MParLoc_event.log', 'a')
            self.file_log.write('##################################################################################################################'+'\n')
            self.file_log.write('##################################################################################################################'+'\n')
            for staz_index in range(0,len(NAME)):
             self.file_log.write('Staz: '+NAME[staz_index]+' pickP:'+str(pick[staz_index])+' pickS: '+str(pickS[staz_index])+'\n')  
            self.file_log.close()
            os.chdir(Path_grid)
            lock_log.release()
            if self.back_az:
             if np.sum((self.BAzsumErr))!=0:
                BAzsumErr=(self.BAzsumErr.copy())
                BAzsumErr=BAzsumErr/np.sum((BAzsumErr))
             else:
                BAzsumErr=1
            else:
                BAzsumErr=1
            if self.diff_Amp:
             if np.sum((self.RsumErr))!=0:
                RsumErr=(self.RsumErr.copy())
                RsumErr=RsumErr/np.sum((RsumErr))
             else:
                RsumErr=1
            else:
                RsumErr=1
            self.TsumErr=np.where(self.TsumErr==1, 0, self.TsumErr)
            if np.sum((self.TsumErr))!=0:
             TsumErr=(self.TsumErr.copy())/np.sum((self.TsumErr))
            else:
             TsumErr=1
            self.TsumErr=np.where(self.TsumErr==0, 1, self.TsumErr)
            self.TsumErr_S_P=np.where(self.TsumErr_S_P==1, 0, self.TsumErr_S_P)			 
            if np.sum((self.TsumErr_S_P))!=0:
             TsumErr_S_P=(self.TsumErr_S_P.copy())/np.sum((self.TsumErr_S_P))
            else:
             TsumErr_S_P=1
            self.TsumErr_S_P=np.where(self.TsumErr_S_P==0, 1, self.TsumErr_S_P)			 
            self.TsumErr_S=np.where(self.TsumErr_S==1, 0, self.TsumErr_S)
            if np.sum((self.TsumErr_S))!=0:
             TsumErr_S=(self.TsumErr_S.copy())/np.sum((self.TsumErr_S))
            else:
             TsumErr_S=1
            self.TsumErr_S=np.where(self.TsumErr_S==0, 1, self.TsumErr_S)
            if self.s_phases:
              self.Locmatrix=TsumErr*TsumErr_S_P*TsumErr_S*RsumErr*BAzsumErr
            else:
              self.Locmatrix=TsumErr*RsumErr*BAzsumErr			 
            self.Locmatrix=self.Locmatrix
            LocX, LocY, LocZ=(np.where(self.Locmatrix==np.max(self.Locmatrix))[0:3])
            P_teo_p=[]
            S_teo_p=[]
            SNR_P=[]
            SNR_SX=[]
            SNR_SY=[]           
            AIN=[]
            if len(LocX)==1:
             index_LocZ=LocZ
             index_LocY=LocY
             index_LocX=LocX
             Locx=self.xnode[LocX,LocY,LocZ]
             Locy=self.ynode[LocX,LocY,LocZ]			
             Locz=self.znode[LocX,LocY,LocZ]
             Lonev,Latev=map(Locx,Locy,inverse=True)
             for i in range(0,len(NAME)):
               os.chdir(Path_grid)
               f = open(Model_path_high+r'/arrival_matrix'+NAME[i]+'.bin', "rb")
               floats = array(self.grid_precision)
               floats.fromfile(f,NodeX*NodeY*NodeZ)
               floats = np.reshape(floats,(NodeX,NodeY,NodeZ))
               self.arrival_matrix=(floats)
               f = open(Model_path_high+r'/arrival_matrixS'+NAME[i]+'.bin', "rb")
               floats = array(self.grid_precision)
               floats.fromfile(f,NodeX*NodeY*NodeZ)
               floats = np.reshape(floats,(NodeX,NodeY,NodeZ))
               self.arrival_matrix_S=(floats)
               dist=calc_dist(Latev,Lonev,staz_lat[i],staz_lon[i],flattening_of_planet=0,radius_of_planet_in_km=6372.797)
               evdp=-self.znode[LocX,LocY,LocZ]/1000
               if evdp < 0:
                deptau=[0]
               else:
                deptau=evdp
               arrivals = model.get_travel_times(source_depth_in_km=deptau[0], distance_in_degree=dist,phase_list=('tts','ttS')) 
               S_teo_p.append(self.arrival_matrix_S[LocX,LocY,LocZ][0])
               arrivals = model.get_travel_times(source_depth_in_km=deptau[0], distance_in_degree=dist,phase_list=('ttp','ttP')) 
               AIN.append(arrivals[0].incident_angle)
               P_teo_p.append(self.arrival_matrix[LocX,LocY,LocZ][0])	
               os.chdir(Path_OUT)
             if second_loc==0:
              Path_grid=Model_path_high
              grd=glob.glob(Model_path_high+r'/*buf')
              grd = NLLGrid(grd[0])
              NodeX=grd.nx#km
              NodeY=grd.ny #km
              NodeZ=grd.nz #km
              dX=grd.dx#km
              dY=grd.dy #km
              dZ=grd.dz #km
              self.xnode, self.ynode, self.znode = np.meshgrid(np.linspace(Xo*1000,(Xo+NodeX*dX-dX)*1000,NodeX),np.linspace(Yo*1000,(Yo+NodeY*dY-dY)*1000,NodeY), -np.linspace(Zo*1000,(Zo+NodeZ*dZ-dZ)*1000,NodeZ))
             T=[]
             pesi=[]
             for i in np.where(np.array(pick)!=0)[0]:
              T.append(pick[i]-np.min(np.array(pick)[np.where(np.array(pick)!=0)[0]])-P_teo_p[i])
              pesi.append(1/SNR[i]**2)
             To=np.min(np.array(pick)[np.where(np.array(pick)!=0)[0]])+(np.sum(np.array(T)*np.array(pesi))/np.sum(np.array(pesi)))
             for i in range(0,len(NAME)):
               P_teo.append(To+P_teo_p[i])
               S_teo.append((To+S_teo_p[i]))
             if second_loc==1:
              os.chdir(Path_grid)
              os.chdir(Path_OUT)
              for i in range(0,len(NAME)):
               ############################################################
               #SNR_P
               Trace_snr1=self.tr[i].copy()
               Trace_snr2=self.tr[i].copy()
               Pa_noise=(np.max(np.abs(Trace_snr1.trim(P_teo[i]-0.7,P_teo[i]-0.2).detrend('constant').data))*self.logger[i]/self.sensor[i])*100
               Pa=(np.max(np.abs(Trace_snr2.trim(P_teo[i],P_teo[i]+0.5).detrend('constant').data))*self.logger[i]/self.sensor[i])*100
               rowP=[NAME[i],"SNR_P:",Pa/Pa_noise]
               SNR_P.append(rowP)
               ##################
               #SNR_S componente X
               Trace_snr3=self.trX[i].copy()
               Trace_snr4=self.trX[i].copy()
               Sa_noiseX=(np.max(np.abs(Trace_snr3.trim(S_teo[i]-0.7,S_teo[i]-0.2).detrend('constant').data))*self.logger[i]/self.sensor[i])*100
               SaX=(np.max(np.abs(Trace_snr4.trim(S_teo[i],S_teo[i]+0.5).detrend('constant').data))*self.logger[i]/self.sensor[i])*100
               rowSX=[NAME[i],"SNR_SX:",SaX/Sa_noiseX]
               SNR_SX.append(rowSX)
               ##################
               #SNR_S componente Y
               Trace_snr5=self.trY[i].copy()
               Trace_snr6=self.trY[i].copy()
               Sa_noiseY=(np.max(np.abs(Trace_snr5.trim(S_teo[i]-0.7,S_teo[i]-0.2).detrend('constant').data))*self.logger[i]/self.sensor[i])*100
               SaY=(np.max(np.abs(Trace_snr6.trim(S_teo[i],S_teo[i]+0.5).detrend('constant').data))*self.logger[i]/self.sensor[i])*100
               rowSY=[NAME[i],"SNR_SY:",SaY/Sa_noiseY]
               SNR_SY.append(rowSY)
              os.chdir('..')
             lock_log.acquire()	
             os.chdir(Path_grid)	
             os.chdir(Path_OUT)							
             self.file_log=open('MParLoc_event.log', 'a')
             self.file_log.write('##################################################################################################################'+'\n')
             self.file_log.write(' ORIGIN TIME:  '+str(To)+'\n') 			 
             for i in np.where(np.array(pick)!=0)[0]:
               RES_P[i]=(P_teo[i]-pick[i])
               self.file_log.write(NAME[i]+'  Res_P:  '+str(RES_P[i])+' Error Pick: '+str(SNR[i])+'\n')

             for i in np.where(np.array(pickS)!=0)[0]:  
               RES_S[i]=((S_teo[i]-pickS[i]))
               self.file_log.write(NAME[i]+'  Res_S:  '+str(RES_S[i])+' Error Pick: '+str(SNR_S[i])+'\n')

             ######################
             self.file_log.write('##################################################################################################################'+'\n')
             for line in SNR_P:
              self.file_log.write(str(line[0])+' '+str(line[1])+' '+str(line[2])+'\n')
             for line in SNR_SX:
              self.file_log.write(str(line[0])+' '+str(line[1])+' '+str(line[2])+'\n')
             for line in SNR_SY:
              self.file_log.write(str(line[0])+' '+str(line[1])+' '+str(line[2])+'\n')
             ######################
             self.file_log.close()
             lock_log.release()
             os.chdir(Path_grid)	
             if len(Locx)>0:
               if len(xev)>0:
                ErrLoc=(np.sqrt((Locx-xev)**2+(Locy-yev)**2)/1000)	
                ErrLocx=(Locx-xev)/1000
                ErrLocy=(Locy-yev)/1000
                ErrLocz=(Locz-zev)/1000
                if len(ErrLocx)==1:
                 lock_log.acquire()
                 os.chdir(Path_grid)
                 os.chdir(Path_OUT)
                 self.file_log=open('MParLoc_event.log', 'a')
                 self.file_log.write('LOC_error:'+str(ErrLoc)+' '+'LOC_error_X: '+str(ErrLocx)+' '+'LOC_error_Y: '+str(ErrLocy)+'LOC_error_Z: '+str(ErrLocz)+'\n')  
                 self.file_log.close()
                 os.chdir(Path_grid)
                 lock_log.release()
               new_sum=False
               if second_loc==1:
                self.mappa.plot_location()
               if second_loc==0:
                self.Final_loc()
                second_loc=1
               if len(str(Locx))>5:
                change_factor=True
            else:
              os._exit(0)
           time.sleep(0.1)	  
	   

         def Final_loc(self):
           global data
           global AccZ
           global AccN
           global AccE		  
           global pick
           global xx
           global yy
           global Baz
           global new_sum
           global SNR
           global SNR_S
           SNR_S=(np.ones(shape=(len(NAME))))
           for staz_index in range(0,len(NAME)):
              if self.Pv[staz_index]==0 and pick[staz_index]!=0:
               accH=data[staz_index].copy()
               accL=data[staz_index].copy()
               acc=data[staz_index].copy()
               accE=dataX[staz_index].copy()
               accN=dataY[staz_index].copy()
               Trace_orig=acc.copy()
               Pa_noise=(np.max(np.abs(Trace_orig.trim(pick[staz_index]-(self.snr_wind_p+0.2),pick[staz_index]-0.2).detrend('constant').data))*self.logger[staz_index]/self.sensor[staz_index])*100
               Pa=(np.max(np.abs(acc.trim(pick[staz_index],pick[staz_index]+self.snr_wind_p).detrend('constant').data))*self.logger[staz_index]/self.sensor[staz_index])*100
               SNR[staz_index]=(Pa/Pa_noise)
               self.Pv[staz_index]=(np.max(np.abs(Trace_orig.slice(pick[staz_index]-3,pick[staz_index]+2).detrend('constant').filter("bandpass",freqmin=0.075,freqmax=25).data))*self.logger[staz_index]/self.sensor[staz_index])*100
               if SNR[staz_index] <= self.Tabella_P[0][0]:
                pick[staz_index]=self.Tabella_P[0][1]
               if SNR[staz_index] > self.Tabella_P[0][0]  and SNR[staz_index] <= self.Tabella_P[1][0]:
                SNR[staz_index]=self.Tabella_P[1][1]
               if SNR[staz_index] > self.Tabella_P[1][0] and SNR[staz_index] <= self.Tabella_P[2][0]:
                SNR[staz_index]=self.Tabella_P[2][1]
               if SNR[staz_index] > self.Tabella_P[2][0] and SNR[staz_index] <= self.Tabella_P[3][0]:
                SNR[staz_index]=self.Tabella_P[3][1]
               if SNR[staz_index] > self.Tabella_P[3][0]:
                SNR[staz_index]=self.Tabella_P[4][1]
              if  pickS[staz_index]!=0:
               accX=dataX[staz_index].copy()
               accY=dataY[staz_index].copy()
               Pa_noise=(np.max(np.abs(accX.slice(pickS[staz_index]-(self.snr_wind_S+0.2),pickS[staz_index]-0.2).detrend('constant').data))*self.logger[staz_index]/self.sensor[staz_index])*100
               Pa=(np.max(np.abs(accX.slice(pickS[staz_index],pickS[staz_index]+0.5).detrend('constant').data))*self.logger[staz_index]/self.sensor[staz_index])*100
               SNRX=(Pa/Pa_noise)
               Pa_noise=(np.max(np.abs(accY.slice(pickS[staz_index]-(self.snr_wind_S+0.2),pickS[staz_index]-0.2).detrend('constant').data))*self.logger[staz_index]/self.sensor[staz_index])*100
               Pa=(np.max(np.abs(accY.slice(pickS[staz_index],pickS[staz_index]+self.snr_wind_S).detrend('constant').data))*self.logger[staz_index]/self.sensor[staz_index])*100
               SNRY=(Pa/Pa_noise)
               SNR_S[staz_index]=np.max([SNRX,SNRY])
               if SNR_S[staz_index]<=self.Tabella_S[0][0]:
                pickS[staz_index]=self.Tabella_S[0][1]
               if SNR_S[staz_index] >self.Tabella_S[0][0]  and SNR_S[staz_index] <=self.Tabella_S[1][0]:
                SNR_S[staz_index]=self.Tabella_S[1][1]
               if SNR_S[staz_index] >self.Tabella_S[1][0] and SNR_S[staz_index] <=self.Tabella_S[2][0]:
                SNR_S[staz_index]=self.Tabella_S[2][1]
               if SNR_S[staz_index] >self.Tabella_S[2][0] and SNR_S[staz_index] <=self.Tabella_S[3][0]:
                SNR_S[staz_index]=self.Tabella_S[3][1]
               if SNR_S[staz_index] >self.Tabella_S[3][0]:
                SNR_S[staz_index]=self.Tabella_S[4][1]
           end_loc=False
           self.TsumErr=(np.ones(shape=(NodeX,NodeY,NodeZ)))
           self.TsumErr_S_P=(np.ones(shape=(NodeX,NodeY,NodeZ)))
           self.TsumErr_S=(np.ones(shape=(NodeX,NodeY,NodeZ)))
           if self.back_az:
            self.BAzsumErr=(np.ones(shape=(NodeX,NodeY,NodeZ)))
           if self.diff_Amp:
            self.RsumErr=(np.ones(shape=(NodeX,NodeY,NodeZ)))
           iter=9
           XMin=0
           XMax=NodeX
           YMin=0
           YMax=NodeY	
           ZMin=0
           ZMax=NodeZ
           Xvect_pre=[]
           Yvect_pre=[]		   
           Zvect_pre=[]	
           for iter_i in range(1,8,3):	
            iter_pre=copy.copy(iter)
            iter=8-iter_i
            if iter_i !=1:
              if np.sum((self.TsumErr[Xindexs,Yindexs,Zindexs]))!=0:
               TsumErr=(self.TsumErr.copy())
               TsumErr[Xindexs,Yindexs,Zindexs]=TsumErr[Xindexs,Yindexs,Zindexs]/np.sum((TsumErr[Xindexs,Yindexs,Zindexs]))
              else:
               TsumErr=(self.TsumErr.copy())	
              if np.sum((self.TsumErr_S_P[Xindexs,Yindexs,Zindexs]))!=0:
               TsumErr_S_P=(self.TsumErr_S_P.copy())
               TsumErr_S_P[Xindexs,Yindexs,Zindexs]=TsumErr_S_P[Xindexs,Yindexs,Zindexs]/np.sum((TsumErr_S_P[Xindexs,Yindexs,Zindexs]))
              else:
               TsumErr_S_P=1
              if np.sum((self.TsumErr_S[Xindexs,Yindexs,Zindexs]))!=0:
               TsumErr_S=(self.TsumErr_S.copy())
               TsumErr_S[Xindexs,Yindexs,Zindexs]=TsumErr_S[Xindexs,Yindexs,Zindexs]/np.sum((TsumErr_S[Xindexs,Yindexs,Zindexs]))
              else:
               TsumErr_S=1
              if self.back_az:
               if np.sum((self.BAzsumErr[Xindexs,Yindexs,Zindexs]))!=0:
                BAzsumErr=(self.BAzsumErr.copy())
                BAzsumErr[Xindexs,Yindexs,Zindexs]=BAzsumErr[Xindexs,Yindexs,Zindexs]/np.sum((BAzsumErr[Xindexs,Yindexs,Zindexs]))
               else:
                BAzsumErr=1
              else:
                BAzsumErr=1
              if self.diff_Amp:
               if np.sum((self.RsumErr[Xindexs,Yindexs,Zindexs]))!=0:
                RsumErr=(self.RsumErr.copy())
                RsumErr[Xindexs,Yindexs,Zindexs]=RsumErr[Xindexs,Yindexs,Zindexs]/np.sum((RsumErr[Xindexs,Yindexs,Zindexs]))
               else:
                RsumErr=1	
              else:
               RsumErr=1					
              if self.s_phases:
               self.Locmatrix=TsumErr*TsumErr_S_P*TsumErr_S*BAzsumErr*RsumErr
              else:
               self.Locmatrix=TsumErr*BAzsumErr*RsumErr
              self.Locmatrix=np.where(self.Locmatrix==1, 0, self.Locmatrix)
              if np.sum(self.Locmatrix)==0:
               Xindex, Yindex, Zindex=(np.where(self.Locmatrix==0)[0:3])
              else:
               Xindex, Yindex, Zindex=(np.where(self.Locmatrix>(np.max(self.Locmatrix)/100))[0:3])
              if np.min(Xindex)!=np.max(Xindex):
               XMin= np.min(Xindex)
               XMax= np.max(Xindex)
              else:
               if np.min(Xindex)-iter_pre>=0:
                XMin= np.min(Xindex)-iter_pre
               else:
                XMin=0
               if np.max(Xindex)+iter_pre<=NodeX:
                XMax= np.max(Xindex)+iter_pre
               else:
                XMax=NodeX
              if np.min(Yindex)!=np.max(Yindex):
               YMin= np.min(Yindex)
               YMax= np.max(Yindex)
              else:
               if np.min(Yindex)-iter_pre>=0:
                YMin= np.min(Yindex)-iter_pre
               else:
                YMin=0
               if np.max(Yindex)+iter_pre<=NodeY:
                YMax= np.max(Yindex)+iter_pre
               else:
                YMax=NodeY
              if np.min(Zindex)!=np.max(Zindex):
               ZMin= np.min(Zindex)
               ZMax= np.max(Zindex)
              else:
               if np.min(Zindex)-iter_pre>=0:
                ZMin= np.min(Zindex)-iter_pre
               else:
                ZMin=0
               if np.max(Zindex)+iter_pre<=NodeZ:
                ZMax= np.max(Zindex)+iter_pre
               else:
                ZMax=NodeZ
              Xvect_pre.extend(list(Xvect))
              Yvect_pre.extend(list(Yvect))
              Zvect_pre.extend(list(Zvect))
              Xvect=np.arange(XMin,XMax,iter)
              Yvect=np.arange(YMin,YMax,iter)
              Zvect=np.arange(ZMin,ZMax,iter) 
            else:
             Xvect=np.arange(XMin,XMax,iter)
             Yvect=np.arange(YMin,YMax,iter)
             Zvect=np.arange(ZMin,ZMax,iter)
             if np.max(Xvect)>NodeX:
              Xvect[-1]=NodeX
             if np.max(Yvect)>NodeY:
              Yvect[-1]=NodeY
             if np.max(Zvect)>NodeZ:
              Zvect[-1]=NodeZ	 
            Xindexs,Yindexs,Zindexs=np.meshgrid(Xvect,Yvect,Zvect)
            self.arrival_matrix=None
            self.arrival_matrix_S=None
            self.Rnode=None
            self.arrival_matrix=[]
            self.arrival_matrix_S=[]
            self.Rnode=[]
            Xindexs_i,Yindexs_i,Zindexs_i=np.where(self.TsumErr[Xindexs,Yindexs,Zindexs]==1)[0:3]
            Xindexs=Xindexs[Xindexs_i,Yindexs_i,Zindexs_i]
            Yindexs=Yindexs[Xindexs_i,Yindexs_i,Zindexs_i]
            Zindexs=Zindexs[Xindexs_i,Yindexs_i,Zindexs_i]
	
            os.chdir(Path_grid)
            for indxc in range(0,len(NAME)):
             if pick[indxc]!=0:
              f = open(Model_path_high+r'/arrival_matrix'+NAME[indxc]+'.bin', "rb")
              floats = array(self.grid_precision)
              floats.fromfile(f,NodeX*NodeY*NodeZ)
              floats = np.reshape(floats,(NodeX,NodeY,NodeZ))[Xindexs,Yindexs,Zindexs]
              self.arrival_matrix.append(floats)
             else:
              self.arrival_matrix.append([])	 
             if pickS[indxc]!=0:
              f = open(Model_path_high+r'/arrival_matrixS'+NAME[indxc]+'.bin', "rb")
              floats = array(self.grid_precision)
              floats.fromfile(f,NodeX*NodeY*NodeZ)
              floats = np.reshape(floats,(NodeX,NodeY,NodeZ))[Xindexs,Yindexs,Zindexs]
              self.arrival_matrix_S.append(floats)             
             else:
              self.arrival_matrix_S.append([])
             if self.Pv[indxc]!=0:
              f = open(Model_path_high+r'/Rnode'+NAME[indxc]+'.bin', "rb")
              floats = array(self.grid_precision)
              floats.fromfile(f,NodeX*NodeY*NodeZ)
              floats = np.reshape(floats,(NodeX,NodeY,NodeZ))[Xindexs,Yindexs,Zindexs]
              self.Rnode.append(floats)
             else:
              self.Rnode.append([])
            os.chdir(Path_OUT)			

            for staz_index in range(0,len(NAME)):
              if pick[staz_index]!=0 or pickS[staz_index]!=0:
                if Baz[staz_index][0]==-9999 and pick[staz_index]!=0 and self.back_az:
                 #########################BAZ ESTIMATION#######################################
                 acc=data[staz_index].copy()
                 accE=dataX[staz_index].copy()
                 accN=dataY[staz_index].copy()
                 Baz_m, BAZ_error =BAZ(acc,accE,accN,staz_index)
                 Baz[staz_index][0]=Baz_m
                 if (Baz_m)!=0:
                  os.chdir(Path_grid)
                  f = open(Model_path_high+r'/BAznode'+NAME[indxc]+'.bin', "rb")
                  floats = array(self.grid_precision)
                  floats.fromfile(f,NodeX*NodeY*NodeZ)
                  floats = np.reshape(floats,(NodeX,NodeY,NodeZ))[Xindexs,Yindexs,Zindexs]
                  self.BAznode=floats
                  os.chdir(Path_OUT)
                  Residui=np.minimum((((((self.BAznode-Baz[staz_index][0])))))*2,(((((self.BAznode-(Baz[staz_index][0]+360))))))*2)
                  Residui=np.minimum(Residui,((((((self.BAznode+360)-(Baz[staz_index][0]))))))**2)
                  self.BAzsumErr[Xindexs,Yindexs,Zindexs]=(self.BAzsumErr[Xindexs,Yindexs,Zindexs]*np.exp(-Residui/(2*((self.back_az_error))**2))) 
                if pickS[staz_index]!=0 and pick[staz_index]!=0 :
                  self.TsumErr_S_P[Xindexs,Yindexs,Zindexs]=(self.TsumErr_S_P[Xindexs,Yindexs,Zindexs]*(1/((SNR[staz_index]+SNR_S[staz_index])))*np.exp(-((((self.arrival_matrix_S[staz_index]-self.arrival_matrix[staz_index])-((pickS[staz_index]-pick[staz_index]))))**2/(((SNR[staz_index]+SNR_S[staz_index])*2)**2))))
                for i in range(0,len(NAME)):
                  if self.Pv[i]!=0 and self.Pv[staz_index]!=0 and self.diff_Amp:
                    self.RsumErr[Xindexs,Yindexs,Zindexs]=self.RsumErr[Xindexs,Yindexs,Zindexs]*np.exp(-(((self.Rnode[staz_index]-self.Rnode[i])-(+self.diff_Amp_b*math.log10(self.Pv[staz_index])-self.diff_Amp_b*math.log10(self.Pv[i])))**2/(2*(self.diff_Amp_error)**2)))
                  if pickS[staz_index]!=0 and pickS[i]!=0 and i!=staz_index and pick[staz_index]==0 :             
                    self.TsumErr_S[Xindexs,Yindexs,Zindexs]=(self.TsumErr_S[Xindexs,Yindexs,Zindexs]*(1/((SNR_S[staz_index]+SNR_S[i])))*np.exp(-(((self.arrival_matrix_S[staz_index]-self.arrival_matrix_S[i])-(pickS[staz_index]-pickS[i]))**2/((((SNR_S[staz_index]+SNR_S[i]))*2)**2))))
                  if pick[staz_index]!=0 and pick[i]!=0 and i!=staz_index :
                    self.TsumErr[Xindexs,Yindexs,Zindexs]=(self.TsumErr[Xindexs,Yindexs,Zindexs]*(1/((SNR[staz_index]+SNR[i])))*np.exp(-(((self.arrival_matrix[staz_index]-self.arrival_matrix[i])-(pick[staz_index]-pick[i]))**2/(((SNR[staz_index]+SNR[i])*2)**2))))
           new_sum=True
           end_loc=True
		   
		   
def Launcher(Network_folder,Event_folder):
 global win
 global Network_path
 global Event_path
 Event_path=Event_folder
 Network_path=Network_folder
 win=MainCode() 