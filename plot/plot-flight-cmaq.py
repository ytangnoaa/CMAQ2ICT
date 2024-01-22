#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
import sys
from datetime import datetime
import cartopy.crs as ccrs
import matplotlib.colors as colors
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

PLOTMAP=1	
a=pd.read_csv('comp-dc8-1.index')
nvar,x=a.shape
if x != 3:
 print('wrong index file',x)
 sys.exit(1)
 
if len(sys.argv) != 4:
  print('need flight start_index end_index output_figure')
  sys.exit(1)


sflight=sys.argv[1]
eflight=sys.argv[2]

pdf=bpdf.PdfPages(sys.argv[3])
np=1
p={}
k={}

for x in range(int(sflight),int(eflight)+1):
 
 b=pd.read_csv('../data/dc8-1m-firevoc-'+str(x).zfill(2)+'.dat')
 c=pd.read_csv('../data/dc8-1m-wrfcmaq-'+str(x).zfill(2)+'.dat')

 cdate=datetime.strptime('2019'+str(b['JDAY'][0]),'%Y%j')
 dstring=cdate.strftime('%m/%d/%Y')
 print('now process flight ',x,'  ',dstring)

# plot Flight map
 if PLOTMAP == 1:
   ax=plt.axes(projection=ccrs.PlateCarree())
   plt.axis([-130,-90,25,50])
   cs3=ax.scatter(b['Longitude (Degs)'],b['Latitude (Degs)'],c=b['Alt_GPS (m)'],s=0.5,marker='D',
    norm=colors.Normalize(vmax=b['Alt_GPS (m)'].max()))
   cbar=plt.colorbar(cs3,fraction=0.07)
   cbar.set_label('Altitude (m)',size=18)
   cbar.ax.tick_params(direction='in',labelsize=10)
 
   ax.set_title('FireX-AQ DC-8 Flight on '+dstring,fontsize=16)
   ax.add_feature(cfeature.BORDERS)
   ax.add_feature(cfeature.COASTLINE)
   ax.add_feature(states_provinces, edgecolor='gray')
   gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
   gl.xlabels_top = False
   gl.ylabels_right = False
   gl.xformatter = LONGITUDE_FORMATTER
   gl.yformatter = LATITUDE_FORMATTER
   gl.ylabel_style = {'size':10, 'weight':'bold'}
   gl.xlabel_style = {'size':10, 'weight':'bold'}
   pdf.savefig(plt.gcf().number)
   plt.clf()

# plot species
 for i in range(nvar):
  if i%6 == 0:
   fpage,ax=plt.subplots(3,2,figsize=(8.5,11))
   plt.subplots_adjust(wspace=0.6,hspace=0.4)
   p[i],p[i+1],p[i+2],p[i+3],p[i+4],p[i+5]=ax[0,0],ax[0,1],ax[1,0],ax[1,1],ax[2,0],ax[2,1]
  
  print('plot ',a['obs_name'][i],a['model_name'][i])
  
  k[0],=p[i].plot(b['UTC (hour)'],pd.to_numeric(b[a['obs_name'][i]],errors='coerce'),'-o',c='k',label='Obs')
  if 'model_J-' in a['model_name'][i]:
   b[a['model_name'][i]]=b[a['model_name'][i]]/60 # /min -> /s
   c[a['model_name'][i]]=c[a['model_name'][i]]/60
  elif 'model_AOE' in a['model_name'][i]:
   b[a['model_name'][i]]=b[a['model_name'][i]]*1000 # /km -> /Mm
   c[a['model_name'][i]]=c[a['model_name'][i]]*1000
  k[1],=p[i].plot(b['UTC (hour)'],b[a['model_name'][i]],'r-',label='GFS-CMAQ')
  k[2],=p[i].plot(c['UTC (hour)'],c[a['model_name'][i]],'b',ls='-.',label='WRF-CMAQ')

  par1=p[i].twinx()
  k[3],=par1.plot(b['UTC (hour)'],b['Alt_GPS (m)'],'y--',label='Altitude')
  
  p[i].set_xlabel('UTC (hour)',fontsize=14)
  p[i].set_ylabel(a['varname'][i],fontsize=14)
  par1.set_ylabel('Altitude (m)',fontsize=14)
  p[i].set_title('FireX-AQ DC-8 Flight '+dstring,fontsize=14)
  
  lines=[k[0],k[1],k[2],k[3]]
  p[i].legend(lines,[l.get_label() for l in lines],fontsize=12)
  
  if i%6 ==5:
   pdf.savefig(plt.gcf().number)
   np=np+1
   plt.clf()
 
 if i%6 != 5:  # end of one flight
   pdf.savefig(plt.gcf().number)
   np=1
   plt.clf()
   
pdf.close()
