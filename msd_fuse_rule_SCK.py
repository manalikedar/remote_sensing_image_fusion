import numpy as np
from numpy import linalg as la
from numpy import matlib
import scipy.ndimage as sc
from scipy import signal
def fuse_band(band1,band2,k):
#activity calculation is coefficient based and activity measure is abs(band)#
#----fusion rule is maximum of activity ----------#
    if k==0:
       a1=np.abs(band1) #---activity measure of band1 (SAR MSD Band)
       a2=np.abs(band2) #--- activity measure of band2 (MS MSD Band)
       temp=np.where(a1<a2)
       band1[temp]=band2[temp]
    return band1

    if k==1:
    #activity calculation is window based and activity measure is abs(band) over a window#
    #----fusion rule is maximum of activity ----------#
      a1=np.abs(band1) #---activity measure of band1(SAR)
      a2=np.abs(band2) #--- activity measure of band2(MS)
      nhood=3
      a11=sc.maximum_filter(a1,nhood,mode='constant') #----max filter of scipy.ndimage
      a22=sc.maximum_filter(a2,nhood,mode='constant')
      temp=np.where(a11<a22)
      #----consistency verification over 5*5 neighbourhood----#
      bd_map=np.zeros(band1.shape)
      bd_map[temp]=1 #----decision map value 0 shows that MSD coeff. is from band 1 and 1 shows that MSD coeff. is from band 2
      bd_map1=sc.median_filter(bd_map,5,mode='constant')
      fused_band=np.zeros(band1.shape)
      temp1=np.where(bd_map1==0)
      fused_band[temp1]=band1[temp1]
      temp2=np.where(bd_map1==1)
      fused_band[temp2]=band2[temp2]
    return fused_band

    if k==2:
        nhood=5
        var_band1=sc.filters.generic_filter(band1,np.var,size=(5,5)) #------local variance of band1(SAR band)-----#
        var_band2=sc.filters.generic_filter(band2,np.var,size=(5,5)) #------local variance of band2 (MS band)-----#
        temp=np.where(var_band1<var_band2)
        #----consistency verification over 5*5 neighbourhood----#
        bd_map=np.zeros(band1.shape)
        bd_map[temp]=1 #----decision map value 0 shows that MSD coeff. is from band 1 and 1 shows that MSD coeff. is from band 2
        bd_map1=sc.median_filter(bd_map,5,mode='constant')
        fused_band=np.zeros(band1.shape)
        temp1=np.where(bd_map1==0)
        fused_band[temp1]=band1[temp1]
        temp2=np.where(bd_map1==1)
        fused_band[temp2]=band2[temp2]
    return fused_band
        
    


