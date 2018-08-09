import numpy as np
from osgeo import gdal
import sys
import pywt
from functionFile import readFile
from functionFile import writeFile
from functionFile import localEnergy
from functionFile import cA
from functionFile import CC
from functionFile import hist_match
from functionFile import maxRule
from msd_fuse_rule_SCK import fuse_band

# Reading data
MS = 'F:\\SY_MTech\\Manali\\imagery_HV\\Registration\\Reprojection_registered\\New folder\\band5_gamma.tif'
#SAR = 'F:\\SY_MTech\\Manali\\imagery_HV\\HV_clip9Nov\\HV_clip9Nov_gamma3_1.tif'
SAR = 'F:\\SY_MTech\\Manali\\imagery_HV\\HV_clip9Nov\\HV_clip9Nov_e_lee8_1.tif'
(x_MS,y_MS,transform_MS,projection_MS,data_MS,met_MS) = readFile(MS)                            #-----Reading MS data
(x_SAR,y_SAR,transform_SAR,projection_SAR,data_SAR,met_SAR) = readFile(SAR)                     #-----Reading SAR data
# Convert data to float64
data_MS = np.float64(data_MS)
data_SAR = np.float64(data_SAR)

print (data_MS.shape)
# wavelet decomposition
coeffs = pywt.wavedec2(data_MS, 'db16',level=3)
[cA,(cH3,cV3,cD3),(cH2,cV2,cD2),(cH,cV,cD)]=coeffs                                              #-----detail and approx coeff of MS data
data_SAR=hist_match(data_SAR,data_MS)                                                           #-----histogram matched SAR data                      
coeffs = pywt.wavedec2(data_SAR, 'db16',level=3)
[cA_SAR,(cH3_SAR,cV3_SAR,cD3_SAR),(cH2_SAR,cV2_SAR,cD2_SAR),(cH_SAR,cV_SAR,cD_SAR)]=coeffs      #-----detail and approx coeff of SAR data

##--------------------------------------------------------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------------------------------------------------------
##Fusion using Rule1 - Maximum Rule
cH3n= maxRule(cH3,cH3_SAR)                                                                      
cV3n= maxRule(cV3,cV3_SAR)
cD3n= maxRule(cD3,cD3_SAR)
cH2n = maxRule(cH2,cH2_SAR)
cV2n= maxRule(cV2,cV2_SAR)
cD2n= maxRule(cD2,cD2_SAR)
cHn= maxRule(cH,cH_SAR)
cVn= maxRule(cV,cV_SAR)
cDn= maxRule(cD,cD_SAR)
cAn=((cA)+(cA_SAR))/2

coeffs = [cAn,(cH3n,cV3n,cD3n),(cH2n,cV2n,cD2n),(cHn,cVn,cDn)]                                  #-----fused coeffs
data = pywt.waverec2(coeffs,'db16',mode='symmetric')                                            #-----wavelet reconstruction of MS image
file = 'F:\\SY_MTech\\Manali\\imagery_HV\\Fusion\\Reprojected\\SAR_Image_Elee_db16\\R1_cA0.5\\Band5_gamma_L3_R1_cA0.5.tif'
writeFile(file,transform_MS,projection_MS,met_MS,data)                                          #-----write Fused image
print('Rule 1 Complete')

##--------------------------------------------------------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------------------------------------------------------
#Fusion using Rule2 - Simple Averaging
cH3= (cH3+cH3_SAR)/2
cV3= (cV3+cV3_SAR)/2
cD3= (cD3+cD3_SAR)/2
cH2 = (cH2+cH2_SAR)/2
cV2= (cV2+cV2_SAR)/2
cD2= (cD2+cD2_SAR)/2
cH= (cH+cH_SAR)/2
cV= (cV+cV_SAR)/2
cD= (cD+cD_SAR)/2
cAn=((cA)+(cA_SAR))/2

coeffs = [cAn,(cH3,cV3,cD3),(cH2,cV2,cD2),(cH,cV,cD)]                                            #-----fused coeffs
data = pywt.waverec2(coeffs,'db16',mode='symmetric')                                             #-----wavelet reconstruction of MS image
file = 'F:\\SY_MTech\\Manali\\imagery_HV\\Fusion\\Reprojected\\SAR_Image_Elee_db16\\R2_cA0.5\\Band5_gamma_L3_R2_cA0.5.tif'
writeFile(file,transform_MS,projection_MS,met_MS,data)                                          #-----write Fused image
print('Rule 2 Complete')

##--------------------------------------------------------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------------------------------------------------------
#Fusion using Rule3 - Fusion based on Local Energy of 3*3 neighbourhood
cH_new= localEnergy(cH,cH_SAR)        
cV_new= localEnergy(cV,cV_SAR)
cD_new= localEnergy(cD,cD_SAR)
cH2_new= localEnergy(cH2,cH2_SAR)        
cV2_new= localEnergy(cV2,cV2_SAR)
cD2_new= localEnergy(cD2,cD2_SAR)
cH3_new= localEnergy(cH3,cH3_SAR)        
cV3_new= localEnergy(cV3,cV3_SAR)
cD3_new= localEnergy(cD3,cD3_SAR)
cA_MS = (cA+cA_SAR)/2

# MultiLevel Wavelet Reconstruction
coeffst_FS = [cA_MS,(cH3_new,cV3_new,cD3_new),(cH2_new,cV2_new,cD2_new),(cH_new,cV_new,cD_new)] #-----fused coeffs
rec_FS = pywt.waverec2(coeffst_FS,'db16')                                                    #-----wavelet reconstruction of MS image

#Write Image
filename = 'F:\\SY_MTech\\Manali\\imagery_HV\\Fusion\\Reprojected\\SAR_Image_Elee_db16\\R3\\Band5_gamma_L3_R3_cA0.5.tif'
writeFile(filename,transform_MS,projection_MS,met_MS,rec_FS)                                #-----write Fused image
print('Rule 3 Complete')

##--------------------------------------------------------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------------------------------------------------------
#Fusion using Rule4 - Activity Measure based on maxima over 3*3 nhood
cH_new= fuse_band(cH_SAR,cH,1)        
cV_new= fuse_band(cV_SAR,cV,1)
cD_new= fuse_band(cD_SAR,cD,1)
cH2_new= fuse_band(cH2_SAR,cH2,1)        
cV2_new= fuse_band(cV2_SAR,cV2,1)
cD2_new= fuse_band(cD2_SAR,cD2,1)
cH3_new= fuse_band(cH3_SAR,cH3,1)        
cV3_new= fuse_band(cV3_SAR,cV3,1)
cD3_new= fuse_band(cD3_SAR,cD3,1)
cA_MS = (cA+cA_SAR)/2

# MultiLevel Wavelet Reconstruction
coeffst_FS = [cA_MS,(cH3_new,cV3_new,cD3_new),(cH2_new,cV2_new,cD2_new),(cH_new,cV_new,cD_new)] #-----fused coeffs
rec_FS = pywt.waverec2(coeffst_FS,'db16')                                                    #-----wavelet reconstruction of MS image

#Write Image
filename = 'F:\\SY_MTech\\Manali\\imagery_HV\\Fusion\\Reprojected\\SAR_Image_Elee_db16\\R4\\Band5_gamma_L3_R4_cA0.5.tif'
writeFile(filename,transform_MS,projection_MS,met_MS,rec_FS)                                #-----write Fused image
print('Rule 4 Complete')
