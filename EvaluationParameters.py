import numpy as np
from osgeo import gdal
import sys
import pywt
from functionFile import CC
from FusionRules import CCC
from FusionRules import readFile
from fus_eval_ref_SCK import SAM
#from functionFile import SAM
from FusionRules import entropy
from FusionRules import bias
from FusionRules import relVar
from fus_eval_ref_SCK import rel_sdd
MS2 = 'F:\\SY_MTech\\Manali\\imagery_HV\\Reprojected\\band2\\AW-NF45S-105-055-23nov13-BAND2.tif'
(x_MS2,y_MS2,transform_MS2,projection_MS2,data_MS2,met_MS2) = readFile(MS2)                            #-----Reading original band2 data
MS3 = 'F:\\SY_MTech\\Manali\\imagery_HV\\Reprojected\\band3\\AW-NF45S-105-055-23nov13-BAND3.tif'
(x_MS3,y_MS3,transform_MS3,projection_MS3,data_MS3,met_MS3) = readFile(MS3)                            #-----Reading original band3 data
MS4 = 'F:\\SY_MTech\\Manali\\imagery_HV\\Reprojected\\band4\\AW-NF45S-105-055-23nov13-BAND4.tif'
(x_MS4,y_MS4,transform_MS4,projection_MS4,data_MS4,met_MS4) = readFile(MS4)                            #-----Reading original band4 data
MS5 = 'F:\\SY_MTech\\Manali\\imagery_HV\\Reprojected\\band5\\AW-NF45S-105-055-23nov13-BAND5.tif'
(x_MS5,y_MS5,transform_MS5,projection_MS5,data_MS5,met_MS5) = readFile(MS5)                            #-----Reading original band5 data

##data_MS2 = np.float64(data_MS2)
##data_MS3 = np.float64 (data_MS3)
##data_MS4 = np.float64(data_MS4)
##data_MS5 = np.float64 (data_MS5)

f2 ='F:\\SY_MTech\\Manali\\imagery_HV\\Fusion\\Reprojected\\SAR_Image_Elee_db16\\R4\\RP_Band2_gamma_L3_R4_cA0.5.tif'
(x_f2,y_f2,transform_f2,projection_f2,data_f2,met_f2) = readFile(f2)                            #-----Reading fused MS data
f3 ='F:\\SY_MTech\\Manali\\imagery_HV\\Fusion\\Reprojected\\SAR_Image_Elee_db16\\R4\\RP_Band3_gamma_L3_R4_cA0.5.tif'
(x_f3,y_f3,transform_f3,projection_f3,data_f3,met_f3) = readFile(f3)                            #-----Reading fused MS data
f4 ='F:\\SY_MTech\\Manali\\imagery_HV\\Fusion\\Reprojected\\SAR_Image_Elee_db16\\R4\\RP_Band4_gamma_L3_R4_cA0.5.tif'
(x_f4,y_f4,transform_f4,projection_f4,data_f4,met_f4) = readFile(f4)                            #-----Reading fused MS data
f5 ='F:\\SY_MTech\\Manali\\imagery_HV\\Fusion\\Reprojected\\SAR_Image_Elee_db16\\R4\\RP_Band5_gamma_L3_R4_cA0.5.tif'
(x_f5,y_f5,transform_f5,projection_f5,data_f5,met_f5) = readFile(f5)                            #-----Reading fused MS data

##data_f2 = np.float64(data_f2)
##data_f3 = np.float64 (data_f3)
##data_f4 = np.float64(data_f4)
##data_f5 = np.float64 (data_f5)

## Bias
print('Rel Variance')
b2 = relVar(data_MS2,data_f2)
print(b2)
b3 = relVar(data_MS3,data_f3)
print(b3)
b4 = relVar(data_MS4,data_f4)
print(b4)
b5 = relVar(data_MS5,data_f5)
print(b5)


print('Correlation Coeficients')
cc = CCC(data_MS2,data_f2[0:2332,0:2199])
print(cc)
cc = CCC(data_MS3,data_f3[0:2332,0:2199])
print(cc)
cc = CCC(data_MS4,data_f4[0:2332,0:2199])
print(cc)
cc = CCC(data_MS5,data_f5[0:2332,0:2199])
print(cc)


##sam = SAM(data_MS2,data_MS3,data_MS4,data_MS5,data_f2,data_f3,data_f4,data_f5)
##print('SAM',sam)


print('entropy')
e = entropy(data_f2)
print (e)
e = entropy(data_f3)
print (e)

e = entropy(data_f4)
print (e)

e = entropy(data_f5)
print (e)


print('Std Deviation Calculation')
s = rel_sdd(data_MS2,data_f2[0:2332,0:2199])
print(s)
s = rel_sdd(data_MS3,data_f3[0:2332,0:2199])
print(s)
s = rel_sdd(data_MS4,data_f4[0:2332,0:2199])
print(s)
s = rel_sdd(data_MS5,data_f5[0:2332,0:2199])
print(s)


