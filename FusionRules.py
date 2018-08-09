import numpy as np
from osgeo import gdal
import sys
import pywt
import scipy as sc
from scipy import stats
from sklearn.preprocessing import normalize
from numpy import linalg as la
import math
def readFile(filename):
    filehandle = gdal.Open(filename)
    band1 = filehandle.GetRasterBand(1)
    geotransform = filehandle.GetGeoTransform()
    geoproj = filehandle.GetProjection()
    band1data=band1.ReadAsArray()
    xsize=filehandle.RasterXSize
    ysize=filehandle.RasterYSize
    met = filehandle.GetMetadata()
    return xsize,ysize,geotransform, geoproj, band1data,met

def writeFile(filename,geotransform,geoprojection,metadata,data):
    (x,y) = data.shape
    format = "GTiff"
    driver = gdal.GetDriverByName(format)
    #dst_datatype = gdal.GDT_Byte
    dst_datatype=gdal.GDT_UInt16
    #print (dst_datatype)
    dst_ds = driver.Create(filename,y,x,1,dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(geoprojection)
    dst_ds.SetMetadata(metadata)
    return 1

def hist_match(s,t):
    return ((np.std(t)/np.std(s))*(s-np.mean(s)))+ np.mean(t)       #returns histogram matched matrix


# Function to calculate local energy s1,s2 and matching function M
# Returns fused coefficients of detail band
def localEnergy(cX_MS,cX_SAR):
    MS=np.abs(cX_MS)
    SAR = np.abs(cX_SAR)
    [a,b]=np.shape(cX_MS)
    s1=np.zeros((a,b))                                                                  #local energy of each pixel of MS image
    s2=np.zeros((a,b))                                                                  #local energy of each pixel of SAR image    
    c=np.zeros((a,b))
    M=np.zeros((a,b))                                                                   #Matching function 
    cX_MS_new=np.zeros((a,b))
    cX_MS_pad=np.pad(MS, 1,'constant')
    cX_SAR_pad=np.pad(SAR, 1,'constant')
    for i in range(1,a+1):
        for j in range(1,b+1):
            s1[i-1,j-1]=0
            s2[i-1,j-1]=0
            c[i-1,j-1]=0
            M[i-1,j-1]=0
            for m in range(-1,1):
                for n in range(-1,1):
                    s1[i-1,j-1]=s1[i-1,j-1]+cX_MS_pad[i+m,j+n]**2
                    s2[i-1,j-1]=s2[i-1,j-1]+cX_SAR_pad[i+m,j+n]**2
                    c[i-1,j-1] = c[i-1,j-1]+(cX_MS_pad[i+m,j+n]*cX_SAR_pad[i+m,j+n])
            M[i-1,j-1]=(2*c[i-1,j-1])/(s1[i-1,j-1]+s2[i-1,j-1])
            # Fused coeff based on Matching function
            if M[i-1,j-1]<=0.5:
                if s1[i-1,j-1]<s2[i-1,j-1]:
                    cX_MS_new[i-1,j-1]= cX_SAR[i-1,j-1]
                else:
                    cX_MS_new[i-1,j-1]= cX_MS[i-1,j-1]
            else:
                wmin = (M[i-1,j-1]-0.5)/(2*(1-0.5))
                wmax = 1-wmin
                if s1[i-1,j-1]>= s2[i-1,j-1]:
                     cX_MS_new[i-1,j-1]=wmax*cX_MS[i-1,j-1] + wmin*cX_SAR[i-1,j-1]
                else:
                     cX_MS_new[i-1,j-1]=wmin*cX_MS[i-1,j-1] + wmax*cX_SAR[i-1,j-1]
    return cX_MS_new

# Function to calculate fused coeff of approximate band
# Returns fused coeff of approximate band
def cA(cA_MS,cA_SAR):
    [a,b]=np.shape(cA_MS)
    s1=np.zeros((a,b))
    s2=np.zeros((a,b))
    c=np.zeros((a,b))
    M=np.zeros((a,b))
    cA_MS_new=np.zeros((a,b))
    cA_MS_pad=np.pad(cA_MS, 1,'constant')
    cA_SAR_pad=np.pad(cA_SAR, 1,'constant')
    for i in range(1,a+1):
        for j in range(1,b+1):
            s1[i-1,j-1]=0
            s2[i-1,j-1]=0
            c[i-1,j-1]=0
            M[i-1,j-1]=0
            for m in range(-1,1):
                for n in range(-1,1):
                    s1[i-1,j-1]=s1[i-1,j-1]+cA_MS_pad[i+m,j+n]**2
                    s2[i-1,j-1]=s2[i-1,j-1]+cA_SAR_pad[i+m,j+n]**2
                    c[i-1,j-1] = c[i-1,j-1]+(cA_MS_pad[i+m,j+n]*cA_SAR_pad[i+m,j+n])
            M[i-1,j-1]=(2*c[i-1,j-1])/(s1[i-1,j-1]+s2[i-1,j-1])
            if M[i-1,j-1]>0.7:
                cA_MS_new[i-1,j-1] = 0.5*(cA_MS[i-1,j-1]+cA_SAR[i-1,j-1])
            else:
                if s1[i-1,j-1]>= s2[i-1,j-1]:
                    cA_MS_new[i-1,j-1]=cA_MS[i-1,j-1]
                else:
                    cA_MS_new[i-1,j-1]=cA_SAR[i-1,j-1]
    return cA_MS_new

#Fusion rule based on Maxima
def maxRule(m,s):
    a1=np.abs(s) #---activity measure of band1 (SAR MSD Band)
    a2=np.abs(m) #--- activity measure of band2 (MS MSD Band)
    temp=np.where(a1>a2)
    m[temp]=s[temp]
    return m

def CCC(o,f):
    o=np.float64(o)
    f=np.float64(f)
    t1=o-np.mean(o)
    t2=f-np.mean(f)
    temp2=np.sum(np.multiply(t1,t2))
    t3=np.sum(np.square(t1))
    t4=np.sum(np.square(t2))
    temp3=np.sqrt(t3*t4)
    cc=temp2/temp3
    return cc

def bias(orig,fused):
##    orig = np.absolute(orig)
##    fused = np.absolute(fused)
    return ((orig.mean()-fused.mean())/orig.mean())

def relVar(orig,fused):
    return (orig.var()-fused.var())/orig.var()

def entropy(m):
    m = m.flatten()
    p_data = np.bincount(m)/len(m)
    entropy = sc.stats.entropy(p_data)
    return entropy
 
