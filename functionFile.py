import numpy as np
from osgeo import gdal
import sys
import pywt
import scipy as sc
import pandas as pd
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
    print (dst_datatype)
    dst_ds = driver.Create(filename,y,x,1,dst_datatype)
    dst_ds.GetRasterBand(1).WriteArray(data)
    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(geoprojection)
    dst_ds.SetMetadata(metadata)
    return 1

def hist_match(s,t):
    return ((np.std(t)/np.std(s))*(s-np.mean(s)))+ np.mean(t)


def localEnergy(cX_MS,cX_SAR):
    MS=np.abs(cX_MS)
    SAR = np.abs(cX_SAR)
    [a,b]=np.shape(cX_MS)
    s1=np.zeros((a,b))
    s2=np.zeros((a,b))
    c=np.zeros((a,b))
    M=np.zeros((a,b))
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
            if M[i-1,j-1]<=0.5:
                #count=count+1
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
##    print(M.max(),M.min(),"\n")
    return cX_MS_new


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

def maxRule(m,s):
    a1=np.abs(s) #---activity measure of band1 (SAR MSD Band)
    a2=np.abs(m) #--- activity measure of band2 (MS MSD Band)
    temp=np.where(a1>a2)
    m[temp]=s[temp]
    return m
 

def bias(orig,fused):
##    orig = np.absolute(orig)
##    fused = np.absolute(fused)
    return ((orig.mean()-fused.mean())/orig.mean())

def relVar(orig,fused):
    return (orig.var()-fused.var())/orig.var()

def CC(o,f):
    [a,b] = np.shape(o)
    x = np.mean(o)
    y = np.mean(f)
    nm = 0
    dn1 = 0
    dn2 = 0
    for i in range(0,a-1):
        for j in range(0,b-1):
            nm = nm + ((o[i,j]-x)*(f[i,j]-y))
            dn1 =  dn1 + ((o[i,j]-x)**2)
            dn2 = dn2 + ((f[i,j]-y)**2)
    cc = nm/((dn1*dn2)**(0.5))
    return cc

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

def SAM(o2,o3,o4,o5,f2,f3,f4,f5):
    [a,b]=o2.shape
    o2 = o2/np.linalg.norm(o2)
    o3 = o3/np.linalg.norm(o3)
    o4 = o4/np.linalg.norm(o4)
    o5 = o5/np.linalg.norm(o5)
    f2 = f2/np.linalg.norm(f2)
    f3 = f3/np.linalg.norm(f3)
    f4 = f4/np.linalg.norm(f4)
    f5 = f5/np.linalg.norm(f5)
    sam = np.zeros([a,b])
    for i in range(0,a):
        for j in range(0,b):
            s = ((o2[i,j]*f2[i,j])+(o3[i,j]*f3[i,j])+(o4[i,j]*f4[i,j])+(o5[i,j]*f5[i,j])/(np.sqrt((o2[i,j]**2)+(o3[i,j]**2)+(o4[i,j]**2)+(o5[i,j]**2))*np.sqrt((f2[i,j]**2)+(f3[i,j]**2)+(f4[i,j]**2)+(f5[i,j]**2))))
            if s>=-1 and s<=1:
                sam[i,j] = np.arccos(s)
            else:
                s ==np.nan_to_num(s)
                sam[i,j] = np.arccos(s)
    sam_g = np.mean(sam)
    return sam_g

def entropy(m):
    m = m.flatten()
    p_data = np.bincount(m)/len(m)
    entropy = sc.stats.entropy(p_data)
    return entropy
    


    
