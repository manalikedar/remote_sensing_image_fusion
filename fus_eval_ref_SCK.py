import numpy as np
from numpy import linalg as la
from numpy import matlib
import math
np.seterr(divide='ignore', invalid='ignore')
def cor_coef(i1,i2):
    i1=np.float64(i1)
    i2=np.float64(i2)
    t1=i1-np.mean(i1)
    t2=i2-np.mean(i2)
    temp2=np.sum(np.multiply(t1,t2))
    t3=np.sum(np.square(t1))
    t4=np.sum(np.square(t2))
    temp3=np.sqrt(t3*t4)
    cc=temp2/temp3
    return cc

def ERGAS(i2_ref,i3_ref,i4_ref,i5_ref,i2_f,i3_f,i4_f,i5_f):
    
    [m,n]=i2_ref.shape
    m=np.float64(m)
    n=np.float64(n)
    i_ref=np.concatenate(([i2_ref],[i3_ref],[i4_ref],[i5_ref]))
    i_fuse=np.concatenate(([i2_f],[i3_f],[i4_f],[i5_f]))
    i_ref=np.float64(i_ref)
    i_fuse=np.float64(i_fuse)
    [u,v,w]=i_ref.shape
    ERG=0
    RMSE=np.zeros(u)
    for k in range (0,u):
        temp=i_ref[k,:,:]-i_fuse[k,:,:]
        temp=np.square(temp)
        RMSE[k]=np.sum(temp)
        RMSE[k]=np.sqrt(RMSE[k]/(m*n))
        temp=(np.square(RMSE[k]))/(np.square(np.mean(i_ref[k,:,:])))
        ERG=ERG+temp
    h=48.5 #-----resolution of SAR image / reference MS image 
    l=65.34 #------resolution of resampled MS image
    L=4.0 #--------------No. of fused bands
    ERG=np.sqrt(ERG/L)
    ERG=(100.0)*(h/l)*ERG
    return ERG,RMSE
def rel_bias(i1,i2):
    i1=np.float64(i1)
    i2=np.float64(i2)
    temp1=np.mean(i1) #--------mean of reference MS band
    temp2=np.mean(i2) #-------- mean of fused MS band
    bias=np.abs(temp1-temp2)
    rel_bias=bias/temp1
    return rel_bias

def rel_var(i1,i2):
    i1=np.float64(i1)
    i2=np.float64(i2)
    temp1=np.var(i1) #--------variance of reference MS band
    temp2=np.var(i2) #-------- variance of fused MS band
    var=np.abs(temp1-temp2)
    rel_var=var/temp1
    return rel_var

def rel_sdd(i1,i2):
    i1=np.float64(i1)
    i2=np.float64(i2)
    temp=i1-i2 #---- difference image (reference minus fused band)
    sd=np.std(temp)#-----standard deviation of difference image
    rel_sdd=sd/np.mean(i1)
    return rel_sdd
def SAM(i2_ref,i3_ref,i4_ref,i5_ref,i2_f,i3_f,i4_f,i5_f):
    i_ref=np.concatenate(([i2_ref],[i3_ref],[i4_ref],[i5_ref]))
    i_fuse=np.concatenate(([i2_f],[i3_f],[i4_f],[i5_f]))
    i_ref=np.float64(i_ref)
    i_fuse=np.float64(i_fuse)
    [u,v,w]=i_ref.shape
    s2=np.zeros([v,w])
    s3=np.zeros([v,w])
    for x in range(0,v):
        for y in range(0,w):
            temp1=np.matrix(i_ref[:,x,y]) #-----vector formed from reference bands, for pixel location (x,y)
            temp2=np.matrix(i_fuse[:,x,y])#-----vector formed from fused bands, for pixel location (x,y)  
            temp5=(la.norm(temp1))*(la.norm(temp2)) #----- product of magnitudes of two vectors
            #s2[x,y]=temp1*np.transpose(temp2)#----dot product
            s2[x,y]=np.dot(temp1,np.transpose(temp2))
            s3[x,y]=np.float32(np.divide(s2[x,y],temp5))
            s3[x,y]=np.nan_to_num(s3[x,y])
            s3[x,y]=math.acos(s3[x,y]) #----spectral angle at pixel (x,y) in radians
            s3[x,y]=(s3[x,y]*180)/math.pi
   
    spec_ang=np.mean(s3) #-----global (average) spectral angle
    return spec_ang
def q_calc(temp_ref,temp_f): #----user defined function to calculate local Q index
    mu_ref=np.mean(temp_ref)
    mu_f=np.mean(temp_f)
    v_ref=np.var(temp_ref)
    v_f=np.var(temp_f)
    covar=np.cov(temp_ref,temp_f)
    c=covar[0,1]
    q=(2*mu_ref*mu_f)*(c)
    q=q/((v_ref+v_f)*((mu_ref**2)+(mu_f**2)))
    return q
def uiqi(i2_ref,i3_ref,i4_ref,i5_ref,i2_f,i3_f,i4_f,i5_f):
    i2_ref1=np.concatenate((i2_ref[0:5,:],i2_ref),axis=0)
    i2_ref1=np.concatenate((i2_ref1,i2_ref1[:,2190:2199]),axis=1)
    i3_ref1=np.concatenate((i3_ref[0:5,:],i3_ref),axis=0)
    i3_ref1=np.concatenate((i3_ref1,i3_ref1[:,2190:2199]),axis=1)
    i4_ref1=np.concatenate((i4_ref[0:5,:],i4_ref),axis=0)
    i4_ref1=np.concatenate((i4_ref1,i4_ref1[:,2190:2199]),axis=1)
    i5_ref1=np.concatenate((i5_ref[0:5,:],i5_ref),axis=0)
    i5_ref1=np.concatenate((i5_ref1,i5_ref1[:,2190:2199]),axis=1)
    i2_f1=np.concatenate((i2_f[0:5,:],i2_f),axis=0)
    i2_f1=np.concatenate((i2_f1,i2_f1[:,2190:2199]),axis=1)
    i3_f1=np.concatenate((i3_f[0:5,:],i3_f),axis=0)
    i3_f1=np.concatenate((i3_f1,i3_f1[:,2190:2199]),axis=1)
    i4_f1=np.concatenate((i4_f[0:5,:],i4_f),axis=0)
    i4_f1=np.concatenate((i4_f1,i4_f1[:,2190:2199]),axis=1)
    i5_f1=np.concatenate((i5_f[0:5,:],i5_f),axis=0)
    i5_f1=np.concatenate((i5_f1,i5_f1[:,2190:2199]),axis=1)
    [u,v]=i2_ref1.shape
    r1=np.int32(u/16)
    c1=np.int32(v/16)
    #s2=np.zeros([v,w])
    q2=np.zeros([r1,c1])
    q3=np.zeros([r1,c1])
    q4=np.zeros([r1,c1])
    q5=np.zeros([r1,c1])
    for x in range(0,r1):
         for y in range(0,c1):
             t1=16*x
             t2=t1+16
             t3=16*y
             t4=t3+16
             temp_ref2=i2_ref1[t1:t2,t3:t4]
             temp_f2=i2_f1[t1:t2,t3:t4]
             temp_ref3=i3_ref1[t1:t2,t3:t4]
             temp_f3=i3_f1[t1:t2,t3:t4]
             temp_ref4=i4_ref1[t1:t2,t3:t4]
             temp_f4=i4_f1[t1:t2,t3:t4]
             temp_ref5=i5_ref1[t1:t2,t3:t4]
             temp_f5=i5_f1[t1:t2,t3:t4]
             q2[x,y]=q_calc(temp_ref2,temp_f2) #-----calculation of q index over 16*16 block for band 2
             q3[x,y]=q_calc(temp_ref3,temp_f3)
             q4[x,y]=q_calc(temp_ref4,temp_f4)
             q5[x,y]=q_calc(temp_ref5,temp_f5)
    q2=np.nan_to_num(q2)
    q3=np.nan_to_num(q3)
    q4=np.nan_to_num(q4)
    q5=np.nan_to_num(q5)
    q2_g=np.mean(q2) #----global q index for band 2
    q3_g=np.mean(q3)
    q4_g=np.mean(q4)
    q5_g=np.mean(q5)
    return q2_g,q3_g,q4_g,q5_g
        
        
