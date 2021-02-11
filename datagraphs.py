#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 10:27:30 2021

@author: alexwardley
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats

hdulist = fits.open("/Users/alexwardley/Downloads/A1_mosaic.fits")
hdulist.info()
fits.getheader("/Users/alexwardley/Downloads/A1_mosaic.fits",0)
file = pd.read_csv("Catalogue_r_ap10_thresh3516.txt",header=None)


pos=[]
b=[]
m_err=[]
m=[]
for i in range(len(file[0])): #extracts file data into lists
    p=(file[0][i],file[1][i])
    pos.append(p)
    b.append(file[2][i])
    m_er=(file[4][i],file[5][i])
    m_err.append(m_er)
    m.append(file[3][i])
m=np.sort(m) #sorts into ascending order
b=np.sort(b)


N=[]
for n in range(1,len(m)):#ignore the first m
    M=m[:n]
    N1=len(M) # the number of m greater than the chose m value
    N.append(N1) # creates a list of the number above
logN = np.log10(N) # logs the N values
m=m[1:]#removes first m value


#plt.plot(m,logN, 'x') # plots all the data points


lowlim=np.where(m<=12.3)[0][-1]#x value upper limit for fit
highlim=np.where(m>=16)[0][0]#lower limit

N_err=np.sqrt(N) # error in the N values 
yerr_pos=np.subtract(np.log10(N+N_err),logN) # error in the y values in the positve y-axis
yerr_neg=np.subtract(logN,np.log10(N-N_err)) # error in the y values in the negativew y-axis
yerr=np.array([yerr_neg,yerr_pos])


m_=m[lowlim:highlim]#limits data to linear trend region
logN_=logN[lowlim:highlim]
N=N[lowlim:highlim]





#plt.errorbar(m,logN,yerr=yerr,fmt='none')
yerr_neg=yerr_neg[lowlim:highlim]

fit=np.polyfit(m_,logN_,deg=1,w=1/yerr_neg,cov=True)
slope,intercept,r_val,p_val,std_err=stats.linregress(m_,logN_)
plt.plot(m,logN,'x')
plt.plot(m_,np.poly1d(fit[0])(m_))#plots linear fit
plt.xlabel('Magnitude')
plt.ylabel('log10(N(m))')
plt.legend(('Data points','Line of Best Fit'))
plt.figure(1)
plt.grid()
plt.show()



no=500

l=np.empty([len(logN_),no])
for i in range(len(logN_)):
    for j in range(no):
        t=np.random.normal(logN_[i],yerr_neg[i],no)
        l[i,j]=t[j]
slopes=[]
for k in range(no):
    logN_2=l[:,k]
  #  print(logN_2)
    N_err=np.sqrt(N)
    yerr_neg_=np.subtract(logN_2,np.log10(N-N_err))
    fit2=np.polyfit(m_,logN_2,deg=1,w=1/yerr_neg_,cov=True)
    slope2,intercept2,r_val2,p_val2,std_err2=stats.linregress(m_,logN_2)
    slopes.append(slope2)
    plt.figure(2)
    plt.grid()
    plt.plot(m_,np.poly1d(fit2[0])(m_))
    plt.xlabel('Magnitude')
    plt.ylabel("logN(m)")
    plt.show()
print('grad=',slope)
print(np.mean(slopes))
print((r_val)**2)
print('std_err',std_err)
#r=np.abs(np.subtract(slope,np.mean(slopes)))
#print(r)

print((max(slopes)-min(slopes))/2)
print(np.std(slopes))


#%%
bin_centers = []

hist, edges = np.histogram(m,30) #  puts the data into bins with the number of bins base on the number in the bracket   
for i in range(1,len(edges)):
    midpoint =((edges[i-1]+edges[i])/2)
    bin_centers.append(midpoint)# gives the centers of the bins 

N3 = np.cumsum(hist) # gives the number of magnitudes in each bin plus the number of magnitudes in the previous bins
logN3 = np.log10(N3)

N_err=np.sqrt(N3) # error in the N values 
yerr_pos_=np.subtract(np.log10(N3+N_err),logN3) # error in the y values in the positve y-axis
yerr_neg_=np.subtract(logN3,np.log10(N3-N_err)) # error in the y values in the negativew y-axis
yerr_=np.array([yerr_neg_,yerr_pos_])

plt.figure(3)
plt.errorbar(bin_centers,logN3,yerr=yerr_,fmt='none')
data, =plt.plot(bin_centers,logN3, 'x')
fit3=np.polyfit(bin_centers,logN3,deg=1,cov=True)
slope3,intercept3,r_val3,p_val3,std_err3=stats.linregress(bin_centers,logN3)



lowlim=np.where(np.array(bin_centers)<=12.3)[0][-1]#x value upper limit for fit
highlim=np.where(np.array(bin_centers)>=16)[0][0]#lower limit
bin_centers=bin_centers[lowlim:highlim]#limits data to linear trend region
logN3=logN3[lowlim:highlim]

plt.xlabel('Magnitude')
plt.ylabel('log10(N(m))')

bestfit, = plt.plot(bin_centers,np.poly1d(fit2[0])(bin_centers))
#leg = plt.legend(handles=[bestfit], loc='lower left')
#plt.gca().add_artist(leg)
plt.legend(('Data points','Line of Best fit '))
plt.grid()
plt.show()
#%%
plt.figure(4)
#plt.plot(bin_centers,logN3, 'x',label='bins')
plt.plot(bin_centers,np.poly1d(fit2[0])(bin_centers),label='bins linear fit')
plt.plot(m_,np.poly1d(fit[0])(m_),label='data linear fit')
#plt.plot(m,logN, 'x',label='data')
plt.legend(loc='upper left')

#%%


