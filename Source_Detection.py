#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 13:38:30 2021

@author: alexwardley
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.stats import norm
import scipy.stats as stats
from numpy import unravel_index
import math

#open file and get data
hdulist = fits.open("/Users/aliyah/Downloads/A1_mosaic.fits")
hdulist.info()
#hdulist[0].header()
image_data = hdulist[0].data
hdulist.close()
x_values = image_data
#print(image_data[:,1])



x_values = x_values[ x_values <= 3600]
x_value = x_values[ x_values <= 3450]
xv = x_value[ x_value >= 3390 ]
n,bins,patches= plt.hist(x_values,bins=3600)
mu, sigma= norm.fit(xv )
y=norm.pdf( bins, mu , sigma)
plt.figure(1)
plt.plot(bins, 10500000 * y)
plt.xlim([3300,3600])
plt.show()
print(mu, sigma)


#mean background count


plt.figure(2)
plt.imshow(image_data, cmap='gray')
plt.colorbar()

#%%
xmask = np.ma.make_mask(image_data,copy=True,shrink=True,dtype=np.bool)
xmask[:]=False#all values remain
xmask[0:206,:]=True
xmask[4405:,:]=True
xmask[:,0:206]=True
xmask[:,2364:]=True

#%%
#masking the large star
xmask[:,1426:1451]=True#vertical bleeding
xmask[426:483,1104:1652]=True#horizontal bleeding 1
xmask[314:371,1018:1704]=True#horizontal bleeding 2
xmask[1391:1475,266:218]=True#horizontal bleeding3
#other horizontal bleeding is masked by edge masking

xmask[2223:2356,901:909]=True#bleeding of smaller stars(vertical)
xmask[3203:3419,772:780]=True
xmask[2704:2836,969:978]=True
xmask[3708:3803,2132:2137]=True

r_bstar=219#raduis of the large star
bstar_pos=(3213,1449)#position of the large star
for x in range(-219,220):#loop for x coordinates in range of the negative of the radius to its positive
    for y in range(-219,220):#loop for y coordinates in the same range, both loops create a square of points
        d_r=dr=np.sqrt(np.add(np.abs(x)**2,np.abs(y)**2))#checking if the length to (x,y) lies within the circle of the stars radius
        if d_r<=r_bstar:
            xmask[bstar_pos[0]+x,bstar_pos[1]+y]=True#masking stars within the radius of the star
image_data = np.ma.masked_array(image_data, mask=xmask)

    #%%        
r_star=np.array([33,19,43,31,18,26,27,23,16,32,21,19,25,17,17,24,19,15,16,17,14,20,20,19,15,37,19,15,15])#radii of foreground stars
star_pos=[(2288, 906),(1495, 634),(3321, 777),(2778,975),(1938,677),(4099,562),(1427,2091),(4033,1458),(1654,969),(3758,2135),(4399,1315),(4330,1365),(2980,1419),(3373,896),(3846,2281),(2310,2132),(844,211),(1561,1934),(289,1799),(2431,360),(3598,367),(2301,448),(578,1775),(2264,705),(4315,579),(3297,2261),(2963,1495),(4223,583),(3207,857)]#star positions
for i,j in zip(r_star,star_pos):#loop for each star in the lists
    for x in range(-i,i+1):#loop for x coords in range of the negative of stars radius to its positive
        for y in range(-i,i+1):#loop for y coords in the same range, both loops create a square of points
            d_r=np.sqrt(np.add(np.abs(x)**2,np.abs(y)**2))#length to each (x,y)
            if d_r<=i:#check if length is within stars radius
                k=[j[0]+x,j[1]+y]#k is the coordinate of the point in the image
                if xmask[k[0],k[1]]==False:
                    xmask[k[0],k[1]]=True#masking stars within the radius of the star
image_data = np.ma.masked_array(image_data, mask=xmask)




#%%


plt.figure(3)
plt.imshow(image_data, cmap='gray')#plotting the data to produce an image
plt.colorbar()

x_val=[]
y_val=[]
x_val.append(3213)# adding the x coordinate of the large star to the list
y_val.append(1449)#adding its y coordinate
for i in star_pos:# appending the x,y coordinates for each star to lists
    x_val.append(i[0])
    y_val.append(i[1])
plt.figure(3)
plt.plot(y_val,x_val,'x')#plots the positions of the masked stars 





#%%
r_ap=6 #aperture radius
r_bkg = 9#1.5*r_ap #radius of circle used to determine local backgound count
thresh= mu + 8*sigma

max_pos = np.unravel_index(np.argmax(image_data, axis=None), image_data.shape)#coordinates of max
max_= image_data[max_pos]  #maximum pixel value
brightnesses=[]
positions=[]
b_err=[]
while max_>thresh:
    max_pos = np.unravel_index(np.argmax(image_data, axis=None), image_data.shape)#coordinates of max
    max_= image_data[max_pos]
    positions.append(max_pos)
    br=[]#list of pixels values within radius
    br1=[]
    for x in range(-r_ap,r_ap+1):
            for y in range(-r_ap,r_ap+1):
                dr=np.sqrt(np.add(np.abs(x)**2,np.abs(y)**2))
                if dr<=r_ap:
                    j=[max_pos[0]+x,max_pos[1]+y]
                    if xmask[j[0],j[1]]==False:
                        br.append(image_data[j[0],j[1]])
                        xmask[j[0],j[1]]=True
    image_data = np.ma.masked_array(image_data, mask=xmask)
    for x1 in range(-r_bkg, r_bkg+1):#x  coords for annulus
        for y1 in range(-r_bkg, r_bkg+1):#y coords
                dr_b=np.sqrt(np.add(np.abs(x1)**2,np.abs(y1)**2))
                if dr_b<=r_bkg:#checking if point is within annulus outer radius
                    j=[max_pos[0]+x1,max_pos[1]+y1]
                    if xmask[j[0],j[1]]==False:

                 #       xmask[j[0],j[1]]=True
                         br1.append((image_data[j[0],j[1]]))
                         print(image_data[j[0],j[1]])
                    #    print(j[0],j[1])
              #      else:
               #         xmask[j[0],j[1]]==True
    backbright = np.median(br1)
    print(backbright)
 #   print(max_pos)
    pixels=len(br) #no. of pixels in radius 
    bright=np.subtract(sum(br),(pixels*backbright))
    print(bright)
    brightnesses.append(bright)
    b_err_=np.sqrt(sum(br))
    bbkg_err = (np.std(br1)/(np.sqrt(len(br1))))
    print(bbkg_err)
 #   bright=[]
    b_error = np.sqrt((b_err_)**2+(bbkg_err)**2)
    b_err.append(b_error)
    
#    backbright = np.median(br1)                    
#    pixels=len(br) #no. of pixels in radius 
#    bright=np.subtract(sum(br),(pixels*backbright))
#    brightnesses.append(bright)
#    b_err_=np.sqrt()
#    stats.stdev()
#    b_lbkg=np.sqrt(backbright)
#    bright=[]
#    for i in br:
 #       j=i-backbright
  #      bright.append(j)
   # brightnesses.append(sum(bright))
  #  err_b=[](k))
    #brightnesses_err=[]
  #  for k in br:
  #      err_b.append(np.sqrt
  #  total_err_b=sum(err_b)
  #  brightnesses_err.append(total_err_b)
    image_data = np.ma.masked_array(image_data, mask=xmask)
data=np.array([positions,brightnesses]).T    
no_of_sources = len(positions)

#%%
#brightnesses=[-1,-2-3,1,2,3,4]
#b_err=[1,1,1,1,1,1,1]
#positions=[0,1,2,3,4,5,6]



#brightnesses=np.asarray([brightnesses])

pos_x=[]
pos_y=[]
for i in positions:
    pos_x.append(i[0])
    pos_y.append(i[1])
    
    
negative= np.where(np.array(brightnesses)<=0)
#%%
b=np.delete(brightnesses,negative)
b_err=np.delete(b_err,negative)
pos_x=np.delete(pos_x,negative)
pos_y=np.delete(pos_y,negative)
    
positions=np.array([pos_x,pos_y]).T


#%%

Zpi = hdulist[0].header['MAGZPT'] 
m=[]
for i in range(len(b)):
    m_ = Zpi - 2.5*np.log10(b[i])
    m.append(m_)
    
#m_err_pos =np.subtract(np.log(b+b_err),np.log(b))
#m_err_neg =np.subtract(np.log(b-b_err),np.log(b))
#m_err = [2.5*m_err_pos,2.5*m_err_neg]

m_err=[]
for i in range(len(b)):
    mep=np.subtract(np.log10(b[i]+b_err[i]),np.log10(b[i]))
    men=np.subtract(np.log10(b[i]),np.log10(b[i]-b_err[i]))
    me=[2.5*mep,2.5*men]
    m_err.append(me)

print(thresh)
#%%
plt.figure(3)
#image_data = hdulist[0].data
plt.imshow(image_data, cmap='gray')
plt.colorbar()
x=[]
y=[]
#%%
for i in positions:
    y.append(i[0])
    x.append(i[1])
plt.plot(x,y,'.', markersize=3)

#%%
file=open('Catalogue_r_ap'+ str(int(r_ap))+'_thresh'+str(int(thresh))+'.txt','w')
#file=open('Catalogue_r_ap'+ str(int(r_ap))+'.txt','w')
for index in range(0,len(positions)):
    file.write(str(positions[index][0])+','+str(positions[index][1]) + ", " +str(brightnesses[index]) +',' + str(m[index])+',' + str(m_err[index]) + " " + "\n")
#file.close()


