import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
import glob

files = sorted(glob.glob(f'../CH2CN-/*nm/*/CH2CN-*PES_qA.dat'))


wms = 5
bes = 0.1

wm = np.arange(651,770,bes)
#Set up 532.1nm as reference
x,y = np.loadtxt(files[1],unpack=True)
xmin = np.arange(9010,x[0],wms)
xmax = np.arange(x[-1],19090,wms)
grid = np.arange(x[0],x[-1],wms)
fd = interp1d(x,y)
yd = fd(grid)
eBE5 = np.concatenate([xmin,grid,xmax])
spec5 = np.concatenate([xmin*0,yd,xmax*0])
f532 = interp1d(eBE5,spec5)

x2,y2 = np.meshgrid(eBE5,wm, sparse=True)
#plt.plot(eBE5,f532(eBE5))
#plt.show()

d2 = eBE5
l2 = ['head']
#Test on other spectra
for i in range(3,122):
    x,y = np.loadtxt(files[i],unpack=True)
    l = str(files[i]).split("CH2CN-")[2]
    lab = l.split("nm")[0]
    #print(lab)
    xmin = np.arange(9000,x[0],wms)
    xmax = np.arange(x[-1],20000,wms)
    grid = np.arange(x[0],x[-1],wms)
    fd = interp1d(x,y)
    yd = fd(grid)
    eBE = np.concatenate([xmin,grid,xmax])
    spec = np.concatenate([xmin*0,yd,xmax*0])
    f = interp1d(eBE,spec)
    fy = f(eBE5)-f532(eBE5)
    #fy = f(eBE5)
    yy = []
    ofs = float(lab)*0.1*0
    for j in fy:
        if j <0:
            yy.append(float(0)+ofs)
        if j >=0:
            yy.append(float(j)+ofs)
    plt.plot(eBE5,yy)
    d2 = np.row_stack((d2,yy))
    l2.append(str(lab))
plt.show()

#Can now start interpolating in the wm y direction
#for i in range(0,np.size(eBE5)):
#    if eBE5[i] > 14486:
#        if eBE5[i] < 14487:
#            print(i)

n = np.shape(d2)
print(n)
l2 = l2[1:]
d2 = d2[1:,::]

fig=plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

d3 = wm
nj = int(n[1]/10)
count = -10
for j in range(0,n[1]):
    if j%nj==0:
        count +=10
        #print('--'+str(j/100)+'%')
        print('--'+str(count)+'%')
    y = d2[:,j]
    x = [float(i) for i in l2]
    f = interp1d(x,y)
    yf = f(wm)
    d3 = np.column_stack((d3,yf))

d3 = d3[::,1:]

print('Size of arrays')
print('wm: '+str(np.size(wm)))
print('eBE5: '+str(np.size(eBE5)))
print('Z: '+str(np.shape(d3)))

cc = plt.pcolor(eBE5,wm,d3,cmap='jet',vmax=0.2)
#plt.plot(x,y)
#plt.plot(wm,yf)
ax1.set_xlim(12000,15500)
ax1.set_xlabel(r'eBE (cm$^{-1}$)',fontsize=14)
ax1.set_ylabel(r'wavelength (nm)', fontsize=14)
#for i in range(0,n[1]):
#    y = d2[:,i]
#    x = [float(i) for i in l2]
#    plt.plot(x,y)
new_ticks = [12480,13000,13800,14375,15141]
new_labels = [r'$0^0_0$',r'$5^1_1$',r'$5^2_0$',r'$5^3_1$',r'$5^4_0$']

ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_ticks)
ax2.set_xticklabels(new_labels)
ax2.set_xlabel('Vibrational mode',fontsize=14)
#plt.plot(x,y)
plt.savefig('color-map.pdf',dpi=400,bbbox_inches='tight')
plt.show()
#for i in files:
    #print(i)
    #x,y = np.loadtxt(i,unpack=True)
    #plt.plot(x,y)

#plt.show()
