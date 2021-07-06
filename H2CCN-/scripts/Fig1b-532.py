import numpy as np
import matplotlib.pyplot as plt

width = 3.487*2
#height = width/1.618 #golden ratio
height = width/1.4
one_col = (width, height)
fig,ax = plt.subplots(figsize=one_col)

x,y = np.loadtxt('../data/CH2CN-6N3CMI1536PESBE.dat',unpack=True)
VMI = np.loadtxt('../data/CH2CN-6N3CMb512.txt')
rows, cols = VMI.shape
trim = 46
#r2 = rows//2 - trim
#c2 = cols//2 - trim
r2 = rows//2
c2 = cols//2


VMI = VMI[trim:-trim, trim:-trim]



plt.plot(x,y,lw=1.5,label='532.1 nm PES')
plt.xlim(11500,16000)
plt.annotate(r'$9^0_1$',(12035,0.0688),ha='center',fontsize='12')
plt.annotate(r'$0^0_0$',(12615,0.9921),ha='center',fontsize='12')
plt.annotate(r'$6^2_1$',(12724,0.08695),ha='center',fontsize='12')
plt.annotate(r'$5^1_1$',(13006,0.3404),ha='center',fontsize='12')
plt.annotate(r'$5^1_0$',(13142,0.2317),ha='center',fontsize='12')
plt.annotate(r'$4^1_0$',(13378,0.0899),ha='center',fontsize='12')
plt.annotate(r'$5^2_0$',(13822,0.2921),ha='center',fontsize='12')
plt.annotate(r'$3^1_0$',(13922,0.1835),ha='center',fontsize='12')
plt.annotate(r'$5^3_1$',(14339,0.1593),ha='center',fontsize='12')
plt.annotate(r'$5^3_0$',(14503,0.1593),ha='center',fontsize='12')
plt.annotate(r'$5^4_0$',(15147,0.0869),ha='center',fontsize='12')
plt.annotate(r'$3^2_0$',(15374,0.0839),ha='center',fontsize='12')

plt.xticks([12000,13000,14000,15000,16000],fontsize=12)
plt.yticks([])
plt.setp(ax.spines.values(), linewidth=1.5)
ax.tick_params(width=1.5)
plt.ylabel('intensity (arb. u.)',fontsize=14)
plt.xlabel('electron binding energy (cm$^{-1}$)',fontsize=14)
plt.minorticks_on()
fig.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.98)
plt.legend()

#sub_axes=plt.axes([.65,.6,.35,.35])
#sub_axes.imshow(VMI, extent=[-r2, r2, -c2, c2], vmin=0, cmap='jet', vmax=VMI.max()/1)
#sub_axes.set_xticks([0, 256])
#sub_axes.set_xticklabels(["0", "256"], fontsize="smaller")
#sub_axes.set_yticks([])
#sub_axes.set_yticklabels([])

plt.savefig('Fig1b-532.pdf')
plt.show()
