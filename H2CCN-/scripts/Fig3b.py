import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

xa,ya = np.loadtxt("2.csv",unpack=True,delimiter=',')
f = interp1d(xa,ya,kind='cubic')
xn = np.arange(0.21,0.79,0.01)
yn = 46600*(xn-0.5)**2+12268

height=3.487*2/1.4
width = height/1.4
one_col = (width,height)


fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(one_col),gridspec_kw={'height_ratios': [2,1]})

ax2.plot(xn,f(xn)-64,'k')
ax.plot(xn,yn,'k')
ax.plot(xn,yn-100,'--',color='silver')
ax2.set_ylim(0,400)
ax.set_ylim(12000,15450)

ax2.annotate(r'$\tilde{X}^1A_1$',(0.2,320))
ax.annotate(r'$\tilde{X}^2B_1$',(0.2,13200))
ax.annotate('DBS',(0.2,12710),color='grey')
ax2.plot((0.28,0.72),(90,90),'k',lw='1')
ax2.plot((0.225,0.775),(220,220),color='k',lw='1')
ax2.annotate(r'$1$',(0.78,220))
ax2.annotate(r'$0$',(0.78,90))

ax.plot((0.422,0.578),(12558,12558),color='k',lw='1')
ax.plot((0.422,0.578),(12458,12458),'--',lw='1',color='silver')
ax.plot((0.359,0.641),(13217,13217),color='k',lw='1')
ax.plot((0.359,0.641),(13117,13117),'--',lw='1',color='silver')
ax.plot((0.314,0.686),(13898,13898),color='k',lw='1')
ax.plot((0.314,0.686),(13798,13798),'--',lw='1',color='silver')
ax.plot((0.28,0.72),(14573,14573),color='k',lw='1')
ax.plot((0.28,0.72),(14473,14473),'--',lw='1',color='silver')
ax.plot((0.248,0.752),(15261,15261),color='k',lw='1')
ax.plot((0.248,0.752),(15161,15161),'--',lw='1',color='silver')
ax.annotate(r'$0$',(0.78,12558))
ax.annotate(r'$1$',(0.78,13217))
ax.annotate(r'$2$',(0.78,13898))
ax.annotate(r'$3$',(0.78,14573))
ax.annotate(r'$4$',(0.78,15261))

#Transitions Below
ax.arrow(0.45,90,0,15071,length_includes_head=True,head_width=0.01,head_length=50,color='C4')
ax2.arrow(0.45,90,0,15071,length_includes_head=True,head_width=0.01,head_length=50,color='C4')
ax.arrow(0.5,15161,0,-588,length_includes_head=True,head_width=0.01,head_length=50,color='C4')
ax2.arrow(0.5,15161,0,-588,length_includes_head=True,head_width=0.01,head_length=50,color='C4')
ax.arrow(0.55,14867,0.1,0,length_includes_head=True,head_width=50,head_length=0.01,color='C0')
ax.annotate(r'$e^-$',(0.59,14900),color='C0',ha='left',fontsize=8)
ax.annotate(r'637 cm$^{-1}$',(0.55,14687),ha='left',color='C0',fontsize=8)
ax2.annotate(r'$hv$ = 15,126 cm$^{-1}$',(0.5,400),ha='left',color='C4',fontsize=8)
ax2.set_yticks([0,200,400])
ax.set_yticks([12000,13000,14000,15000])
ax.annotate(r'Potential Energy (cm$^{-1})$',xy=(-5.8, 0.0), rotation=90, xycoords='axes fraction', xytext=(-0.30, 0.2),fontsize=12, ha='center', va='center')
ax.annotate(r'$v_5$',xy=(-5.8, 0.0), rotation=0, xycoords='axes fraction', xytext=(0.95, 1.05),fontsize=12, ha='center', va='center')

ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off')  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
d = .015
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
fig.subplots_adjust(left=0.25, bottom=0.05, right=0.95, top=0.95)
plt.xticks([], [])
plt.savefig("Fig3b.pdf",dpi=400,bbbox_inches="tight")
plt.show()

