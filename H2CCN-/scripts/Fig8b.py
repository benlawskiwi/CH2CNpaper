from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#width=3.33*2
#height=width/1.4
#one_col=(width,2*height)
#fig,ax = plt.subplots(figsize=one_col)

height=3.487*2
width = height/1.4
one_col = (width,height)
fig,ax = plt.subplots(figsize=one_col)

kb = 8.617333262145e-5
ev = 8065.74

#Anion Ground State
Add = 9.29431
Bdd = 0.338427
Cdd = 0.327061
Djdd = 1.65e-7
Djkdd = 1.297e-5
Dkdd = 1.0028e-3

#Dipole Bound State
Ad = 9.51035  #Lineberger
Bd = 0.341049
Cd = 0.328764
Djd = 2.19e-7
Djkd = 1.499e-5
Dkd = 0.7546e-3

#Neutral State (Microwave study - Yamamoto)
m = 1/29979.2458
An = 284981*m
Bn = 10426.765*m
Cn = 9876.035*m
Djn = 0.0040515*m
Djkn = 0.41573*m
Dkn = 23.536*m

#Detachment Threshold
D0 = 129      #Neumark
DBS = 39
EA = 12468

#Photon Energy
hv = 1e7/801.567
tol = 0.3

#Set J
J = 20

def Gauss (E0,FWHM,E):
    alpha = FWHM/2/sqrt(log(2.0))
    return exp(-((E-E0)/alpha)**2)

def En (J,K,A,B,C,Dj,Djk,Dk):
    Bb = (B+C)/2
    E = Bb*J*(J+1)+(A-Bb)*K**2-Dj*J**2*(J+1)**2-Djk*J*(J+1)*K**2-Dk*K**4
    return E

def Boltz (E,T,J):
    P = (2*J+1)*exp(-(E/ev)/(T*kb))
    return P

#First set two anion example levels

ae1 = En(23,1,Add,Bdd,Cdd,Djdd,Djkdd,Dkdd)
ae2 = En(33,1,Add,Bdd,Cdd,Djdd,Djkdd,Dkdd)
ae0 = En(0,0,Add,Bdd,Cdd,Djdd,Djkdd,Dkdd)
ae3 = En(8,0,Add,Bdd,Cdd,Djdd,Djkdd,Dkdd)

plt.plot((1,3.5),(ae1,ae1),'k')
plt.plot((1,3.5),(ae2,ae2),'k')
plt.plot((1,3.5),(ae0,ae0),'k')
plt.plot((1,3.5),(ae3,ae3),'C7--')
plt.annotate('1, 23',(3.55,ae1),va='center')
plt.annotate('1, 33',(3.55,ae2),va='center')
plt.annotate('0, 0',(3.55,ae0),va='center')
plt.annotate('0, 8',(3.55,ae3+10),color='C7',va='center')
plt.annotate("K'', J''",(3.55,ae2+55),va='center')


de1 = En(24,0,Ad,Bd,Cd,Djd,Djkd,Dkd)+1000-DBS
de2 = En(32,2,Ad,Bd,Cd,Djd,Djkd,Dkd)+1000-DBS
de0 = En(0,0,Ad,Bd,Cd,Djd,Djkd,Dkd)+1000-DBS
de3 = En(8,1,Ad,Bd,Cd,Djd,Djkd,Dkd)+1000-DBS

plt.plot((2.5,3.5),(de1,de1),'C3')
plt.plot((2.5,3.5),(de2,de2),'C3')
plt.plot((2.5,3.5),(de0,de0),'C3')
plt.plot((2.5,3.5),(de3,de3),'C7--')
plt.annotate('0, 24',(3.55,de1))
plt.annotate('2, 32',(3.55,de2))
plt.annotate('0, 0',(3.55,de0))
plt.annotate('1, 8',(3.55,de3),color='C7')
plt.annotate("K, J",(3.55,de2+55),va='center')

ne1 = En(21,0,An,Bn,Cn,Djn,Djkn,Dkn)+1000
ne2 = En(31,0,An,Bn,Cn,Djn,Djkn,Dkn)+1000
ne0 = En(0,0,An,Bn,Cn,Djn,Djkn,Dkn)+1000
ne4 = En(22,0,An,Bn,Cn,Djn,Djkn,Dkn)+1000
ne5 = En(32,0,An,Bn,Cn,Djn,Djkn,Dkn)+1000

plt.plot((1,2),(ne1,ne1),'C2')
plt.plot((1,2),(ne2,ne2),'C2')
plt.plot((1,2),(ne0,ne0),'C2')
plt.annotate('0, 21',(0.5,ne1))
plt.annotate('0, 31',(0.5,ne2))
plt.annotate('0, 0',(0.5,ne0))
plt.annotate("K', J'",(0.50,ne2+55),va='center')

#Absorption arrows
plt.plot((2.8,2.8),(ae1,de1),'C4')
plt.plot((2.8,2.77),(de1,de1-10),'C4')
plt.plot((2.8,2.83),(de1,de1-10),'C4')

plt.plot((3.0,3.0),(ae2,de2),'C4')
plt.plot((3.0,2.97),(de2,de2-10),'C4')
plt.plot((3.0,3.03),(de2,de2-10),'C4')

plt.plot((3.3,3.3),(ae3,de3),'C7')
plt.plot((3.3,3.27),(de3,de3-10),'C7')
plt.plot((3.3,3.33),(de3,de3-10),'C7')

#Autodetachment arrows
plt.plot((2.5,2),(de1,ne1),'--C4')
plt.plot((2,2.03),(ne1,ne1+10),'--C4')
plt.plot((2,2.03),(ne1,ne1-10),'--C4')

plt.plot((2.5,2),(de2,ne2),'--C4')
plt.plot((2,2.03),(ne2,ne2+10),'--C4')
plt.plot((2,2.03),(ne2,ne2-10),'--C4')

#Direct detachment
e1 = ae0+1e7/804.058-EA+1000
e2 = ae1+1e7/804.058-EA+1000
e3 = ae2+1e7/804.058-EA+1000

plt.plot((1.1,1.1),(ae0,e1),'C9')
plt.plot((1.1,1.13),(e1,e1-10),'C9')
plt.plot((1.1,1.07),(e1,e1-10),'C9')

plt.plot((1,1.2),(e1,e1),'C9--')

#Labels
plt.annotate(r'$^RP_1$',(3.06,848),ha='left',color='C4')
plt.annotate(r'$\lambda = 804.058$nm',(2.72,607),ha='right',color='C4')
plt.annotate(r'$^PR_1$',(2.72,557),ha='right',color='C4')
plt.annotate(r'$^RQ_0$',(3.36,588),ha='left',color='C7')

plt.annotate(r'$^QN_0$',(2.25,de1+10),ha='center',color='C4')
plt.annotate(r'$^OP_2$',(2.25,de2+10),ha='center',color='C4')

plt.annotate(r'Direct Detachment',(1.16,800),color='C9')
plt.annotate(r'$\lambda=804.058$nm',(1.16,750),color='C9')

plt.annotate(r'Autodetachment',(2.25,1400),ha='center',color='C4')

plt.annotate('Anion',(2.25,-50),ha='center')
plt.annotate('Neutral',(1.5,1450),ha='center',color='C2')
plt.annotate('Dipole Bound State',(3.0,1450),ha='center',color='C3')

#DBS arrow
plt.plot((1.7,1.7),(ne0,de0),'k')
plt.plot((1.7,1.73),(ne0,ne0-7),'k')
plt.plot((1.7,1.67),(ne0,ne0-7),'k')
plt.plot((1.7,1.73),(de0,de0+7),'k')
plt.plot((1.7,1.67),(de0,de0+7),'k')
plt.annotate(r'$39$cm$^{-1}$',(1.75,(de0)))

#Plot settings
plt.xlim(0.45,4.15)
plt.ylim(-100,1500)
plt.xticks([])
ytick = [0,200,400,1000+(12400-EA),1000-DBS,1000,1000+(12600-EA),1000+(12800-EA)]
ylab = ['0','200','400','12400','DBS','EA','12600','12800']
plt.setp(ax.spines.values(), linewidth=1.5)

plt.yticks(ytick,ylab)
plt.ylabel(r'Rotational Energy (cm$^{-1})$',fontsize=14)
fig.subplots_adjust(left=0.20, bottom=0.02, right=0.95, top=0.98)
#plt.title('Detachment pathways at $\lambda=804.058$nm')
plt.savefig('Fig8b.pdf',dpi=400,bbbox_inches='tight')

plt.show()
