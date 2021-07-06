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


ak = []
ae = []

for Kdd in range (0,5):
    Edd = En(J,Kdd,Add,Bdd,Cdd,Djdd,Djkdd,Dkdd)
    ak.append(Kdd)
    ae.append(Edd)

for i in ak:
    plt.plot((1,3.5),(ae[i],ae[i]),'k')
    plt.annotate(str(ak[i]),(3.55,ae[i]),va='center')
plt.annotate("K''",(3.55,ae[4]+15))


dk = []
de = []

for Kd in range (0,5):
    Ed = En(J,Kd,Ad,Bd,Cd,Djd,Djkd,Dkd)+400-DBS
    dk.append(Kd)
    de.append(Ed)

for i in dk:
    plt.plot((2.5,3.5),(de[i],de[i]),'C3')
    plt.annotate(str(dk[i]),(3.55,de[i]),va='center')

plt.annotate("K",(3.55,de[4]+15))

nk = []
ne = []

for Kn in range (0,5):
    Ene = En(J,Kn,An,Bn,Cn,Djn,Djkn,Dkn)+400
    nk.append(Kn)
    ne.append(Ene)

for i in dk:
    plt.plot((1,2.0),(ne[i],ne[i]),'C2')
    plt.annotate(str(nk[i]),(0.70,ne[i]),va='center')

plt.annotate("K'",(0.70,ne[4]+15))

plt.plot((3,3),(ae[2],de[3]),'C4')
plt.plot((3,2.97),(de[3],de[3]-7),'C4')
plt.plot((3,3.03),(de[3],de[3]-7),'C4')

plt.plot((2.5,2),(de[3],ne[1]),'--C4')
plt.plot((2,2.03),(ne[1],ne[1]+7),'--C4')
plt.plot((2,2.04),(ne[1],ne[1]-2),'--C4')

e1 = ae[1]+1e7/801.567-EA+400
e2 = ae[2]+1e7/801.567-EA+400
e3 = ae[3]+1e7/801.567-EA+400

plt.plot((1.1,1.1),(ae[1],e1),'C9')
plt.plot((1.1,1.13),(e1,e1-7),'C9')
plt.plot((1.1,1.07),(e1,e1-7),'C9')

plt.plot((1.3,1.3),(ae[2],e2),'C9')
plt.plot((1.3,1.33),(e2,e2-7),'C9')
plt.plot((1.3,1.27),(e2,e2-7),'C9')

plt.plot((1.5,1.5),(ae[3],e3),'C9')
plt.plot((1.5,1.53),(e3,e3-7),'C9')
plt.plot((1.5,1.47),(e3,e3-7),'C9')

plt.plot((1.2,1.2),(e1,ne[0]),'--C9')
plt.plot((1.2,1.23),(ne[0],ne[0]+7),'--C9')
plt.plot((1.2,1.17),(ne[0],ne[0]+7),'--C9')

plt.plot((1.4,1.4),(e2,ne[1]),'--C9')
plt.plot((1.4,1.43),(ne[1],ne[1]+7),'--C9')
plt.plot((1.4,1.37),(ne[1],ne[1]+7),'--C9')

plt.plot((1.6,1.6),(e3,ne[2]),'--C9')
plt.plot((1.6,1.63),(ne[2],ne[2]+7),'--C9')
plt.plot((1.6,1.57),(ne[2],ne[2]+7),'--C9')

plt.annotate(r'$\Delta K = -1$',(0.9,600),color='C9')
plt.annotate(r'$\Delta K = -2$',(2.1,600),color='C4')
plt.annotate(r'$\lambda = 801.567$nm',(1.64,380),color='C9')
plt.annotate(r'$\lambda = 801.567$nm',(3.14,380),color='C4')

plt.plot((1.7,1.7),(ne[0],de[0]),'k')
plt.plot((1.7,1.73),(ne[0],ne[0]-7),'k')
plt.plot((1.7,1.67),(ne[0],ne[0]-7),'k')
plt.plot((1.7,1.73),(de[0],de[0]+7),'k')
plt.plot((1.7,1.67),(de[0],de[0]+7),'k')
plt.annotate(r'$39$cm$^{-1}$',(1.75,(ne[0]+de[0])/2))
plt.annotate('Anion',(2.25,125),ha='center')
plt.annotate('Neutral',(1.5,698),ha='center',color='C2')
plt.annotate('Dipole Bound State',(3.0,698),ha='center',color='C3')


plt.xlim(0.5,4)
plt.ylim(110,720)
plt.xticks([])
plt.setp(ax.spines.values(), linewidth=1.5)
ytick = [ae[0],ae[0]+50,ae[0]+100,ae[0]+150,ae[0]+400-68,400-DBS+ae[0],400+ae[0],ae[0]+400-68+100,ae[0]+600-68]
ylab = ['0','50','100','150','12,400','DBS','EA','12,500','12,600']
plt.yticks(ytick,ylab,fontsize=12)
plt.ylabel(r'Rotational Energy (cm$^{-1})$',fontsize=14)
#plt.title('Energy diagram for J=20')
fig.subplots_adjust(left=0.20, bottom=0.02, right=0.95, top=0.98)
plt.savefig('Fig7b.pdf',dpi=400,bbbox_inches='tight')
plt.show()
    
