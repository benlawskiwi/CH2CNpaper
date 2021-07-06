from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

kb = 8.617333262145e-5
ev = 8065.74

x,y = loadtxt('../data/CH2CN-799-512IPESBE.dat',unpack=True)

width = 3.487*2
height = width/1.4
one_col=(width,height)

fig,ax = plt.subplots(figsize=one_col)

#Fitting Parameters
T = 150
amp = 1
EA = 12468
#DBS = 39
DBS = 39.08
FWHM = 1 
err = -1.073
#err = 0
#ofs = 0
ofs = -1.2

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
#D0 = 0

#Photon Energy
hv = 1e7/801.567
tol = 0.3

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

def spectrum(Ev,T,FWHM,amp,ofs,DBS,EA):
    #Set up spectrum
    #Ev = arange(EA-200,EA+200,0.1)
    spec = zeros(size(Ev))
    spec2 = zeros(size(Ev))
    linelist = {}
    linelist2 = {}
    nu = EA-DBS+err
    nu = 12428.665
    #Selection rules/intensities                    (Jdd,Kdd) -> (Jd,Kd)
    for Jdd in range(0,52):
        for dJ in (-1,0,1):
            Jd = Jdd + dJ                           #dJ=0,=/-1
            if Jd <0 :continue
            for Kdd in range (0,Jdd+1):
                Edd = En(Jdd,Kdd,Add,Bdd,Cdd,Djdd,Djkdd,Dkdd)
                I = Boltz(Edd,T,Jdd)
                if Kdd%2 !=0:
                    I *=3                           #Spin statistics 3:1
                for dK in (-1,1):                   #dK=+/-1
                    Kd = Kdd + dK
                    if Kd < 0: continue
                    if Kd > Jd : continue
                    if dJ == 1 and dK == 1: HL = (Jdd+2+Kdd)*(Jdd+1+Kdd)/((Jdd+1)*(2*Jdd+1))
                    if dJ == 1 and dK == -1: HL = (Jdd+2-Kdd)*(Jdd+1-Kdd)/((Jdd+1)*(2*Jdd+1))
                    if dJ == 0 and dK == 1: HL = (Jdd+1+Kdd)*(Jdd-Kdd)/(Jdd*(Jdd+1))
                    if dJ == 0 and dK == -1: HL = (Jdd+1-Kdd)*(Jd+Kdd)/(Jdd*(Jdd+1))
                    if dJ == -1 and dK == 1: HL = (Jdd-1-Kdd)*(Jdd-Kdd)/(Jdd*(2*Jdd+1))
                    if dJ == -1 and dK == -1: HL = (Jdd-1+Kdd)*(Jdd+Kdd)/(Jdd*(2*Jdd+1))
                    I *= HL #Honl-London Factors - perpendicular transition of a prolate symmetric top
                    Ed = En(Jd,Kd,Ad,Bd,Cd,Djd,Djkd,Dkd)
                    if Ed < D0: continue  #These lines won't autodetach
                    dE = Ed-Edd+nu
                    if dE < hv-tol : continue
                    if dE > hv+tol : continue
                    Jarr = arange(-5,2,1)
                    Jamp=1
                    for dJn in Jarr:
                        Jn = Jd+dJn
                        if dJn ==0 :
                            I *= 1      #propensity rule - lowest delta J have largest intensity
                        if dJn ==-1:
                            I*=2
                        if dJn <=-3:
                            Jamp=0.2
                        if dJn >-3:
                            Jamp = 1
                        for dKn in (-2,0):
                            Kn = Kd + dKn
                            if Kn < 0 : continue
                            if Kn > Jn : continue
                            Edn = En(Jn,Kn,An,Bn,Cn,Djn,Djkn,Dkn)
                            #KE = Edn+DBS-Ed
                            KE = Ed-(Edn+DBS)
                            #if KE <= 25 and KE > 0:
                            #    I *=1/5*KE**(1/2)
                            BE = hv-KE+ofs
                            linelist[(Jd,Jn,Kd,Kn,dJn,dKn)] = (BE,amp*I)
                            #linelist[(Jd,Jn,Kd,Kn,Ed,Edn)] = (BE,amp*I)
                            spec += Jamp*amp*I*Gauss(BE,FWHM,Ev)
    return (spec,linelist,Ev)

def fit(Ev,EA):            #Change amp to whatever parameter(s) you want to fit. Currently set for amp, T, FWHM, but easy to extend to others
    spec,linelist,Ev = spectrum(Ev,T,FWHM,amp,ofs,DBS,EA)
    return spec

#Curve Fit
popt, pcov = curve_fit(fit,x,y,bounds=(12460,12480))
print(popt)
print(sqrt(diag(popt)))


spec,linelist,Ev = spectrum(x,T,FWHM,amp,ofs,DBS,EA)
linea = []
posa = []
#print('(Jdd,Jd,Kdd,Kd) (E,I)')
for line,posAmp in sorted(linelist.items(),key=lambda x: x[1][0]):
    #print(line,posAmp)
    linea.append(line)
    posa.append(posAmp)

savetxt('linelist.dat',column_stack((linea,posa)))

#normalise
spec *= 1/max(spec)

#normalise fit
#ynnn = fit(x,*popt)
#normnnn = 1/max(ynnn)

#plt.plot(x,y,label='801.567nm PES')
plt.plot(Ev,spec,'C1',label='Autodetachment Model')
plt.plot(x,y,'C0',label='801.567nm PES')
#plt.plot(x,fit(x,*popt)*normnnn,'C3')
#plot params
#for i in range (0,size(popt)):
#    x0 = 12350
#    y0 = 500-50*i
#    perr = sqrt(diag(pcov[i]))
#    perrs=round(perr[0][0],3)
#    st = 'param = '+str(round(popt[i],3))+' +/- '+str(perrs)
#    plt.annotate(st,(x0,y0))

#savetxt('threshold-set.dat',column_stack((x,fit(x,*popt))))

plt.xlabel(r'Electron Binding Energy (cm$^{-1}$)',fontsize=14)
plt.ylabel('Intensity (arb. u.)',fontsize=14)
plt.xlim(12300,12500)
plt.yticks([])
plt.xticks([12300,12350,12400,12450,12500],fontsize=12)
plt.legend(loc=2)
plt.annotate(r'$^oP_3$',(12425,0.66),color='C0',fontsize=12)
plt.annotate(r'$^oQ_3$',(12428,0.986),color='C0',fontsize=12)
plt.annotate(r'$^oR_3$',(12452,0.448),color='C0',fontsize=12)
plt.annotate(r'$^oP_3: \Delta K=-2, \Delta J-1, K_1 \leftarrow K_3$',(12311,0.65))
plt.annotate(r'$^oQ_3: \Delta K=-2, \Delta J=0, K_1 \leftarrow K_3$',(12311,0.60))
plt.annotate(r'$^oR_3: \Delta K=-2, \Delta J+1, K_1 \leftarrow K_3$',(12311,0.55))
plt.annotate(r'$|\Delta J| > 1$',(12384,0.24),color='C0')
#plt.title('autodetachment 801.567nm resonance')

plt.setp(ax.spines.values(), linewidth=1.5)

plt.savefig('Fig7a.pdf',pdf=400,bbbox_inches='tight')
plt.show()

#plt.plot(x,y)
#plt.plot(x,test(x,*popt))
#plt.show()
