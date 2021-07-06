from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

kb = 8.617333262145e-5
ev = 8065.74

x,y = loadtxt('../data/CH2CN-804.dat',unpack=True)

width = 3.487*2
height = width/1.4
one_col=(width,height)

fig,ax = plt.subplots(figsize=one_col)


#Fitting Parameters
T = 150
amp = 1
EA = 12468
#EA = 12481.3
#DBS = 39
DBS = 39.08
#DBS = 53
FWHM = 1 
ofs = -1.073-4
#ofs = -1.2

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
#D0 = 129      #Neumark
D0=39

#Photon Energy
#hv = 1e7/802.280
hv = 1e7/804.058
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
    nu = EA-DBS
    nu = 12428.665
    #Selection rules/intensities                    (Jdd,Kdd) -> (Jd,Kd)
    for Jdd in range(0,152):
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
                    if Ed < D0:  #These lines won't autodetach
                        #print('no auto!')
                        continue
                    dE = Ed-Edd+nu
                    if dE < hv-tol : continue
                    if dE > hv+tol : continue
                    Jarr = arange(-4,2,1)
                    Jamp=1
                    for dJn in Jarr:
                        Jn = Jd+dJn
                        if dJn ==0 :
                            I *= 1      #propensity rule - lowest delta J have largest intensity
                        if dJn ==-1:
                            I*=2
                        if dJn == -4:
                            Jamp =0.1
                            print('yes')
                        if dJn != -4:
                            Jamp=1
                            print('no')
                        for dKn in (-2,0):
                            Kn = Kd + dKn
                            if Kn < 0 : continue
                            if Kn > Jn : continue
                            Edn = En(Jn,Kn,An,Bn,Cn,Djn,Djkn,Dkn)
                            #KE = Edn+DBS-Ed
                            KE = Ed-(Edn+DBS)
                            #if KE <= 25 and KE > 0:
                            #    I *=1/5*KE**(1/2)
                            if I < 1e-5: continue
                            if KE <0: continue
                            #if dJn == -3:
                            #    I *=4.0
                            BE = hv-KE+ofs
                            linelist[(Jdd,Kdd,Jd,Jn,Kd,Kn,dJn,dKn)] = (BE,amp*I*Jamp)
                            #linelist[(Jd,Jn,Kd,Kn,Ed,Edn)] = (BE,amp*I)
                            spec += Jamp*amp*I*Gauss(BE,FWHM,Ev)
    return (spec,linelist,Ev)

def fit(Ev,ofs):            #Change amp to whatever parameter(s) you want to fit. Currently set for amp, T, FWHM, but easy to extend to others
    spec,linelist,Ev = spectrum(Ev,T,FWHM,amp,ofs,DBS,EA)
    return spec

#Curve Fit
#popt, pcov = curve_fit(fit,x,y,bounds=(-6,-2))
#print(popt)
#print(sqrt(diag(popt)))


spec,linelist,Ev = spectrum(x,T,FWHM,amp,ofs,DBS,EA)
linea = []
posa = []
#print('(Jdd,Jd,Kdd,Kd) (E,I)')
for line,posAmp in sorted(linelist.items(),key=lambda x: x[1][0]):
    print(line,posAmp)
    linea.append(line)
    posa.append(posAmp)

savetxt('linelist.dat',column_stack((linea,posa)))

#normalise
spec *= 1/max(spec)

#normalise fit
#ynnn = fit(x,*popt)
#normnnn = 1/max(ynnn)

plt.plot(Ev,spec,'C1',label='Autodetachment Model')
plt.plot(x,y,'C0',label='804.058nm PES')
#plt.plot(x,fit(x,*popt)*normnnn,'C3')
#plot params
#for i in range (0,size(popt)):
#    x0 = 12350
#    y0 = 0.588
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
plt.annotate(r'$^QN_2$',(12411,0.237),color='C0',fontsize=12,ha='right')
plt.annotate(r'$^QN_0$',(12428,0.986),color='C0',fontsize=12)
plt.annotate(r'$^oP_2$',(12416,0.484),color='C0',fontsize=12,ha='center')
#plt.annotate(r'$^QM_0$',(12411,0.309),color='C0',fontsize=12,ha='right')
plt.annotate(r'$^oO_2$',(12395,0.150),color='C0',fontsize=12,ha='center')
plt.annotate(r'$^ON_2$',(12375,0.150),color='C0',fontsize=12,ha='center')
plt.annotate(r'$^QN_0: \Delta K=0, \Delta J-3, K_0 \leftarrow K_0$',(12311,0.65))
plt.annotate(r'$^OP_2: \Delta K=-2, \Delta J-1, K_0 \leftarrow K_2$',(12311,0.60))
plt.annotate(r'$^QN_2: \Delta K=0, \Delta J+-3, K_2 \leftarrow K_2$',(12311,0.55))
#plt.title('autodetachment 804.058nm resonance')
plt.legend(loc=2)
plt.setp(ax.spines.values(), linewidth=1.5)

#Manual curvefit
#ofsa = arange(-2.8,-2.0,0.1)
#for i in range(0,size(ofsa)):
#    ynn=fit(x,ofsa[i])
#    ynn *= 1/max(ynn)
#    plt.plot(x,ynn,label=ofsa[i])

#plt.legend()
plt.savefig('Fig8a.pdf',dpi=400,bbbox_inches='tight')
plt.show()

#plt.plot(x,y)
#plt.plot(x,test(x,*popt))
#plt.show()
