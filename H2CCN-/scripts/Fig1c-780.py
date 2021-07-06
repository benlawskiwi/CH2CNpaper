from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

kb = 8.617333262145e-5
ev = 8065.74

width = 3.487*2
#height = width/1.618 #golden ratio
height = width/1.4
one_col = (width, height)
fig,ax = plt.subplots(figsize=one_col)


x,y = loadtxt('../data/CH2CN-4F1CMI1536PESBE.dat',unpack=True)
#x,y = loadtxt('CH2CN-test.dat',unpack=True)

VMI = loadtxt('../data/CH2CN-6N3CMb512.txt')
rows, cols = VMI.shape
trim = 46
#r2 = rows//2 - trim
#c2 = cols//2 - trim
r2 = rows//2
c2 = cols//2


VMI = VMI[trim:-trim, trim:-trim]


#x,y = loadtxt('ANUscan.dat',delimiter=',',unpack=True)
#x = 1e7/x
#Remove background
#step = 10
#n = round(size(x)/step)
#xx = []
#yy = []
#for i in range(0,n-1):
#    s = slice(step*i,step*i+step,1)
#    a = mean(x[s])
#    b = min(y[s])
#    xx.append(a)
#    yy.append(b)
#bkg = mean(yy)
#y -= bkg
#print(bkg)

#Fitting Parameters
#T = 150
T = 199
#T = 2.7
#amp = 1
amp = 0.0008816272559812547
EA = 12468
#EA = 12481.3
DBS = 39
#FWHM = 2 
FWHM = 6.1
#ofs = -1.073
#ofs = +8.2
ofs = 13.3
ofs = 1.2

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

#Detachment Threshold
D0 = 129      #Neumark


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

def spectrum(Ev,T,FWHM,amp,ofs):
    #Set up spectrum
    #Ev = arange(EA-200,EA+200,0.1)
    spec = zeros(size(Ev))
    linelist = {}
    #nu = EA-DBS+ofs
    nu = EA
    #Selection rules/intensities                    (Jdd,Kdd) -> (Jd,Kd)
    for Jdd in range(0,32):
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
                    #if Ed < D0: continue  #These lines won't autodetach
                    dE = Ed-Edd+nu+ofs
                    linelist[(Jdd,Jd,Kdd,Kd,dJ,dK)] = (dE,amp*I)
                    spec += amp*I*Gauss(dE,FWHM,Ev)
    return (spec,linelist,Ev)

def fit(Ev,T):            #Change amp to whatever parameter(s) you want to fit. Currently set for amp, T, FWHM, but easy to extend to others
    spec,linelist,Ev = spectrum(Ev,T,FWHM,amp,ofs)
    return spec

#Curve Fit
#popt, pcov = curve_fit(fit,x,y)
#print(popt)


spec,linelist,Ev = spectrum(x,T,FWHM,amp,ofs)
linea = []
posa = []
#print('(Jdd,Jd,Kdd,Kd) (E,I)')
for line,posAmp in sorted(linelist.items(),key=lambda x: x[1][0]):
    print(line,posAmp)
    linea.append(line)
    posa.append(posAmp)

#savetxt('linelistFull.dat',column_stack((linea,posa)))

#Peak Labels
#pts = [12351,12368,12385,12402,12420,12438,12456,12475,12494,12514]
#pi = [5,4,3,2,1,0,1,2,3,4]
#pf = [4,3,2,1,0,1,2,3,4,5]
#ypt = fit(x,*popt)
#for i in range (0, size(pts)):
#    subr = logical_and(x>pts[i]-5,x<pts[i]+5)
#    yr = ypt[subr]
#    xr = x[subr]
#    ymax = max(yr)
#    xmax = xr[argmax(yr)]
#    lab = r"K''"+str(pi[i])+r"$\rightarrow$ K'"+str(pf[i])
#    plt.annotate(lab,(xmax,ymax+10),ha='center')

plt.plot(x,y,'C0',label='781.15nm PES')
plt.plot(Ev,spec,'C1',label='Rotational Model')
#plt.plot(x,y,'C0',label='781.15nm PES')
#plt.plot(x,fit(x,*popt),'C3')
#plot params
#for i in range (0,size(popt)):
#    x0 = 12350
#    y0 = 500-50*i
#    perr = sqrt(diag(pcov[i]))
#    perrs=round(perr[0][0],3)
#    st = 'param = '+str(round(popt[i],3))+' +/- '+str(perrs)
#    plt.annotate(st,(x0,y0))

#savetxt('linelist.dat',column_stack((x,fit(x,*popt))))

#print('this is your offset')
#print(popt[0])

plt.setp(ax.spines.values(), linewidth=1.5)
plt.xticks([12300,12400,12500,12600,12700],fontsize=12)
plt.yticks([])
ax.tick_params(width=1.5)

fig.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.98)
plt.xlabel(r'electron binding energy (cm$^{-1}$)',fontsize=14)
plt.ylabel('intensity (arb. u.)',fontsize=14)
plt.xlim(12300,12700)
plt.annotate("",
            xy=(12468, 0.84), xycoords='data',
            xytext=(12468, 1.0), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"),
            )
#plt.legend()
plt.minorticks_on()

sub_axes=plt.axes([.70,.70,.25,.25])
sub_axes.imshow(VMI, extent=[-r2, r2, -c2, c2], vmin=0, cmap='jet', vmax=VMI.max()/1)
sub_axes.set_xticks([0, 256])
sub_axes.set_xticklabels(["0", "256"], fontsize="smaller")
sub_axes.set_yticks([])
sub_axes.set_yticklabels([])


plt.savefig('Fig1c-781.pdf',dpi=400,bbbox_inches='tight')
plt.show()

#plt.plot(x,y)
#plt.plot(x,test(x,*popt))
#plt.show()
