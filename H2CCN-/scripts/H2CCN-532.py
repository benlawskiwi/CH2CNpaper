import numpy as np
import matplotlib.pyplot as plt
import os
import glob 
from scipy.interpolate import UnivariateSpline
from lxml import etree
from lmfit import minimize, Parameters, report_fit, fit_report

def Gd(pars, mode, vib):
    omegae = pars['omegae'+str(mode)].value
    if mode == 5:
        omegaexe = pars['omegaexe'+str(mode)].value
        return omegae*(vib+1/2) - omegaexe*(vib+1/2)**2
    else:
        return omegae*(vib+1/2)

def Gdd(pars, mode, vib):
    return pars['omegaedd'+str(mode)].value * (vib+1/2)


def gauss(x, x0, FWHM):
    sigma = FWHM/2/np.sqrt(np.log(2))
    return np.exp(-((x-x0)/sigma)**2)/sigma/np.sqrt(np.pi)


def pgo_spectrum(FWHM, eBEb, PESb, epgo, pgo):
    spect = np.zeros_like(PESb)

    subrpgo = np.logical_and(epgo >= eBEb[0], epgo <= eBEb[-1])

    # Gaussian function for pgo each transition
    for e, amp in zip(epgo[subrpgo], pgo[subrpgo]):
        subr = np.logical_and(eBEb > e-2*FWHM, eBEb < e + 2*FWHM)
        spect[subr] += gauss(eBEb[subr], e, FWHM)*amp

    return spect


# xml parameter substitution and least-squares functions -----------
def setparam(param, value, tree):
    # set xml parameter, multiple parameter has int at end
    occurrence = 1
    if param[-1].isdigit():
        occurrence = int(param[-1]) + 1
        param = param[:-1]

    seen = 0
    for x in root.iter('Parameter'):
        if x.get('Name') == param:
            seen += 1
            if occurrence == seen:
                x.set('Value', str(value))
                return

def evaluate_pgo(pars):
    global pC0, pC1, pTemperature

    print('evaluate pgo')
    setparam('Gaussian', 0, tree)
    setparam('Origin', 0, tree)
    setparam('Fmin', -250, tree)
    setparam('Fmax', 250, tree)
    setparam('A0', pars['A0'].value, tree)
    setparam('B0', pars['B0'].value, tree)
    setparam('C0', pars['C0'].value, tree)
    setparam('A1', pars['A1'].value, tree)
    setparam('B1', pars['B1'].value, tree)
    setparam('C1', pars['C1'].value, tree)
    setparam('Temperature', pars['Temperature'].value, tree)
    tree.write('ffit.pgo', xml_declaration=True, encoding="utf-8")

    os.system('/usr/local/bin/pgo --plot ffit.pgo ffit.dat')
    pC0 = pars['C0'].value
    pC1 = pars['C1'].value
    pTemperature = pars['Temperature'].value

    epgo, pgo = np.loadtxt("ffit.dat", unpack=True)
    peaks = pgo > 0
    epgo = epgo[peaks]
    pgo = pgo[peaks]

    return epgo, pgo


def residual(pars, eBE, PES, tree, plot=False):
    global epgo, pgo, pC0, pC1, pTemperature

    model = np.zeros_like(eBE)

    profile = {}
    # for each mode evaluate the pgopher profile
    for band, pkeBE in transitions.items():
        combination = band.split('_')
        
        nu = 0 
        for b in combination:
            mode = int(b[:-2]); vd = int(b[-2]); vdd = int(b[-1])

            Ed = Gd(pars, mode, vd) - Gd(pars, mode, 0)
            Edd = Gdd(pars, mode, vdd) - Gdd(pars, mode, 0)
            nu += Ed - Edd

        # mode eBE
        eBEn = nu + pars['offset'].value
        subr = np.logical_and(eBE>eBEn-200, eBE<eBEn+200)

        FWHM = pars['FWHM'].value
        FWHMreduce = pars['slope'].value*(eBEn-EA)/1000
        if FWHMreduce > 0:
            FWHM -= FWHMreduce
            if FWHM <= 0:
                FWHM = 0.1
        amp = pars['amp'+band].value

        #if pC0 != pars['C0'].value or pC1 != pars['C1'].value or\
        #   pTemperature != pars['Temperature'].value or pC0 == 0.0:
        #    epgo, pgo = evaluate_pgo(pars)
          
        #spectrum = pgo_spectrum(FWHM, eBE[subr], PES[subr], epgo+eBEn, pgo)*amp
        # stg - 4/6/19 - use Gaussian function for the whole band 
        spectrum = gauss(eBE[subr], eBEn, FWHM)*amp

        model[subr] += spectrum
        
        if plot:
            profile[band] = (eBE[subr], spectrum)

    if plot:
        return model, profile
    else:
        return model - PES

def tolatex(bands):
    band = bands.split('_')
    mainbr = True
    if len(band) > 1 or band[-1] == '1':
        mainbr = False
    lbl = ''
    for b in band:
       lbl += rf'${b[:-2]}^{b[-2]}_{b[-1]}$' 
       if mainbr and b[-1] == '1': mainbr = False
    return lbl, mainbr

# main --------------------------
tree = etree.parse('ch2cn-fit.pgo')
root = tree.getroot()

anion = os.getcwd().split('/')[-1]

files = sorted(glob.glob(f'*nm/*/{anion}*PES_qA.dat'))
if len(files) > 1:
    print()
    j = []
    for i, fn in enumerate(files):
        deblobbed = 'b' in fn
        if not deblobbed and fn.replace('PES', 'bPES') in files:
            continue    # skip if deblobed version exists
        j.append(i)
        indent = ''
        print(f'  {indent}[{len(j)-1}] {fn}')

    num = input('\nVMI file number: ')
    fn = files[j[int(num)]]
else:
    fn = files[0]
print(f'  {fn}\n')

eBE, PES = np.loadtxt(fn, unpack=True)

EA = 12468  # cm-1
EAint = int(EA)

# drop data with no information
subE = eBE > 11800
eBE = eBE[subE]
PES = PES[subE]

# Neumark table II - eBE position is not used
transitions = {
 '901': 12000, '621': 12000, '611': 12000,
 '501': 12000,
 '000': 12468, '610': 12888, '511': 12997,
 '510': 13137, '920': 13200, '620': 13279,
 '410': 13495, '630': 13000, '520': 13808,
 '310': 13907, '520_610': 14224, '531': 14386,
 '530': 14483, '310_510': 14600, '540': 15171,
 '310_520': 15291, '320': 15360, '550': 15500
}

# drop transitions beyond threshold
transitions = {k:v for k, v in transitions.items()\
               if v >= eBE[0] and v <= eBE[-1]+20}


# initial values 
pars = Parameters()

for b in transitions.keys():
    dummy, mainbr = tolatex(b)
    print(b)
    pars.add('amp'+b, value = 20 if mainbr else 1, vary=True, min=0,
             max = 200 if mainbr else 100)

strFWHM = input('\nFWHM  < 0 vary [5]: ')
if len(strFWHM) > 1:
    FWHM = abs(float(strFWHM))    
else:
    FWHM = 5

pars.add('FWHM', value=FWHM, min=0.1, max=200, vary=True)
pars.add('slope', value=4, vary=False, min=0, max=50)

pars.add('offset', value = EA, min=EA-5, max=EA+5, vary=True)
pars.add('bkg', value = 0.004, vary=False)
pars.add('omegae0', value = 0.0, vary=False)
pars.add('omegae3', value = 1453, min=1300, max=1500, vary=False)
pars.add('omegae4', value = 1023, min=1000, max=1100, vary=False)
pars.add('omegae5', value = 634+36, min=600, max=700, vary=True)
pars.add('omegae6', value = 415, min=390, max=430, vary=False)
pars.add('omegae9', value = 361, min=300, max=400, vary=False)

pars.add('omegaexe5', value = -1, min=-30, max=30, vary=False)

pars.add('omegaedd0', value = 0, vary=False)
pars.add('omegaedd3', value = 1419, min=1400, max=1440, vary=False)
pars.add('omegaedd4', value = 1061, min=1040, max=1080, vary=False)
pars.add('omegaedd5', value = 126, min=50, max=150, vary=False)
pars.add('omegaedd6', value = 590, min=523, max=700, vary=False)
pars.add('omegaedd9', value = 418, min=400, max=440, vary=False)

pars.add('A0', value = 9.29431, min=9, max=9.6, vary=False)
pars.add('B0', value = 0.338427, min=0.31, max=0.36, vary=False)
pars.add('C0', value = 0.32761, min=0.30, max=0.35, vary=False)
pars.add('A5', value = 8, min=5, max=15, vary=True) # v5 umbrellla

pars.add('A1', value = 9.506, min=9.2, max=9.8, vary=False)
pars.add('B1', value = 0.347799, min=0.31, max=0.37, vary=False)
pars.add('C1', value = 0.329429, min=0.3, max=0.6, vary=False)

pars.add('Temperature', value = 237-50, min=150, max=350, vary=False)

# keep track: only call pgopher if parameters C0, C1, FWHM change
pC0 = 0.0
pC1 = 0.0
pTemperature = 0.0
epgo = None
pgo = None

# least squares fit --------------------------------
result = minimize(residual, pars, args=(eBE, PES, tree))

report = fit_report(result)  # report_fit(result.params)
print(report)

reportfn = fn[:-10] + 'report.txt'
np.savetxt(reportfn, report.split('\n'), fmt='%s')
print(f'fit report saved to {reportfn}')

model, profile = residual(result.params, eBE, PES, tree, plot=True)
fnmodel = fn[:-10]+'model.dat'
np.savetxt(fnmodel, np.column_stack((eBE, model)))
print(f'model fit saved to {fnmodel}')

# plot -----------
plt.plot(eBE, PES, '-', color='k', ms=1, label='PES')
plt.plot(eBE, model, '--', color='C2', label='Model')


for i, (pband, (peBE, pgauss)) in enumerate(profile.items()):
    if len(pgauss) == 0:
        continue 

    pbandstr, mainbr = tolatex(pband)
    if len(pband) == 3:
        col = int(pband[0]) % 10
    else:
        col = int(pband[:2]) % 10

    #plt.plot(peBE, pgauss, f'-C{col:d}', lw=1, ls='-' if mainbr else '--')

    pkeBE = peBE[pgauss.argmax()]
    pk = pgauss.max()
    ofs = 0.05 + mainbr*0.1 if pband != '000' else 0.02

    #plt.plot((pkeBE, pkeBE), (pk+0.02, pk+ofs), f'--C{col}', lw=1)
    #plt.annotate(pbandstr, (pkeBE, pk+ofs+0.02),
    #             ha='center', color=f'C{col}' if mainbr else 'w',
    #             fontsize='large' if mainbr else 'small') 
    print('band - eBE - pk')
    print(str(pbandstr)+' - '+str(pkeBE)+' - '+str(pk+ofs))

    np.savetxt(fnmodel+pband, np.column_stack((peBE, pgauss)))

plt.annotate(r'$9^0_1$',(12054,0.074), ha='center', color='C2', fontsize='large')
plt.annotate(r'$0^0_0$',(12616,1.000), ha='center', color='C2', fontsize='large')
plt.annotate(r'$6^2_1$',(12710,0.111), ha='center', color='C2', fontsize='large')
plt.annotate(r'$5^1_1$',(13014,0.368), ha='center', color='C2', fontsize='large')
plt.annotate(r'$5^1_0$',(13140,0.272), ha='center', color='C2', fontsize='large')
plt.annotate(r'$5^2_0$',(13814,0.303), ha='center', color='C2', fontsize='large')
plt.annotate(r'$5^3_0$',(14492,0.164), ha='center', color='C2', fontsize='large')
plt.annotate(r'$5^3_1$',(14364,0.164), ha='center', color='C2', fontsize='large')
plt.annotate(r'$5^4_0$',(15170,0.099), ha='center', color='C2', fontsize='large')
plt.annotate(r'$5^5_0$',(15850,0.065), ha='center', color='C2', fontsize='large')
plt.annotate(r'$3^1_0$',(13937,0.213), ha='center', color='C2', fontsize='large')
plt.annotate(r'$3^2_0$',(15375,0.087), ha='center', color='C2', fontsize='large')
plt.annotate(r'$4^1_0$',(13494,0.087), ha='center', color='C2', fontsize='large')

plt.xlim(11800,16200)
plt.legend()
plt.xlabel(r'eBE (cm$^{-1}$)')
plt.ylabel('intensity (arb. u.)')
#plt.title(fn, fontsize='small')
plt.savefig(f'H2CCN-532.pdf', bb_inches='tight')

plt.show()
