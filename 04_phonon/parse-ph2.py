import os
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
from pwtools import thermo, parse, crys, constants, num, io, eos
import sympy as sp
import pandas as pd

# the program assumes orthorhombic symmetry!

# specify parameters

Tmin=1
T1=5
Tmax=500
NT1=11
NT2=100

bar = 6.241509125883258e-07

T = np.unique(np.concatenate((np.linspace(Tmin, T1, NT1),np.array([298.15]),np.linspace(T1,Tmax,NT2),)))#K
print('T [K]:', T)
P = np.array([1e-4]) #GPa
print('P [eV/Ang^3]:', P*100*bar)
volDirs = np.array([0.955, 0.965, 0.975, 0.985, 0.995, 1.005, 1.015, 1.025, 1.035, 1.045, 1.055, 1.065, 1.075, 1.085, 1.095, 1.105, 1.125, 1.135, 1.145])
print('volume scale factors:', volDirs)
# define a volume function of the cell constants
volfunc_ax = lambda x: crys.volume_cc(np.array([x[0], x[1], x[2]] + [90]*3))
# create a dictionary of lists
data = {'volume factor': [], 'electronic energy': [], 'phononic DOS': [], 'cell constants': []}


for volDir in volDirs:
    print('Iteration for volume factor ',volDir)
    # go to the directory of the respective volume
    os.chdir(str(volDir))
    # load the density of states
    dos = np.loadtxt('matdyn2.dos', skiprows=0, usecols=(0,1)) # wave numbers in cm^-1, DOS in cm
    # add it to the list of [Freq, DOS]-arrays
    data['phononic DOS'].append(dos)
    # parse the pw.x output, get the electronic energy in Ry and convert it to eV
    structure = io.read_pw_scf('qe-pw.out')
    etot = structure.get_etot() # convert Rydberg to eV
    data['electronic energy'].append(etot)
    # get the cell constants in Bohr and convert them to Angstrom
    cell = structure.get_cell() # convert Bohr to Angstrom
    data['cell constants'].append(np.array([cell[0][0], cell[1][1], cell[2][2]]))
    # do statistical mechanics
    os.chdir('..')

# convert cell constants and DOS into numpy arrays
data['cell constants'] = np.array(data['cell constants'])
data['phononic DOS'] = np.array(data['phononic DOS'])

# definde the volume function for the orthorhombic unit cell

print('cell constants:',data['cell constants'])

gibbs = thermo.Gibbs(T=T,
                     P=P,
                     temp=None,
                     etot=data['electronic energy'],
                     phdos=data['phononic DOS'],
                     axes_flat=data['cell constants'],
                     volfunc_ax=volfunc_ax,
                     case='1d',
                     skipfreq=True,
                     dosarea=66,
                     integrator=scp.integrate.simps,
                     verbose=True,
                     eps=np.float64(3.3306690738754696e-16),
                     fixnan=True, nanfill=0.0
                    )

f = gibbs.calc_F(calc_all=True)

print('T [K]; F [eV]:')
print(f['/ax0-ax1-ax2/T/F'][0,:])

plt.figure(1)
plt.scatter(f['/ax0-ax1-ax2/V'], f['/ax0-ax1-ax2/Etot'])
plt.xlabel('Volumen [A^3]')
plt.ylabel('Etot in eV')
plt.savefig('parse2/Etot-V.png')

plt.figure(2)
plt.scatter(f['/T/T'], f['/ax0-ax1-ax2/T/F'][0,:])
plt.xlabel('Volumen [A^3]')
plt.ylabel('F in eV')
plt.savefig('parse2/F-V.png')

"""
plt.figure(1)
plt.scatter(g['/ax0-ax1-ax2/V'], g['/ax0-ax1-ax2/Etot'])
plt.xlabel('Volumen [A^3]')
plt.ylabel('Etot in eV')
plt.savefig('Etot-V.png')

# plot E-V for the first T
plt.figure(2)
plt.scatter(g['/ax0-ax1-ax2/V'], g['/ax0-ax1-ax2/T/Evib'][:,9], s=1)
plt.xlabel('Zellvolumen in Angstrom^3')
plt.ylabel('Evib in eV')
plt.savefig('Evib-V.png')

plt.figure(4)
plt.scatter(g['/ax0-ax1-ax2/V'], g['/T/P/ax0-ax1-ax2/G'][49,0])
plt.xlabel('V [Ang^3]')
plt.ylabel('G in eV')
plt.savefig('G-V.png')

plt.figure(5)
plt.scatter(g['/T/T'], g['/#opt/T/P/V'][:,0])
plt.xlabel('T [K]')
plt.ylabel('V in Ang^3')
plt.savefig('V-T.png')

#io.write_h5("Ba(BH4)2_gibbs.h5", g)

#scipy section

# define Cp-Function
T,a,b,c,d,e,f,h,C1,C2 = sp.symbols('T a b c d e f h C1 C2')

spCp = a + b*T + c*T**(-2) + d*T**2 + e*T**3 + f/T + h*T*sp.log(T)

spG = sp.integrate(-spCp/T,T,T)+ C1*T + C2

#S = sp.integrate(Cp, T) + C1 ????

print(spG)

lamCp = sp.lambdify([T,a,b,c,d,e,f,h], spCp)
lamG = sp.lambdify([T,a,b,c,d,e,f,h,C1,C2], spG)

T = np.linspace(Tmin, Tmax, NT) #K

def Cp(T,a,b,c,d,e,f,h):
    return lamCp(T,a,b,c,d,e,f,h)

def G(T,a,b,c,d,e,f,h,C1,C2):
    return lamG(T,a,b,c,d,e,f,h,C1,C2)


print(np.array(g['/#opt/T/P/G'][:, 0])/constants.J_to_eV*constants.avo/2)
popt,pconv = scp.optimize.curve_fit(G, np.array(g['/T/T']), np.array(g['/#opt/T/P/G'][:, 0])/constants.J_to_eV*constants.avo/2)



scpG = scp.interpolate.CubicSpline(T, np.array(g['/#opt/T/P/G'][:, 0])/constants.J_to_eV*constants.avo/2,bc_type='not-a-knot',extrapolate=False)
scpCp = -T*scpG(T,2)

npS = (-1)*np.gradient(np.array(g['/#opt/T/P/G'][:, 0])/constants.J_to_eV*constants.avo/2,T, edge_order=2)

npCp = T*np.gradient(npS,T, edge_order=2)

Cpopt,Cpconv = scp.optimize.curve_fit(Cp, T, npCp)

print(popt)
print(Cpopt)

plt.figure(8, figsize=(7, 4.8))
#plt.plot(g['/T/T'], G(g['/T/T'], *popt), linewidth=0.7, label='Gfit')
plt.scatter(g['/T/T'], np.array(g['/#opt/T/P/G'][:,0])/constants.J_to_eV*constants.avo/2, s=1.5, label='data')
#plt.plot(g['/T/T'], scpG(g['/T/T']), linewidth=0.7, label='cubic spline')
plt.xlabel('T [K]')
plt.ylabel('G [J/mol]')
plt.legend()
plt.savefig('G-T.png', dpi=500)

plt.figure(10)
plt.scatter(T, npS, s=1.5, label='numerical derivative')
plt.xlabel('T [K]')
plt.ylabel('S  [J/K/mol]')
plt.legend()
plt.savefig('S-T.png', dpi=500)

plt.figure(7)
#plt.plot(g['/T/T'], Cp(g['/T/T'],popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6]), linewidth=1, label='Gfit')
#plt.scatter(g['/T/T'], np.array(g['/#opt/T/P/Cp'][:,0])*constants.R/2, s=1.5, label='poly')
plt.scatter(T, npCp, s=1, label='numerical derivative')
plt.plot(T, Cp(T, *Cpopt), linewidth=1, label='Cp-fit')
#plt.plot(g['/T/T'], scpCp, linewidth=1, label='cubic spline')
plt.xlabel('T [K]')
plt.ylabel('Cp [J/K/mol]')
plt.legend()
#plt.ylim(bottom=0)
plt.savefig('Cp-T.png', dpi=500)

Cpdata = {'T [K]': T, 'G [J/mol]': np.array(g['/#opt/T/P/G'][:,0])/constants.J_to_eV*constants.avo/2, 'S [J/mol/K]': npS, 'Cp [J/mol/K]': npCp}

pdata=pd.DataFrame(Cpdata)

pdata.to_csv('Ba(BaH4)2-Cp.csv')
"""
print('done')
