import os
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
from pwtools import thermo, parse, crys, constants, num, io, eos
import sympy as sp
import pandas as pd

# the program assumes orthorhombic symmetry!

# specify parameters

Tmin=5
Tmax=700
NT=279

T = np.linspace(Tmin, Tmax, NT) #K
#T = np.array([298])
P = np.array([1e-4]) #GPa

volDirs = np.array([0.955, 0.965, 0.975, 0.985, 0.995, 1.005, 1.015, 1.025, 1.035, 1.045, 1.055, 1.065, 1.075, 1.085, 1.095, 1.105, 1.125, 1.135, 1.145])
# create a dictionary of lists
data = {'volume factor': [], 'electronic energy': [], 'phononic DOS': [], 'cell constants': []}

for volDir in volDirs:
    # go to the directory of the respective volume
    os.chdir(str(volDir))
    # load the density of states
    dos = np.loadtxt('matdyn2.dos') # wave numbers in cm^-1, DOS in cm
    # add it to the list of [Freq, DOS]-arrays
    data['phononic DOS'].append(dos)
    # parse the pw.x output, get the electronic energy in Ry and convert it to eV
    structure = io.read_pw_scf('qe-pw.out')
    etot = structure.get_etot() # convert Rydberg to eV
    data['electronic energy'].append(etot)
    # get the cell constants in Bohr and convert them to Angstrom
    cell = structure.get_cell() # convert Bohr to Angstrom

    data['cell constants'].append(np.array([cell[0][0], cell[1][1], cell[2][2]]))
    os.chdir('..')

# convert cell constants and DOS into numpy arrays
data['cell constants'] = np.array(data['cell constants'])
data['phononic DOS'] = np.array(data['phononic DOS'])

# definde the volume function for the orthorhombic unit cell
volfunc_ax = lambda x: crys.volume_cc(np.array([x[0], x[1], x[2]] + [90]*3))

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

gibbs.set_fitfunc('C', lambda x,y: num.PolyFit1D(x, y, deg=10)) # Fit Gmin-T
gibbs.set_fitfunc('1d-G', lambda x,y: eos.EosFit(x,y)) # Fit G-V
gibbs.set_fitfunc('alpha', lambda x,y: num.PolyFit1D(x,y, deg=7)) # Fit G-V

g = gibbs.calc_G(calc_all=True)

def Cp(T,a,b,c,d,e,f,h):
    cp = a + b * T + c * T ** (-2) + d * T ** 2 + e * T ** 3 + f / T + h * T * np.log(T)
    return cp

popt,pconv = scp.optimize.curve_fit(Cp, np.array(g['/T/T']), np.array(g['/#opt/T/P/Cp'][:, 0])*constants.R/2)

print(popt)



plt.figure(1)
plt.scatter(g['/ax0-ax1-ax2/V'], g['/ax0-ax1-ax2/Etot'])
plt.xlabel('Volumen [A^3]')
plt.ylabel('Etot in eV')
plt.savefig('parse-polyfit/Etot-V.png')
plt.close(1)

# plot E-V for the first T
plt.figure(2)
plt.scatter(g['/ax0-ax1-ax2/V'], g['/ax0-ax1-ax2/T/Evib'][:,9], s=1)
plt.xlabel('Zellvolumen in Angstrom^3')
plt.ylabel('Evib in eV')
plt.savefig('parse-polyfit/Evib-V.png')
plt.close(2)

plt.figure(3)
plt.plot(g['/T/T'], np.array(g['/#opt/T/P/Cp'][:,0])*constants.R/2, linewidth=1, label='derivative of G')
plt.plot(g['/T/T'], Cp(g['/T/T'], *popt), linewidth=1, label='Cp-fit')
plt.xlabel(r'$T ~ [K]$')
plt.ylabel(r'$C_p ~ [J ~K^{-1} ~mol^{-1}]$')
plt.legend()
plt.savefig('parse-polyfit/Cp-T.png')
plt.close(3)

plt.figure(4)
plt.scatter(g['/ax0-ax1-ax2/V'], g['/T/P/ax0-ax1-ax2/G'][49,0])
plt.xlabel('V [Ang^3]')
plt.ylabel('G in eV')
plt.savefig('parse-polyfit/G-V.png')

plt.figure(5)
plt.scatter(g['/T/T'], g['/#opt/T/P/V'][:,0], )
plt.xlabel('T [K]')
plt.ylabel('V in Ang^3')
plt.savefig('parse-polyfit/V-T.png')

plt.figure(6, figsize=(7,5))
plt.scatter(g['/T/T'], g['/#opt/T/P/G'][:,0])
plt.xlabel('T [K]')
plt.ylabel('G in eV')
plt.savefig('parse-polyfit/G-T.png')
plt.close(6)


'''
pd.DataFrame(np.transpose(np.array([g['/T/T'], g['/#opt/T/P/V'][:,0]]))).to_csv('V-T.csv')

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

print(Cpopt)

plt.figure(8, figsize=(7, 4.8))
plt.scatter(g['/T/T'], np.array(g['/#opt/T/P/G'][:,0])/constants.J_to_eV*constants.avo/2, s=1.5)
plt.xlabel(r'$T ~ [K]$')
plt.ylabel(r'$G ~ [J ~mol^{-1}]$')
plt.savefig('G-T.png', dpi=500)

plt.figure(10)
plt.scatter(T, npS, s=1.5)
plt.xlabel(r'$T ~ [K]$')
plt.ylabel(r'$S ~ [J ~K^{-1} ~mol^{-1}]$')
plt.savefig('S-T.png', dpi=500)

plt.figure(7)
plt.scatter(T, npCp, s=1)
plt.plot(T, Cp(T, *Cpopt), linewidth=1, label='Cp-fit')
plt.xlabel(r'$T ~ [K]$')
plt.ylabel(r'$C_p ~ [J ~K^{-1} ~mol^{-1}]$')
plt.legend()
plt.savefig('Cp-T.png', dpi=500)
'''



Cpdata = {'T [K]': g['/T/T'], 'G [J/mol]': np.array(g['/#opt/T/P/G'][:,0])/constants.J_to_eV*constants.avo/2}

pdata=pd.DataFrame(Cpdata)

pdata.to_csv('parse-polyfit/Ba(BaH4)2-Cp-polyfit.csv')

print('done')
