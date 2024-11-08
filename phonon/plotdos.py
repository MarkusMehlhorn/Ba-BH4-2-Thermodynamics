import os

import numpy as np
from pwtools import thermo
import scipy
from matplotlib import pyplot as plt

volDirs = np.array([0.955, 0.965, 0.975, 0.985, 0.995, 1.005, 1.015, 1.025, 1.035, 1.045, 1.055, 1.065, 1.075, 1.085, 1.095, 1.105, 1.125, 1.135, 1.145])
#print(VolRange)
cm=0.3937

freqs = []
doss= []
for volDir in volDirs:

    os.chdir(str(volDir))
    dospack = np.loadtxt('matdyn2.dos')
    dospackt=np.transpose(dospack)
    #print(dospackt)
    freq=dospackt[0]
    #print(freq)
    dos=dospackt[1]
    #print(dos)
    freqs.append(freq)
    doss.append(dos)

    Int=scipy.integrate.trapz(dos,freq)
    print(Int)
    T=np.linspace(0,500, 101)
    HT = thermo.HarmonicThermo(freq,
                               dos,
                               T=T,
                               temp=None,
                               skipfreq=True,
                               eps=3.3306690738754696e-16,
                               fixnan=True, nanfill=0.0,
                               dosarea=66,
                               integrator=scipy.integrate.simps,
                               verbose=True)
    Cv= HT.cv()

    #print(Cv)
    plt.figure(str(volDir)+'Cv')
    plt.scatter(T, Cv, s=0.1)
    plt.savefig('./plots/Cv.png', dpi=400)
    plt.close(str(volDir)+'Cv')


    fig = plt.figure(str(volDir)+'dos')
    fig, ax = plt.subplots(figsize=(12*cm,6*cm))
    fig.dpi=500
    ax.tick_params(labelsize=7)
    ax.set_xlabel(r'$\lambda$ / cm$^{-1}$', fontsize=9)
    ax.set_xticks(np.arange(0, 2501, 250))
    ax.set_xlim(0, 2500)
    ax.set_ylabel('DOS / cm', fontsize=9)
    ax.set_yticks(np.arange(0, 2.51, 0.25))
    ax.set_ylim(-0.05, 2.0)
    ax.plot(freq, dos, linewidth=0.7)

    plt.savefig('./plots/dos.png', dpi=1000, bbox_inches="tight")
    plt.close(str(volDir)+'dos')

    os.chdir('..')




plt.figure('dos')
delta_Y = 0.2
for i in range(len(doss)):
  y_change = delta_Y * (i+1)
  doss[i] = doss[i] - y_change
  plt.plot(freqs[i], doss[i], label=volDirs[i])
#plt.legend()
plt.savefig('dos.png', dpi=400)

