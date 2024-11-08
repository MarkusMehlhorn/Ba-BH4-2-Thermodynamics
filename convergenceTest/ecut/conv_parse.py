import os
import numpy as np
from ase.units import Ry, eV
from matplotlib import pyplot as plt
import math


parameters = [ f.name for f in os.scandir() if f.is_dir() ]
parameters.remove('.idea')
parameters.remove('template')

parameterrange=[]

for parameter in parameters:
    parameterrange.append(float(parameter))

parameterrange.sort()

E=[]

for parameter in parameterrange:
    os.chdir(str(parameter))
    with open("qe-pw.out", "r") as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if "!" in lines[i]:
            E.append(float(lines[i].split()[4]))
    os.chdir('..')

DeltaE=[]

for i in range(len(parameterrange)):
    if i==0:
        DeltaE.append(-1)
    else:
        DeltaE.append(E[i]-E[i-1])

DeltaE_meVperat = eV*1e3/Ry*np.array(DeltaE)/22

print(parameterrange)
print(E)
print(DeltaE)
print(DeltaE_meVperat)
print(Ry)
print(eV)

lgDeltaE_meVperat=np.log10(-1*DeltaE_meVperat)

table = np.array([parameterrange, DeltaE, DeltaE_meVperat, lgDeltaE_meVperat])
plt.scatter(parameterrange, lgDeltaE_meVperat)
plt.ylabel('$lg(\Delta E \quad [meV/at]^{-1})$')
plt.xlabel('$E_{cut}$ [Ry]')
plt.xlim(left=70)
plt.savefig('conv_ecut.png')
np.savetxt('conv_parameter.txt', np.transpose(table))
