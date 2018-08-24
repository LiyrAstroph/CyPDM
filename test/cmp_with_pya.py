import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming import pyPDM
import sys;
sys.path.append("../")
import cypdm 

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

nbins = 10
covers = 3

cres1 = np.loadtxt("cresult1.txt")
cres2 = np.loadtxt("cresult2.txt")

con = np.loadtxt("con_all.txt")
t0 = con[0, 0]
con[:, 0] -= t0

S = pyPDM.Scanner(minVal=1.0/(50.0*365.0), maxVal=1.0e-2, dVal=1.0e-5, mode="frequency")
P = pyPDM.PyPDM(con[:, 0], con[:, 1])
f1, t1 = P.pdmEquiBinCover(nbins, covers, S)
f2, t2 = P.pdmEquiBin(nbins, S)

pdm = cypdm.PyPDM(con[:, 0], con[:, 1], nbins, covers)
myp = 1.0/f1
myt1 = pdm.getPDM_EquiBinCover(myp)
myt2 = pdm.getPDM_EquiBin(myp)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(211)
plt.plot(1.0/f1, t1, color='k', label="PyAstronomy")
plt.plot(myp, myt1,color='g', marker='o', markersize=2, label="CyPDM")
plt.plot(cres1[:, 0], cres1[:, 1], color='b', label="CPDM")

ax.legend()
ax.set_xlabel(r"Period")
ax.set_ylabel(r"$\Theta_{\rm PDM}$")

ax = fig.add_subplot(212)
plt.plot(1.0/f2, t2, color='k', label="PyAstronomy")
plt.plot(myp, myt2,color='g', marker='o', markersize=2, label="CyPDM")
plt.plot(cres2[:, 0], cres2[:, 1], color='b', label="CPDM")

ax.legend()

ax.set_xlabel(r"Period")
ax.set_ylabel(r"$\Theta_{\rm PDM}$")

plt.show()
