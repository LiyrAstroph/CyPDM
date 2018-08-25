import sys
sys.path.append("../")
import cypdm 
import numpy as np 
#import matplotlib.pyplot as plt 

nbins=10
covers=3

nd = 1000
P = 10.0
jd = np.linspace(0.0, 100.0, nd, dtype=np.double)
fs = np.sin(jd/P * 2.0*np.pi) + np.random.randn(nd)*0.1

pdm = cypdm.PyPDM(jd, fs, nbins, covers)
myp = 1.0/np.linspace(1.0/50.0/365.0, 1.0e-2, 995, dtype=np.double)
myt1 = pdm.getPDM_EquiBinCover(myp)
myt2 = pdm.getPDM_EquiBin(myp)

#myt3 = pdm.getPDM(myp)

#plt.plot(myp, myt1)
#plt.plot(myp, myt2)

#
# another light curve 
#
#con = np.loadtxt("con_all_combine.txt")
#pdm2 = cypdm.PyPDM(con[:, 0], con[:, 1], 10, 3)
#dm2.setScan(12, 4)

#con_com = np.loadtxt("con_all_combine.txt", usecols=(0, 1))
#pdm.setData(con_com[:, 0], con_com[:, 1])
#myt1 = pdm2.getPDM_EquiBinCover(myp)
#myt2 = pdm2.getPDM_EquiBin(myp)

#plt.plot(myp, myt1)
#plt.plot(myp, myt2)

#plt.show()

