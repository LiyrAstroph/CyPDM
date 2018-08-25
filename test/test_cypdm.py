import sys
sys.path.append("../")
import cypdm 
import numpy as np 
#import matplotlib.pyplot as plt 

con = np.loadtxt("con_all.txt", usecols=(0, 1))
t0 = con[0, 0]
con[:, 0] -= t0

pdm = cypdm.PyPDM(con[:, 0], con[:, 1], 10, 3)
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

