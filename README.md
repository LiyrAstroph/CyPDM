# CyPDM

A fast package to apply the phase disperion minimization (PDM) algorithm, based on PyAstronomy module PyPDM (https://github.com/sczesla/PyAstronomy).

The package is written in C and wrapped with Cython, so that it can be called in both C and Python code.

The PDM alogrithm refers to the reference:
http://adsabs.harvard.edu/abs/1978ApJ...224..953S

 
Author:
Yan-Rong Li, liyanrong@mail.ihep.ac.cn

# Compiling

```bash
python setup.py build_ext -i
```
This will generate a library ``cypdm.so``, which can be directly imported in Python.

# Usage

```python
import numpy as np
import cypdm 

con = np.loadtxt("con_all.txt", usecols=(0, 1))
t0 = con[0, 0]
con[:, 0] -= t0

pdm = cypdm.PyPDM(con[:, 0], con[:, 1], 10, 3)
myp = 1.0/np.linspace(1.0/50.0/365.0, 1.0e-2, 995, dtype=np.double)
myt1 = pdm.getPDM_EquiBinCover(myp)
myt2 = pdm.getPDM_EquiBin(myp)
```

This loads data from file "con_all.txt" and calculates PDMs.

# Comparison with PyPDM

PyPDM is distributed in PyAstronomy package, which written in pure Python with severe additional overhead.

![Comparison between CyPDM and PyPDM](https://github.com/liyropt/MyGithubPic/blob/master/cypdm_cmp.jpg)

Running the test codes of CyPDM in the test directory of this package using my own desktop 
```bash
time python test_cypdm.py
```
gives the time statistics
```bash
real	0m0.250s
user	0m0.380s
sys	0m0.394s
```
As a comparison, running the test codes of PyPDM
```bash
time python test_pya.py
```
gives 
```bash
real	0m0.809s
user	0m1.042s
sys	0m0.839s
```
CyPDM speeds up by a factor of 3-4. (Of course, the time statistics depend on data size.)