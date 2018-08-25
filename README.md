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

nd = 1000
P = 10.0
jd = np.linspace(0.0, 100.0, nd, dtype=np.double)
fs = np.sin(jd/P * 2.0*np.pi) + np.random.randn(nd)*0.1

# create a CyPDM instance with 10 bins and 3 covers
pdm = cypdm.CyPDM(jd, fs, 10, 3)
# setup period grid
myp = 1.0/np.linspace(1.0/50.0, 1.0e-2, 10000, dtype=np.double)
# get PDM using equidistant bins and covers
myt1 = pdm.getPDM_EquiBinCover(myp)
# get PDM using equidistant bins (no covers)
myt2 = pdm.getPDM_EquiBin(myp)
# get PDM according to covers, if covers==0, using getPDM_EquiBin(), and if covers!=0, using getPDM_EquiBinCover()
myt3 = pdm.getPDM(myp)


jd_new = np.linspace(0.0, 100.0, nd, dtype=np.double)
fs_new = np.sin(jd/P * 2.0*np.pi) + np.random.randn(nd)*0.1
# set new data 
pdm.setData(jd_new, fs_new)

# set new bins 10 and covers 5
pdm.setScan(12, 5)
```

# Comparison with PyPDM

PyPDM is distributed in PyAstronomy package, which written in pure Python with severe additional overhead.

Running the test code cmp_with_pya.py in the test directory gives the figure shown below.

![Comparison between CyPDM and PyPDM](https://github.com/liyropt/MyGithubPic/blob/master/cypdm_cmp.jpg)

Running the test codes of CyPDM in the test directory of this package using my own desktop 
```bash
time python test_cypdm.py
```
gives the time statistics
```bash
real    0m1.714s
user    0m1.897s
sys 0m0.498s
```
As a comparison, running the test codes of PyPDM
```bash
time python test_pya.py
```
gives 
```bash
real    0m8.394s
user    0m8.634s
sys 0m0.909s
```
CyPDM speeds up by a factor of 3-4. (Of course, the time statistics depend on data size.)