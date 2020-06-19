

on osx with git-installed  eigency
```
pip install git+https://github.com/ericmjonas/eigency.gti
```

```
CC=gcc-6 python setup.py install
```

works. 

With eigency-1.75 from pip, 

```
CC=gcc-6 python setup.py install
```

works. 

To get this to work with clang required adding
```
 '-stdlib=libc++'
```


### Trying to debug the eigency madness

pip install git+https://github.com/wouterboomsma/eigency.git

Build fails with:
/data/jonas/anaconda/lib/python2.7/site-packages/eigency/eigen_3.2.8/Eigen/Core:15:49: fatal error: src/Core/util/DisableStupidWarnings.h: No such file or directory

jonas@c65:/data/jonas/pychirpz$ ls /data/jonas/anaconda/lib/python2.7/site-packages/eigency
conversions_api.h  conversions.so  core.so        __init__.py
conversions.pxd    core.pxd        eigen_3.2.8    __init__.pyc
conversions.pyx    core.pyx        eigency_cpp.h
jonas@c65:/data/jonas/pychirpz$ ls /data/jonas/anaconda/lib/python2.7/site-packages/eigency/eigen_3.2.8/
Eigen
jonas@c65:/data/jonas/pychirpz$ ls /data/jonas/anaconda/lib/python2.7/site-packages/eigency/eigen_3.2.8/Eigen/
Array           Geometry                PardisoSupport   SparseQR
Cholesky        Householder             PaStiXSupport    SPQRSupport
CholmodSupport  IterativeLinearSolvers  QR               StdDeque
Core            Jacobi                  QtAlignedMalloc  StdList
Dense           LeastSquares            Sparse           StdVector
Eigen           LU                      SparseCholesky   SuperLUSupport
Eigen2Support   MetisSupport            SparseCore       SVD
Eigenvalues     OrderingMethods         SparseLU         UmfPackSupport


now we uninstall

pip uninstall eigency


now try installing from the locally cloned dir with python setup.py install
git clone
python setup.py install


Then our package succeeds (and eignecy is installed via an egg)

What if we install from the locally-cloned eigency git repo but via 
pip install eigency


Same error! 






## Optimization
Start: 6.2ms/eval
switch to block operations: 5.3ms/eval 
