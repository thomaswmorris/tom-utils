


# Utilities for time-ordered data. Mostly for ACT. 

import numpy as np
import numpy.linalg as la
import scipy as sp
from scipy import special, stats, signal

msqrt  = lambda M: [np.matmul(u,np.diag(np.sqrt(s))) for u,s,vh in [la.svd(M)]][0]
matern = lambda r,r0,nu : 2**(1-nu)/sp.special.gamma(nu)*sp.special.kv(nu,r/r0+1e-10)*(r/r0+1e-10)**nu


class generator():

    def __init__(self, *params, method='outer', cfun=None):

        self.params = params
        self.method = method
        
        if not np.all([self.params[0].shape==p.shape for p in self.params]):
            raise(Exception('oops'))

        if not len(self.params[0].shape) == 1:
            self.param_shape = self.params[0].shape
            self.params = [p.ravel() for p in self.params]

        if cfun == 'matern': cfun = lambda r,r0,nu : 2**(1-nu)/sp.special.gamma(nu)*sp.special.kv(nu,r/r0+1e-10)*(r/r0+1e-10)**nu

        if method == 'outer':
            P = []; [P.extend(np.meshgrid(_p,_p)) for _p in self.params]
            self.cov = cfun(*P)
            self.gen = msqrt(self.cov)
            del P 



    def generate(self):

        if self.method == 'outer': return np.matmul(self.gen, np.random.standard_normal(self.params[0].shape))






time_ = np.linspace(0,60,600+1)
x_ = 2 * time_
y_ = 3 * time_

cfun = lambda xi,yi,xj,yj : matern(np.sqrt((xi-xj)**2+(yi-yj)**2),1,5/6)


G = generator(x_,y_,cfun=cfun)


print(G.generate().shape)
