from calendar import c
from re import U
import numpy as np
import pandas as pd
from regex import B
from scipy.optimize import linprog

def model(t,a,b,c,d):
    return a + b * t + c * (t**2) + d * (t**3)

def fi(x,i,y,t):
    res = model(t[i],*x) - y[i]
    return 0.5 * (res**2)

def mount_Idelta(f,ind,q,delta):
    I = []
    m = 0
    fq = f[q]
    for i in range(samples):
        if np.abs(fq - f[i]) <= delta:
            I.append(ind[i])
            m += 1
    
    return I, m

df = pd.read_csv("output/original.txt",header=None, sep=" ")

# Set parameters
n = 5
samples = df[0].size
epsilon = 1.e-3
delta = 0.001
theta = 0.5
q = 35

t       = np.zeros(samples)
y       = np.zeros(samples)
indices = np.zeros(samples,dtype=int)
xk      = np.array([-1.,-2.,1.,-1.])
xtrial  = np.zeros(n-1)
faux    = np.zeros(samples)
c       = np.zeros(n)
A       = np.zeros((samples,n))
b       = np.zeros(samples)
x_bounds= np.zeros((n,2))
grad    = np.zeros((samples,n-1))

for i in range(samples):
    t[i] = df[0].values[i]
    y[i] = df[1].values[i]

iter = 0

for i in range(samples):
    faux[i] = fi(xk,i,y,t)

indices = np.argsort(faux)
faux = np.sort(faux)

Idelta, m = mount_Idelta(faux,indices,q,delta)

while True:
    iter += 1

    for i in range(m):
        ti          = t[indices[i]]
        gaux        = model(ti,*xk) - y[Idelta[i]]
        grad[i,0]   = 1.0
        grad[i,1]   = ti
        grad[i,2]   = ti**2
        grad[i,3]   = ti**3
        grad[i,:]   = gaux * grad[i,:]


    # Objective function
    c[n-1] = 1.

    # Constraints matrix
    for i in range(m):
        A[i,:n-1]   = grad[i,:]
        A[i,n-1]    = -1.

    x_bounds[:n-1,0] = np.maximum(-10. - xk,-np.ones(n-1))
    x_bounds[:n-1,1] = np.minimum(10. - xk,np.ones(n-1))
    x_bounds[n-1,0] = -np.inf
    x_bounds[n-1,1] = np.inf

    res = linprog(c, A_ub=A[:m,:], b_ub=b[:m], bounds=x_bounds)

    print(res.fun)

    break