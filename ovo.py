import numpy as np
import pandas as pd
from scipy.optimize import linprog

def model(ti,x):
    return x[0] + (x[1] * ti) + (x[2] * (ti**2)) + (x[3] * (ti**3))

def fi(x,i):
    res = model(t[i],x) - y[i]
    return  0.5 * res**2

def mount_Idelta(f,ind,q,delta):
    I = []
    m = 0
    fq = f[q]

    for i in range(samples):
        if abs(fq - f[i]) <= delta:
            I.append(ind[i])
            m += 1
    
    return I, m

df = pd.read_csv("output/original.txt",header=None, sep=" ")

# Set parameters
n = 5
q = 35
samples = df[0].size
epsilon = 1.e-5
delta = 0.001
theta = 0.5
sigmin = 0.1
sigmax = 0.9
max_iter = 1000
max_int_iter = 100

t       = np.zeros(samples)
y       = np.zeros(samples)
indices = np.zeros(samples,dtype=int)
xtrial  = np.zeros(n-1)
faux    = np.zeros(samples)
c       = np.zeros(n)
A       = np.zeros((samples,n))
b       = np.zeros(samples)
grad    = np.zeros((samples,n-1))
x_bounds= [(None,None)] * n

# Initial solution
xk = np.array([-1.,-2.,1.,-1.])

for i in range(samples):
    t[i] = df[0].values[i]
    y[i] = df[1].values[i]

iter = 0

for i in range(samples):
    faux[i] = fi(xk,i)

indices = np.argsort(faux)
faux = np.sort(faux)
fxk = faux[q]

Idelta, m = mount_Idelta(faux,indices,q,delta)

while True:
    iter += 1

    for i in range(m):
        ti          = t[Idelta[i]]
        gaux        = model(ti,xk) - y[Idelta[i]]
        grad[i,:]   = gaux * np.array([1., ti, ti**2, ti**3])

    # Objective function
    c[-1] = 1.

    # Constraints matrix
    for i in range(m):
        A[i,:n-1]  = grad[i,:]
        A[i,-1]    = -1.

    # Box constraints
    for i in range(n-1):
        x_bounds[i] = (max(-10. - xk[i],-1.),min(10. - xk[i],1.))

    x_bounds[-1] = (None,0.)

    # Solve with linprog
    res = linprog(c, A_ub=A[:m,:], b_ub=b[:m], bounds=x_bounds)

    alpha = 1.
    int_iter = 1

    # Backtracking
    while True:
        xtrial = xk + alpha * res.x[:n-1]

        for i in range(samples):
            faux[i] = fi(xtrial,i)

        indices = np.argsort(faux)
        faux    = np.sort(faux)
        fxtrial = faux[q]
        
        if fxtrial <= fxk + theta * alpha * res.fun: break
        if int_iter >= max_int_iter: break

        alpha = 0.5 * (sigmax - sigmin) * alpha
        int_iter += 1


    print(iter,int_iter,fxk,abs(res.fun))

    if abs(res.fun) <= epsilon: break
    if iter >= max_iter: break

    xk = xtrial
    fxk = fxtrial

    Idelta, m = mount_Idelta(faux,indices,q,delta)

print(xk)