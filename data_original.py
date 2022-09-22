import numpy as np
import random
import matplotlib.pyplot as plt

def model(t,a,b,c,d):
    return a + b * t + c * (t**2) + d * (t**3) 

t_data = np.linspace(-1,3.5,46)

x = np.array((0.,2.,-3.,1.))

y = model(t_data,*x)
y_noise = np.ones(t_data.size)

for i in range(t_data.size):
    y_noise[i] = y_noise[i] * random.uniform(-0.01,0.01)

y_data = y + y_noise

for i in range(6,16):
    y_data[i] = 10.

with open("output/original.txt",'w') as f:
    for i in range(t_data.size):
        f.write("%f %f\n" % (t_data[i],y_data[i]))

plt.plot(t_data,y)
plt.plot(t_data,y_data,"o")
plt.show()