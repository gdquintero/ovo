import numpy as np
import matplotlib.pyplot as plt

def model(t,a,b,c,d):
    return a * np.exp(b * np.sin(c * t)) + d

t_data = np.linspace(0, 1.5 * np.pi,50)

x = np.empty(4)
x = (1.,0.5,2.,-1.)
noise_factor = 0.04
outliers = np.empty(3)
ind = []

outliers[0] = model(np.pi / 4.,*x) + 0.5
outliers[1] = model(3. * np.pi / 4.,*x) - 0.5
outliers[2] = model(5. * np.pi / 4.,*x) + 0.5

y = model(t_data,*x)
rng = np.random.default_rng()
y_noise = noise_factor * rng.normal(size=t_data.size)
y_data = y + y_noise
j = 0

for i in range(1,6,2):
    ind.append(np.argmin(abs(t_data - i * np.pi/4.)))
    y_data[ind[j]] = outliers[j]
    j += 1

with open("output/outlier_indices_wave.txt",'w') as f:
    for i in range(3):
        f.write("%i\n" % ind[i])

with open("output/wave.txt",'w') as f:
    for i in range(t_data.size):
        f.write("%f %f\n" % (t_data[i],y_data[i]))

plt.plot(t_data,y)
plt.plot(t_data,y_data,"o")
plt.show()