import matplotlib.pyplot as plt
import numpy as np

data = "Rg&Energy.dat"
step = []
rg = []
total = []

with open(data,"r") as file:
    for line in file:
        values = line.split()
        step.append(float(values[0]))
        rg.append(float(values[1]))
        total.append(float(values[2]))
        
var_rg = np.var(rg)
var_total = np.var(total)

plt.plot(step,rg,marker='o')
plt.xlabel('Replicas')
plt.ylabel('Radius of Gyration')
plt.title('Radius of Gyration of EWS/FLI1\nVariance: {:.2f}'.format(var_rg))
x_ticks_int = range(int(min(step)), int(max(step))+1)
plt.xticks(x_ticks_int)
x_ticks = np.arange(min(step)-1,max(step)+1,2)
plt.xticks(x_ticks)
# plt.grid(True)
# plt.show()
plt.savefig('2nd_Analysis//Rg_replicas.png',dpi=300)
plt.close()
plt.plot(step,total,marker='o')
plt.xlabel('Replicas')
plt.ylabel('Potential Energy')
plt.title('Potential Energy of EWS/FLI1\nVariance: {:.2f}'.format(var_total))
x_ticks = np.arange(min(step)-1,max(step)+1,2)
plt.xticks(x_ticks)
# plt.show()
plt.savefig('2nd_Analysis//TotalEnergy_replicas.png',dpi=300)