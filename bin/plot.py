import sys
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt(sys.argv[1], delimiter=',')

step = data[:,0]
potential_energy = data[:,1]
temperature = data[:,2]
volume = data[:,3]
clipped_step = (data[:,0] > sys.argv[2])
clipped_data = data[clipped_step]

clipped_pe = np.array(clipped_data[0:,1])
pe_mean = np.mean(clipped_pe, axis=0)
pe_stdev = np.std(clipped_pe, axis=0)


plt.plot(step, potential_energy)
plt.xlabel("Step")
plt.ylabel("Potential energy (kJ/mol)")
plt.savefig('PE_plot.png')
plt.show()
plt.close()
plt.plot(step, temperature)
plt.xlabel("Step")
plt.ylabel("Temperature (K)")
plt.savefig('Temp_plot.png')
plt.show()
plt.close()
plt.plot(step, volume)
plt.xlabel("Step")
plt.ylabel("Volume (nm^3)")
plt.savefig('Vol_plot.png')
plt.show()
plt.close()
