import numpy as np
import matplotlib.pyplot as plt

#use theta and phi to get random distribution of sources across the sky
#theta ranges 0 to pi
#phi ranges 0 to 2pi

phi_array = np.random.uniform (0, 2*np.pi, 1000) 
#print(phi_array)

#print(np.mean(phi_array))


costheta_array = np.random.uniform (-1, 1, 1000)
#print(costheta_array)
#print(np.mean(costheta_array))

theta_array = np.arccos(costheta_array) 
#print(theta_array)
#print(np.mean(theta_array))

#create figure
#myfig,ax = plt.subplots(projection='3d')
myfig = plt.figure()
ax=myfig.add_subplot(projection='3d')


ax.scatter(np.sin(theta_array)*np.cos(phi_array), np.sin(theta_array)*np.sin(phi_array), costheta_array)
plt.show()
