import numpy as np
import matplotlib.pyplot as plt
from new_log_spir.LogSpir import LogSpir, Neutron

log = LogSpir(1, 3, 13, 2, precision=1e-8)
neutron = Neutron(0, 0.2, 1, 0.05)
fig, ax = plt.subplots(1)
theta_end = log.theta_end
theta_range = np.linspace(0, theta_end, 10001)
z, x = log.return_cart_coords(theta_range)
for branch in range(log.branches):#plotting all spirals of the log spiral for testing
    cos = np.cos(branch*theta_end)
    sin = np.sin(branch*theta_end)
    z_r, x_r = cos*z - sin*x, sin*z + cos*x
    ax.plot(z_r, x_r, color='black', label='mirror')
    ax.plot([log.zstart, log.zend], [0, log.xend], label='approx')#approximation of the non-rotated mirror, determines
    #the initial mirror
path, neutron = log.propagate_neutron(neutron)
ax.plot([neutron.z, 6], [neutron.x, neutron.x+neutron.vx/neutron.vz*(6-neutron.z)], label='incoming')
z_path, x_path = zip(*[k[:2] for k in path])
ax.set_aspect('equal')
ax.plot(z_path, x_path, ls='-', marker=' ')
plt.show()