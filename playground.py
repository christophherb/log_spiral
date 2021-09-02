from new_log_spir.LogSpir import LogSpir, Neutron
import numpy as np
import matplotlib.pyplot as plt

p=True
if p:
    log = LogSpir(1, 4, 5, 4, precision=1e-8)
    print(log.theta_end)
    neutron = Neutron(4, 1.5, -1, -1)
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
    ax.plot([neutron.z, neutron.z + log.zend], [neutron.x, neutron.x+neutron.vx/neutron.vz*log.zend], label='incoming')
    z_path, x_path = zip(*[k[:2] for k in path])
    ax.set_aspect('equal')
    ax.plot(z_path, x_path, ls='-', marker=' ')
    plt.show()

p=True
if p:
    log = LogSpir(1, 4, 5, 4, precision=1e-8)
    neutron = Neutron(0, 0.35, 1, -0.1)
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
    ax.plot([neutron.z, neutron.z + log.zend], [neutron.x, neutron.x+neutron.vx/neutron.vz*log.zend], label='incoming')
    z_path, x_path = zip(*[k[:2] for k in path])
    ax.set_aspect('equal')
    ax.plot(z_path, x_path, ls='-', marker=' ')
    plt.show()
