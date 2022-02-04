from new_log_spir.LogSpir import LogSpir, Neutron
import numpy as np
import matplotlib.pyplot as plt

p=False
if p:
    log = LogSpir(1, 3, 5, 4, precision=1e-8)
    print(log.theta_end)
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
    for slope in np.linspace(0, 0.3, 21):
        neutron = Neutron(0, 0, 1, slope)
        path, neutron = log.propagate_neutron(neutron)
        z_path, x_path = zip(*[k[:2] for k in path])
        ax.set_aspect('equal')
        ax.plot(z_path, x_path, ls='-', marker=' ')
        ax.plot([neutron.z, 2*log.zend], [neutron.x, neutron.x+neutron.vx/neutron.vz*(2*log.zend-neutron.z)])
    fig.savefig('/home/cherb/Downloads/plsdelet.pdf', bbox_inches='tight')
    plt.show()

p=True
if p:
    log = LogSpir(1, 3, 7, 2, precision=1e-8)
    neutron = Neutron(0, 0.299755, 761.892137, -225.713987     )
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
    print('this', neutron.vz, neutron.vx)
    ax.plot([neutron.z, neutron.z + log.zend], [neutron.x, neutron.x+neutron.vx/neutron.vz*log.zend], label='incoming')
    z_path, x_path = zip(*[k[:2] for k in path])
    ax.set_aspect('equal')
    ax.plot(z_path, x_path, ls='-', marker=' ')
    plt.show()
