import matplotlib
from new_log_spir.LogSpir import LogSpir, Neutron
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("/home/cherb/LRZSync/Doktorarbeit/Vorlagen_Verschiedenes/stylelibs/christoph.mplstyle")
#from point radially diverging
p=1
if p:
    log = LogSpir(1, 3, 5, 5, phi_rot=0.0965243188196*180/3.1415, precision=1e-6)
    log.double_sided = 1
    print(log.theta_end)
    fig, ax = plt.subplots(1, figsize=(6, 4))

    theta_end = log.theta_end
    theta_range = np.linspace(0, theta_end, 10001)
    z, x = log.return_cart_coords(theta_range)

    for slope in np.linspace(-0.5, 0.5, 51):
        neutron = Neutron(6, -slope*6, -1, slope)
        #print('incoming {} {}'.format(neutron.vz, neutron.vx))
        path, neutron = log.propagate_neutron(neutron)
        z_path, x_path = zip(*[k[:2] for k in path])
        ax.set_aspect('equal')

        if neutron.vz > 0:
            z_path=list(z_path)+[2*log.zend]
            x_path=list(x_path)+[neutron.x+neutron.vx/neutron.vz*(2*log.zend-neutron.z)]
            #ax.plot([neutron.z, 2*log.zend], [neutron.x, neutron.x+neutron.vx/neutron.vz*(2*log.zend-neutron.z)])
        else:
            z_path=list(z_path)+[0]
            x_path=list(x_path)+[neutron.x+neutron.vx/neutron.vz*(0-neutron.z)]
            #ax.plot([neutron.z, 0], [neutron.x, neutron.x+neutron.vx/neutron.vz*(0-neutron.z)])

        ax.plot(z_path, x_path, ls='-', marker=' ')

    for branch in range(log.branches):#plotting all spirals of the log spiral for testing
        cos = np.cos(branch*theta_end)
        sin = np.sin(branch*theta_end)
        z_r, x_r = cos*z - sin*x, sin*z + cos*x

        ax.plot([log.zstart, log.zend], [0, log.xend], label='approx', linestyle='-', marker=' ')#approximation of the non-rotated mirror, determines
        ax.plot(z_r, x_r, color='black', label='mirror', marker=' ', linestyle='-', linewidth=2)
        if log.double_sided:
            cos = np.cos(branch*theta_end)
            sin = np.sin(branch*theta_end)
            z_r, x_r = cos*z - sin*x, sin*z + cos*x
            ax.plot(z_r, -x_r, color='black', label='mirror', marker=' ', linestyle='-', linewidth=2)    
        #the initial mirror
    ax.set_xlabel("$z$ (m)")
    ax.set_ylabel("$x$ (m)")
    ax.set_title(r"$r = z_s \exp(\cotan(\psi)\theta)$  mit $z_s = 1$ m, $\psi=5^\circ$")
    fig.savefig('/home/cherb/Downloads/log_strahlengänge.pdf', bbox_inches='tight')
    plt.show()

#a cover of neutrons coming from the outside
p=0
if p:
    log = LogSpir(1, 3, 5, 5, precision=1e-8)
    print(log.theta_end)
    fig, ax = plt.subplots(1, figsize=(15, 15))

    theta_end = log.theta_end
    theta_range = np.linspace(0, theta_end, 10001)
    z, x = log.return_cart_coords(theta_range)
    for slope in np.linspace(0, 0.5, 501):
        neutron = Neutron(4, slope*4, -1, -slope)
        path, neutron = log.propagate_neutron(neutron)
        z_path, x_path = zip(*[k[:2] for k in path])
        ax.set_aspect('equal')

        if neutron.vz > 0:
            z_path=list(z_path)+[2*log.zend]
            x_path=list(x_path)+[neutron.x+neutron.vx/neutron.vz*(2*log.zend-neutron.z)]
            #ax.plot([neutron.z, 2*log.zend], [neutron.x, neutron.x+neutron.vx/neutron.vz*(2*log.zend-neutron.z)])
        else:
            z_path=list(z_path)+[0]
            x_path=list(x_path)+[neutron.x+neutron.vx/neutron.vz*(0-neutron.z)]
            #ax.plot([neutron.z, 0], [neutron.x, neutron.x+neutron.vx/neutron.vz*(0-neutron.z)])

        ax.plot(z_path, x_path, ls='-', marker=' ')
    for branch in range(log.branches):#plotting all spirals of the log spiral for testing
        cos = np.cos(branch*theta_end)
        sin = np.sin(branch*theta_end)
        z_r, x_r = cos*z - sin*x, sin*z + cos*x

        ax.plot([log.zstart, log.zend], [0, log.xend], label='approx', linestyle='-', marker=' ')#approximation of the non-rotated mirror, determines
        ax.plot(z_r, x_r, color='black', label='mirror', marker=' ', linestyle='-', linewidth=2)
        #the initial mirror
    ax.set_xlabel("$z$ (m)")
    ax.set_ylabel("$x$ (m)")
    ax.set_title(r"$r = z_s \exp(\cotan(\psi)\theta)$  mit $z_s = 1$ m, $\psi=2^\circ$")
    fig.savefig('/home/cherb/Downloads/log_strahlengänge_rückseite_2g.pdf', bbox_inches='tight')
    plt.show()
#a single ray for comparison
p=0
if p:
    log = LogSpir(1, 3, 7, 5, precision=1e-8)
    neutron = Neutron(4, 0 , -1143.311567, 290.297801  )
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
