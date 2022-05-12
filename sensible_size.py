from new_log_spir.LogSpir import LogSpir, Neutron
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("/home/cherb/LRZSync/Doktorarbeit/Vorlagen_Verschiedenes/stylelibs/christoph.mplstyle")
from random import normalvariate, uniform
p=1
test_log = LogSpir(0.15+0.025, 0.3-0.025, 1.45, 6, precision=1e-8, phi_rot=0)

if p:
    log = LogSpir(0.15+0.025, 0.3-0.025, 1.45, 6, precision=1e-8, phi_rot=test_log.theta_end*0.35)
    log.branches = 3
    while True:
        angle = (log.branches-1)*log.phi_rot+log.theta_end
        height = np.sin(angle)*log.zend

        print(angle*180/np.pi, height)
        print(height/log.zend < 0.023/(0.675-0.06))
        if height/log.zend < 0.023/(0.675-0.06):
            log.branches += 3
        else:
            break
    print(log.branches, height/log.zend, 0.02/(0.675-0.06))
    log.double_sided = 1
    log.m = 6.2
    print(log.theta_end)
    fig, ax = plt.subplots(1, figsize=(7, 7))
    v0 = 980
    theta_end = log.theta_end
    theta_range = np.linspace(0, theta_end, 501)
    z, x = log.return_cart_coords(theta_range)
    hist = []
    for slope in np.linspace(-0.0324, 0.0324, 101):
        z0 = 0.675-0.06
        vz = -1
        vx = -slope
        x0 = slope*z0#+(uniform(-0.5, 0.5))*0.002
        vz, vx = v0*vz/(vx**2+vz**2)**0.5, v0*vx/(vx**2+vz**2)**0.5
        neutron = Neutron(z0, x0, vz, vx)
        path, neutron = log.propagate_neutron(neutron)
        z_path, x_path = zip(*[k[:2] for k in path])
        #ax.set_aspect('equal')

        if neutron.vz > 0:
            z_path=list(z_path)+[2*log.zend]
            x_path=list(x_path)+[neutron.x+neutron.vx/neutron.vz*(2*log.zend-neutron.z)]
            ax.plot([neutron.z, 2*log.zend], [neutron.x, neutron.x+neutron.vx/neutron.vz*\
                (2*log.zend-neutron.z)], marker=' ')
        else:
            z_path=list(z_path)+[0]
            x_path=list(x_path)+[neutron.x+neutron.vx/neutron.vz*(0-neutron.z)]
            ax.plot([neutron.z, 0], [neutron.x, neutron.x+neutron.vx/neutron.vz*\
                (0-neutron.z)], marker=' ')
        hist += [x_path[-1]]
        ax.plot(z_path, x_path, ls='-', marker=' ')

    for branch in range(log.branches):#plotting all spirals of the log spiral for testing
        cos = np.cos(branch*log.phi_rot)
        sin = np.sin(branch*log.phi_rot)
        z_r, x_r = cos*z - sin*x, sin*z + cos*x

        #ax.plot([log.zstart, log.zend], [0, log.xend], label='approx', linestyle='-', marker=' ')#approximation of the non-rotated mirror, determines
        ax.plot(z_r, x_r, color='black', label='mirror', marker=' ', linestyle='-', linewidth=2)
        if log.double_sided == 1:
            ax.plot(z_r, -x_r, color='black', label='mirror', marker=' ', linestyle='-', linewidth=2)
        #the initial mirror
    ax.set_xlabel("$z$ (m)")
    ax.set_ylabel("$x$ (m)")
    ax.set_title(r"$r = z_s \exp(\cotan(\psi)\theta)$  mit $z_s = {:.3}$ m, $z_e = {:.3}$ m, $\psi={}^\circ$".format(log.zstart, log.zend, log.psi))
    ax.legend()
    fig.savefig('/home/cherb/Downloads/sensible_size_divergence.pdf', bbox_inches='tight')
    plt.show()

    # print the histogram at the
    fig, ax = plt.subplots(1)
    hist = np.histogram(hist, bins = 1001, range=(-0.013, 0.013))
    ax.plot(np.array(hist[1][:-1]+0.5*(hist[1][1]-hist[1][0]))*100, hist[0], linestyle='-', marker=' ')
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("intensity arb. u.")
    plt.show()
p=0
if p:
    log = LogSpir(0.175, 0.275, 1.55, 6, precision=1e-8)
    log.double_sided = 0
    print(log.theta_end)
    fig, ax = plt.subplots(1, figsize=(7, 7))

    theta_end = log.theta_end
    theta_range = np.linspace(0, theta_end, 1001)
    z, x = log.return_cart_coords(theta_range)
    for slope in np.linspace(-0.03, 0.030, 501):
        neutron = Neutron(0, 0, 1, slope)
        path, neutron = log.propagate_neutron(neutron)
        z_path, x_path = zip(*[k[:2] for k in path])
        ax.set_aspect('equal')

        if neutron.vz > 0:
            z_path=list(z_path)+[2*log.zend]
            x_path=list(x_path)+[neutron.x+neutron.vx/neutron.vz*(2*log.zend-neutron.z)]
            ax.plot([neutron.z, 2*log.zend], [neutron.x, neutron.x+neutron.vx/neutron.vz*(2*log.zend-neutron.z)])
        else:
            z_path=list(z_path)+[0]
            x_path=list(x_path)+[neutron.x+neutron.vx/neutron.vz*(0-neutron.z)]
            ax.plot([neutron.z, 0], [neutron.x, neutron.x+neutron.vx/neutron.vz*(0-neutron.z)])

        ax.plot(z_path, x_path, ls='-', marker=' ')
    for branch in range(log.branches):#plotting all spirals of the log spiral for testing
        cos = np.cos(branch*theta_end)
        sin = np.sin(branch*theta_end)
        z_r, x_r = cos*z - sin*x, sin*z + cos*x

        #ax.plot([log.zstart, log.zend], [0, log.xend], label='approx', linestyle='-', marker=' ')#approximation of the non-rotated mirror, determines
        ax.plot(z_r, x_r, color='black', label='mirror', marker=' ', linestyle='-', linewidth=2)
        #the initial mirror
    ax.set_xlabel("$z$ (m)")
    ax.set_ylabel("$x$ (m)")
    ax.set_title(r"$r = z_s \exp(\cotan(\psi)\theta)$  mit $z_s = {}$ m, $z_e = {}$ m, $\psi={}^\circ$".format(log.zstart, log.zend, log.psi))
    #fig.savefig('/home/cherb/Downloads/sensible_size.pdf', bbox_inches='tight')
    plt.show()