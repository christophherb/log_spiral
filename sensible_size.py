from new_log_spir.LogSpir import LogSpir, Neutron
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("/home/cherb/LRZSync/Doktorarbeit/Vorlagen_Verschiedenes/stylelibs/christoph.mplstyle")

p=1
if p:
    log = LogSpir(0.5, 0.675, 1.55, 4, precision=1e-8)
    print(log.theta_end)
    fig, ax = plt.subplots(1, figsize=(15, 15))

    theta_end = log.theta_end
    theta_range = np.linspace(0, theta_end, 10001)
    z, x = log.return_cart_coords(theta_range)
    for slope in np.linspace(0, 0.030, 51):
        neutron = Neutron(0.7, slope*0.7, -1, -slope)
        path, neutron = log.propagate_neutron(neutron)
        z_path, x_path = zip(*[k[:2] for k in path])
        #ax.set_aspect('equal')
        
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
        
        #ax.plot([log.zstart, log.zend], [0, log.xend], label='approx', linestyle='-', marker=' ')#approximation of the non-rotated mirror, determines
        ax.plot(z_r, x_r, color='black', label='mirror', marker=' ', linestyle='-', linewidth=2)
        #the initial mirror
    ax.set_xlabel("$z$ (m)")
    ax.set_ylabel("$x$ (m)")
    ax.set_title(r"$r = z_s \exp(\cotan(\psi)\theta)$  mit $z_s = {}$ m, $z_e = {}$ m, $\psi={}^\circ$".format(log.zstart, log.zend, log.psi))
    fig.savefig('/home/cherb/Downloads/sensible_size_multimirror.pdf', bbox_inches='tight')
    plt.show()