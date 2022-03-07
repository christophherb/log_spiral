import numpy as np
import matplotlib.pyplot as plt
import random as rn
rad2deg = 180/np.pi
deg2rad = 1/rad2deg
V2Q_conic = 1.58825361e-3

def calc_supermirror_reflectivity(m, q, alpha=2.5, W=0.004):
    Q_c = 0.0218
    arg = 0
    beta = 0
    R_0 = 0.995
    #alpha = 2.5
    #W = 0.004
    weight = 1
    if m >= 10:
        return weight
    q = abs(q)
    if (W==0 and alpha==0):
        m=m*0.9853+0.1978
        W=-0.0002*m+0.0022
        alpha=0.2304*m+5.0944
        beta=-7.6251*m+68.1137
        if (m<=3):
	        alpha=m
	        beta=0
    if W > 0:
        arg = (q - m*Q_c)/W
    else:
        arg = 11

    if (arg > 10 or m <= 0 or Q_c <=0 or R_0 <= 0):
        weight = 0
        return weight


    if (m < 1):
        Q_c *= m
        m=1

    if(q <= Q_c):
        weight = R_0
        return weight



    weight = R_0*0.5*(1 - np.tanh(arg))*(1 - alpha*(q - Q_c) + beta*(q - Q_c)*(q - Q_c))

    return weight

class Neutron:
    def __init__(self, z, x, vz, vx) -> None:
        """initializes a neutron with position and velocity

        Args:
            z (float): z-coordinate of the neutrons origin
            x (float): x-coordinate of the neutrons origin
            vz (float): velocity of the neutron in z-direction
            vx (float): velocity of the neutron in x-direction
        """
        self.z = z
        self.x = x
        self.vz = vz
        self.vx = vx

    def return_coords(self):
        """returns the neutronparameters to be used by other functions

        Returns:
            (tuple): tuple of neutron coordinates (z, x, vz, vx)
        """
        return self.z, self.x, self.vz, self.vx

    def prop_neutron(self, dt):
        """returns the coordinates of the neutron propagated by dt

        Args:
            dt (float): time by which the neutron flies

        Returns:
            (tuple): (z, x, vz, vx) neutron propagated by dt
        """
        return self.z + self.vz*dt, self.x+self.vx*dt, self.vz, self.vx

class LogSpir:
    def __init__(self, zstart: float, zend: float, psi: float, branches: int, precision: float = 1e-7, phi_rot: float = 0) -> None:
        """initializes the logarithmic mirror r(theta) = zstart * exp(k*theta), with k = cotan(psi)
           and r = sqrt(z*z + x*x)

        Args:
            zstart (float): where the logspiral cuts the zaxis.
            zend (float): at which z-value the log spiral ends, important for the linear approx.
            psi (float): the angle in degrees (deg) at which the spiral cuts the zaxis, equals the angle under
            which neutrons originating from the origin hit the mirror.
            branches (int): number of branches from the log spiral, how many mirrors.
            precision (float): precision to which the approximation methods are carried out. Defaults to 1e-7.
        """
        self.zstart = zstart
        self.zend = zend
        self.theta_start = 0#mirror starts at the z-axis or theta = 0
        self.psi = psi
        self.psi_rad = psi*deg2rad
        self.k = 1/np.tan(self.psi_rad) #helper variable showing in the formula
        self.precision=precision
        self.m = 6.2
        #function and derivative are used to find the angle theta at which the z_value of the spiral equals zend
        self.function = lambda theta: np.cos(theta)*self.zstart*np.exp(self.k*theta)-self.zend
        self.derivative = lambda theta: self.zstart*np.exp(self.k*theta)*(np.cos(theta)*self.k-np.sin(theta))
        self.theta_end = self.return_precise_theta_end(self.zend)
        self.phi_rot = self.theta_end if phi_rot == 0 else phi_rot*deg2rad
        self.xend = self.zend*np.tan(self.theta_end)
        self.double_sided = 1
        self.branches = branches

    def return_r(self, theta):
        """returns the distance from the origin to the spiral at a given angle theta according to r(theta) = zstart*exp(k*theta_rad)

        Args:
            theta (float): angle (deg)

        Returns:
            float: distance to the origin (m)
        """
        return self.zstart*np.exp(self.k*theta)

    def return_intersect_lines(self, x0, y0, vx0, vy0, x1, y1, vx1, vy1):
        """returns the intersection of 2, 2D lines, (x0, y0) + t*(vx0, vy0) == (x1, y1) + lam*(vx1, vy1)

        Args:
            x0 (float): [description]
            y0 (float): [description]
            vx0 (float): [description]
            vy0 (float): [description]
            x1 (float): [description]
            y1 (float): [description]
            vx1 (float): [description]
            vx2 (float): [description]
        """
        #find solutions to the equation line0 + linedir*t == log0 + logdir*lam and finally check if the x value makes sense
        try:
            t = (vy1*(x0-x1) + vx1*(y1-y0))/\
                (vx1*vy0 - vx0*vy1)
            x_intersection = x0 + vx0*t
            y_intersection = y0 + vy0*t
        except ZeroDivisionError: #special cases have to be solved in lambda not in t
            try:
                lam = (vx0*(y0-y1) + vy0*(x1-x0))/\
                    (vy1*vx0 - vx1*vy0)
                x_intersection = x1 + vx1*lam
                y_intersection = y1 + vy1*lam
            except ZeroDivisionError:#
                x_intersection = None
                return None
        return x_intersection, y_intersection
        #if x_intersection and (self.logspir.xstart < x_intersection < self.logspir.xend):
        #    return x_intersection
        #return None
    def return_cart_coords(self, theta):
        """returns cartesian coordinates z, x at a given theta value

        Args:
            theta (float): angle (deg)

        Returns:
            tuple: z, x (m, m)
        """
        r = self.return_r(theta)#only theta is needed as the logspiral is defined
        z = np.cos(theta)*r
        x = np.sin(theta)*r
        return z, x

    def return_approx_theta_end(self, zend=None):#those are only important for determining whether a neutron hits the mirror, can also compare to zend?
        """returns a close approx to the theta value at a specific zend
        Args:
            zend (float, optional): z coordinate at theta is determined. Defaults to None, self.zend.
        Returns:
            float: theta_end approximation
        """
        if zend == None:
            zend = self.zend
        prelim_theta = np.log(zend/self.zstart)*1/self.k# r in log spiral is approximated by zend (good enough in most cases)
        return prelim_theta

    def return_precise_theta_end(self, zend, max_iterations=10, precision=None):
        """uses newton raphson to calculate a precise value under which angle (theta) the end (at z=zend) of the log spiral appears
        Args:
            zend (float, optional): z coordinate at which the mirror ends.
            max_iterations (int, optional): maximum number of iterations, if no convergence is reached, return None. Defaults to 10.
            precision (float, optional): Precision of thetha end. Defaults to self.precision.

        Returns:
            float: theta_end such that the z coord at this angle equals zend
        """
        if precision == None:
            precision = self.precision
        theta_0 = self.return_approx_theta_end(zend)
        for _ in range(max_iterations):
            theta_n = theta_0 - self.function(theta_0)/self.derivative(theta_0)#function from above
            if abs(theta_0-theta_n) < precision:
                return theta_n
            theta_0 = theta_n
        return None

    def return_approx_theta_int(self, neutron: Neutron):
        """returns an approximate value for theta, approximating the mirror as a line from beginning of the mirror
        to its end

        Args:
            neutron (Neutron): incoming neutron to be refleted, could be exchanged for a tuple containing the origin and velocity
            of the neutrons

        Returns:
            theta (float): theta value of the intersection
        """
        #1st calc prelim intersection, then tan takes care of the rest
        z0, x0, vz0, vx0 = neutron.return_coords() #would it be better to feed the
        z1, x1, vz1, vx1 = self.zstart, 0, self.zend-self.zstart, self.xend
        zint, xint = self.return_intersect_lines(z0, x0, vz0, vx0, z1, x1, vz1, vx1)
        if zint < self.zstart or zint > self.zend:
            return None
        approx_theta = np.arctan(xint/zint)
        return approx_theta

    def return_precise_theata_int(self, z, x, zdir, xdir, max_iterations=20, precision=None):
        """returns the precise angle under which the neutron hits the logspiral

        Args:
            z (float): z-coordinate of the origin of the neutron
            x (float): x-coordinate of the origin of the neutron
            zdir (float): z-component of the velocity of the neutron
            xdir (float): x-component of the velocity of the neutron
            max_iterations (int, optional): Max iterations of Newton Raphson after which None is returned. Defaults to 10.
            precision (float, optional): Precision after which the function terminates successfully. Defaults to self.precision.

        Returns:
            float: theta angle of intersection
        """
        if precision == None:
            precision = self.precision
        k = self.k
        m = xdir/zdir #neutron.vx/neutron.vz, slope of the neutron path for calculation
        x0 = x+(0-z)/zdir*xdir #x-coordinate of the intersection of the neutron with the z axis

        def function(theta):#line intersects spiral where function(theta) = 0
            return self.zstart*np.exp(k*theta)*(np.sin(theta) - m*np.cos(theta))-x0
        def derivative(theta):#theta derivative of the above function for newton iterations
            return self.zstart*np.exp(k*theta)*(np.cos(theta)*(1-k*m)+np.sin(theta)*(k+m))

        theta0 = self.return_approx_theta_int(Neutron(z, x, zdir, xdir))
        if theta0 is None:
            return None
        for ind in range(max_iterations):
            thetan = theta0 - function(theta0)/derivative(theta0)
            if ind > 3 and (thetan < 0 or thetan > self.theta_end):
                print('inbetween?')
                return None
            if abs(thetan-theta0) < precision:
                if 0 <= thetan < self.theta_end:
                    return thetan
                return None
            theta0 = thetan
        print('non applicable')
        return None

    def return_inters_coords(self, neutron: Neutron):
        """returns the carthesian coordinates of the intersection

        Args:
            tuple (z, x): tuple of z and x coordinates of neutron intersection
        """
        z, x, vz, vx = neutron.return_coords()
        theta_int = self.return_precise_theata_int(z, x, vz, vx)
        if theta_int == None:
            return None
        z_int, x_int = self.return_cart_coords(theta_int)
        return z_int, x_int

    def return_normal_vec(self, theta):
        """returns the normal vector of the spiral at given theta,
        relevant for the calculation of the direction of the reflected vector

        Args:
            theta (float): angle of intersection (deg)

        Returns:
            tuple: nx,ny x and y component of the normal vector
        """
        k = self.k
        norm_prefac = 1/(1 + k**2)**0.5
        cos = np.cos(theta)
        sin = np.sin(theta)
        nz = (cos+k*sin)*norm_prefac
        nx = (sin-k*cos)*norm_prefac
        return nz, nx

    def return_reflected_dir(self, nz, nx, zdir, xdir):
        """returns the direction of the reflected neutron when given the norm-vector of the surface

        Args:
            nz (float): z-component of the normal-vector
            nx (float): x-component of the normal-vector
            zdir (float): z-direction of the incoming vector
            xdir (float): x-direction ot the incoming vector

        Returns:
            (tuple): (rz, rx): tuple comprising the z and x direction of the reflected neutron
        """
        vdotn = nz*zdir + nx*xdir
        weight = calc_supermirror_reflectivity(self.m, vdotn*2*V2Q_conic)
        if rn.random()<weight:
            rz, rx = zdir-2*vdotn*nz, xdir-2*vdotn*nx #classic reflected direction
            return rz, rx
        else:
            return zdir, xdir

    def rotate_vector(self, z, x, theta):
        """quick function to return a vector rotated around the origin by an angle theta

        Args:
            z (float): z-coordinate of the rotated vector
            x (float): x-coordinated of the rotated vector
            theta (float): the angle to rotate by

        Returns:
            (tuple): (z_rotated, x_rotated)
        """
        sin, cos = np.sin(theta), np.cos(theta)
        return z*cos - x*sin, z*sin + x*cos

    def rotate_neutron(self, z0, x0, zdir, xdir, theta_rot):
        """rotates a line defined by (z0, x0) + t*(zdir, xdir) around the origin by an angle of theta_rot

        Args:
            z0 (float): z0
            x0 (float): x0
            zdir (float): direction in z
            xdir (float): direction in x
            theta_rot (float): angle of rotation (rad)

        Returns:
            (tuple): (z0_r, x0_r, zdir_r, xdir_r): tuple representing the rotated line
        """
        sin, cos = np.sin(theta_rot), np.cos(theta_rot)
        rot_z0 = cos*z0 - sin*x0
        rot_x0 = sin*z0 + cos*x0
        rot_zdir = cos*zdir - sin*xdir
        rot_xdir = sin*zdir + cos*xdir
        return rot_z0, rot_x0, rot_zdir, rot_xdir

    def mirror_neutron(self, z0, x0, zdir, xdir, mirrored):
        """inverts neutron at the xy plane

        Args:
            z0 (float): 
            x0 (_type_): _description_
            zdir (_type_): _description_
            xdir (_type_): _description_
            theta_rot (_type_): _description_
        """
        return z0, x0*mirrored, zdir, xdir*mirrored

    def calc_interaction_time(self, z0, x0, vz, vx, theta_int):
        """calculates and returns the time till interaction of the neutron with the logspir returns None if no valid time is found

        Args:
            z0 (float): z-pos of neutron
            x0 (float): x-pos of neutron
            vz (float): z-vel of neutron
            vx (float): x-vel of neutron

        Returns:
            float: time until interaction None if none is found
        """
        if theta_int:
            z_int, _ = self.return_cart_coords(theta_int)
            t = (z_int - z0)/vz
            if t > 0.25*(self.zstart*self.phi_rot)/(vx**2+vz**2)**0.5:#need to make sure this is positive
                return t
        return None # else we dont return anything


    def return_timebranch_first_interaction(self, neutron: Neutron):
        """returns the index and time of the first first possible interaction of neutron and mirror

        Args:
            neutron (Neutron): incoming neutron

        Returns:
            (tuple): (branch, time, theta) a tuple containing the branch, time and theta
        """
        z0, x0, zdir, xdir = neutron.return_coords()
        #to calculate the intersection of the neutron with each branch of the polarizer, rotate the neutron accordingly
        min_time, branchind, theta_min = float('inf'), None, None
        for ind in range(self.branches*(self.double_sided+1)):
            if ind >= self.branches:
                indc = ind - self.branches
                mirrored = -1
                theta_rot = -(indc)*self.phi_rot
            else:
                indc = ind
                mirrored = 1
                theta_rot = -ind*self.phi_rot
            #print("what is this ", theta_rot, mirrored)
            zi, xi, z_diri, x_diri = self.mirror_neutron(z0, x0, zdir, xdir, mirrored)
            zi, xi, z_diri, x_diri = self.rotate_neutron(zi, xi, z_diri, x_diri, theta_rot)#instead of rotating the device clockwise, rotate the neutron in the opposite direction
            theta_int = self.return_precise_theata_int(zi, xi, z_diri, x_diri)
            t = self.calc_interaction_time(zi, xi, z_diri, x_diri, theta_int)
            if t:#only if the new time exists and is greater than 0
                if t < min_time:
                    min_time = t
                    branchind = indc
                    theta_min = theta_int
                    min_mirrored = mirrored
        #now if they exist return them all
        if branchind != None:
            return branchind, min_time, theta_min, min_mirrored
        else:
            return None

    def return_first_interaction(self, neutron: Neutron):
        """calculates the intersection with all spiral branches and returns the neutron after interacting with the first one reached
           if no interaction is possible returns None

        Args:
            neutron (Neutron): neutron

        Returns:
            (Neutron): the reflected Neutron
        """
        branchtime = self.return_timebranch_first_interaction(neutron)
        if not branchtime:
            return None
        #to calculate the intersection of the neutron with each branch of the polarizer,
        branch, time, theta, mirrored = branchtime
        theta_rot = self.phi_rot*branch
        z_int, x_int, vz, vx = neutron.prop_neutron(time) #neutron is propagated
        nz, nx = self.return_normal_vec(theta)

        nz, nx = self.rotate_vector(nz, nx, theta_rot)#normal vector returned has to be rotated according to the branch it originates from

        nz, nx, _, _ = self.mirror_neutron(nz, nx, 0, 0, mirrored)
        rz, rx = self.return_reflected_dir(nz, nx, vz, vx)
        return Neutron(z_int, x_int, rz, rx)# return the reflected neutron

    def propagate_neutron(self, neutron: Neutron, max_interaction=100):
        """given a single incoming neutron the function propagates the neutron through the mirror assembly and returns it and the intermediate postions for debugging

        Args:
            neutron (Neutron): incoming neutron to be propageted
            max_interaction (int, optional): max number of interactions after which the code returns the neutron (or should it abort?). Defaults to 100.

        Returns:
            (tuple): (list of intermediate positions (z, x, vz, vx), neutron)
        """
        positions = [neutron.return_coords()]#positions are only for testing and not needed for McStas
        for interaction in range(max_interaction):
            n_r = self.return_first_interaction(neutron)
            if not n_r:#if no more interactions can be reached the neutrons last init is returned
                return positions, neutron
            neutron = n_r
            positions.append(n_r.return_coords())
        return positions, neutron
