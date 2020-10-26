import numpy as np
import matplotlib.pyplot as plt

class LogSpiral:
    """Allow to calculates parameters for a Log Spiral with two parameters
    """
    def __init__(self, xstart: float, xend: float, psi: float, branches: int):
        """initializes the logarithmic mirror

        Args:
            xstart (float): where the logspiral cuts the x axis
            xend (float): at which x value the log spiral shall end, important for the linear approx
            psi (float): the angle at which the spiral cuts the xaxis
            branches (int): how many identical branches should the log spiral comprise
        """
        self.xstart = xstart
        self.xend = xend
        self.theta_start = 0
        self.theta_end = None
        self.branches = branches
        self.psi = psi
        self.psi_rad = np.pi/180*psi
        self.k = 1/np.tan(self.psi_rad)
        self.function = lambda theta: np.cos(theta)*self.xstart*np.exp(self.k*theta)-self.xend
        self.derivative = lambda theta: self.xstart*np.exp(self.k*theta)*(np.cos(theta)*self.k-np.sin(theta))

    def return_r(self, theta):
        """returns the distance of the spiral at a given angle theta according to r(theta) = xstart*exp(k*theta_rad)

        Args:
            theta (float): angle (deg)

        Returns:
            float: distance to the origin (m)
        """
        return self.xstart*np.exp(self.k*theta)

    def return_cart_coords(self, theta):
        """returns cartesian coordinates x, y at a given theta value

        Args:
            theta (float): angle (deg)

        Returns:
            tuple: x, y (deg,deg)
        """
        r = self.return_r(theta)
        x = np.cos(theta)*r
        y = np.sin(theta)*r
        return x, y

    def return_prelim_theta(self):
        """returns a close approx to the theta_end

        Returns:
            float: theta_end approximation
        """
        prelim_theta = np.log(self.xend/self.xstart)*1/self.k
        return prelim_theta

    def return_precise_theta(self, max_iterations=10, precision=10**-5):
        """uses newton approx to calculate a precise value for theta_end and returns the value

        Args:
            max_iterations (int, optional): maximum number of iterations, if no convergence is reached return None. Defaults to 10.
            precision (float, optional): how precision of thetha end. Defaults to 10**-5.

        Returns:
            float: theta_end such that the x coord approximates xend
        """
        x_0 = self.return_prelim_theta()
        for _ in range(max_iterations):
            x_n = x_0 - self.function(x_0)/self.derivative(x_0)
            if abs(x_0-x_n) < precision:
                return x_n
            x_0 = x_n
        return None

    def return_theta_range(self):
        """return the range of theta

        Returns:
            tuple: theta_start, theta_end
        """
        if self.theta_end:
            return self.theta_start, self.theta_end
        self.update_theta_end()
        return self.theta_start, self.theta_end

    def update_theta_end(self, max_iterations=10, precision=10**-5):
        """updates the classes theta_end parameter

        Args:
            max_iterations (int, optional): maximum number of iterations, if no convergence is reached return None. Defaults to 10.
            precision ([type], optional): how precision of thetha end. Defaults to 10**-5.
        """
        self.theta_end = self.return_precise_theta(max_iterations=max_iterations, precision=precision)

    def return_equiv_line(self, plot=False):
        """returns an approximation of the mirror segment between theta_start and theta_end by a line
        """
        xstart = self.xstart
        ystart = 0
        if self.theta_end == None:
            self.update_theta_end()
        xend,yend = self.return_cart_coords(self.theta_end)
        m = (yend-ystart)/(xend-xstart)
        y0 = yend - m*xend
        if plot:
            _,ax = plt.subplots(1)
            ax.plot([xstart, xend], [ystart, yend], linestyle='-', marker=' ')
            theta_range = np.linspace(0, self.theta_end, 10)
            #print(self.return_cart_coords(theta_range))
            ax.plot(*self.return_cart_coords(theta_range), linestyle='-', marker=' ')
            xlim = ax.get_xlim()
            ax.set_xlim(0,xlim[1])
            ax.set_aspect('equal')
            y = lambda x: x*m +y0
            ax.plot([0,xend],[y(0),y(xend)])
        return m,y0

class Line2D:
    """quick class for 2D line objects
    """
    def __init__(self, m, y0):
        """init line with m and y0 according to y = y0 + m*x

        Args:
            m (float): slope of the line
            y0 (float): y intersection of the line
        """
        self.m = m
        self.y0 = y0

    def return_rotated_line(self, theta: float):
        """returns the slope and intersection of the turned line

        Args:
            theta (float): angle of turning (deg)

        Returns:
            tuple: (m, y0) tuple of slope and intersection of the turned line
        """
        theta *= np.pi/180
        sin = np.sin(theta)
        cos = np.cos(theta)
        m = self.m
        rot_m = (sin + cos * m)/(cos - sin * m)
        rot_y0 = self.y0*(cos + sin*rot_m)
        return rot_m, rot_y0

    def rotate_line(self, theta: float):
        """rotates the line around the origin by an angle theta

        Args:
            theta (float): angle of rotation (deg)
        """
        self.m, self.y0 = self.return_turned_line(theta)


class Intersection:
    """class to determine the intersection of a line and
    a logarithmic spiral
    """
    def __init__(self, logspir: LogSpiral, line: Line2D):
        """initialize the intersection with a line and logspir

        Args:
            logspir (LogSpiral): LogSpiral obeject to intersect with
            line (Line2D): A 2D line m*x + t
        """
        self.logspir = logspir
        self.line = line

    def return_prelim_x(self):
        """returns an approximation for the x value of the intersection
        between the line and the logarithmic spiral

        Returns:
            float: x_intersection
        """
        log_m, log_y0 = self.logspir.return_equiv_line()
        line_m, line_y0 = self.line.m, self.line.y0
        x_intersection = (log_y0-line_y0)/(line_m-log_m)
        return x_intersection

    def return_prelim_theta_int(self):
        """returns the theta calculated from the line approximation
        to be used in a newton refinement

        Returns:
            float: theta_intersection (rad)
        """
        x_intersection = self.return_prelim_x()
        y_intersection = self.line.m*x_intersection + self.line.y0
        r_intersection = (x_intersection**2 + y_intersection**2)**0.5
        theta_prelim = np.log(r_intersection/self.logspir.xstart)/self.logspir.k
        return theta_prelim

    def return_precise_thetaint(self, max_iterations=10, precision=10**-5,p=False):
        """returns a precise value for the intersection between the
        log spir and the ray using newton approx

        Args:
            max_iterations (int, optional): After how many tries
            do we give up. Defaults to 10.
            precision (float, optional): maximum difference between
            consecutive theta values. Defaults to 10**-5.
            p (bool, optional): plot the function which should be minimized. Defaults to False.

        Returns:
            float: theta value of intersection
        """
        k = self.logspir.k
        m = self.line.m
        xstart = self.logspir.xstart
        y0 = self.line.y0

        def function(theta):
            return np.exp(k*theta)*(np.sin(theta) - m*np.cos(theta))-y0/xstart

        def derivative(theta):
            return np.exp(k*theta)*\
                    (np.cos(theta)*(1-k*m)+np.sin(theta)*(k+m))

        theta0 = self.return_prelim_theta_int()
        if p:
            _,ax = plt.subplots(1)
            theta_plot = np.linspace(0,self.logspir.theta_end*2,100000)
            ax.plot(theta_plot, function(theta_plot), linestyle='-', marker=' ')
            ax.plot(theta0, function(theta0))
            ax.plot([theta_plot[0], theta_plot[-1]],[0,0], marker=' ')
        for _ in range(max_iterations):
            thetan = theta0 - function(theta0)/derivative(theta0)
            if p:
                ax.plot(theta_plot, derivative(theta0)*(theta_plot-theta0)+function(theta0), linestyle='-', marker=' ')
                ax.plot(thetan,function(thetan))
            if abs(thetan-theta0) < precision:
                return thetan
            theta0 = thetan
        return None

    def return_scatter_coords(self):
        """returns the carthesian coordinates of the intersection

        Returns:
            tuple: x_intersection, y_intersection
        """
        theta = self.return_precise_thetaint()
        return self.logspir.return_cart_coords(theta)

    def return_normal_vec(self,theta):
        """returns the normal vector of the spiral at given theta,
        relevant for the calculation of the reflection

        Args:
            theta (float): angle of intersection (deg)

        Returns:
            tuple: nx,ny x and y component of the normal vector
        """
        #print('theta', theta*180/np.pi)
        k = self.logspir.k
        prefac = 1/(1 + k**2)**0.5
        nx = (np.cos(theta)+k*np.sin(theta))*prefac
        ny = (np.sin(theta)-k*np.cos(theta))*prefac
        #print('normal', nx, ny)
        return nx, ny

    def return_reflect_dir(self, theta, vx, vy):
        """returns the direction of a ray reflected at theta

        Args:
            theta (float): angle of intersection
            vx (float): x component of incoming velocity vector
            vy (float): y component of incoming velocity vector

        Returns:
            tuple: 1, ry/rx; x and y component of reflected beam normalized to x = 1
        """
        #print('before', vx**2 + vy**2)
        nx, ny = self.return_normal_vec(theta)
        #print('normalize', nx**2 + ny**2)
        vtimesn = nx*vx + ny*vy
        #print('leftscalar', vtimesn/((nx*nx+ny*ny)*(vx*vx+vy*vy))**0.5)
        rx, ry = vx-2*vtimesn*nx, vy-2*vtimesn*ny
        #print('rightscalar', (rx*nx+ry*ny)/((nx*nx+ny*ny)*(rx*rx+ry*ry)**0.5))
        #print('after', rx**2 + ry**2)
        return 1, ry/rx

    def plot_intersection(self):
        """plots the situation and returns the fig,ax
        """
        int_theta = self.return_precise_thetaint()
        if int_theta < 0 or int_theta > self.logspir.theta_end:
            #print('shitty no')
            fig,ax = plt.subplots(1)
            theta_range = self.logspir.return_theta_range()
            theta_range = np.linspace(theta_range[0],theta_range[1],100)
            ax.plot(*self.logspir.return_cart_coords(theta_range),\
             linestyle='-', marker=' ', color='black', label='mirror')
            ax.plot([0, self.logspir.xend],[self.line.y0, self.line.y0 + self.logspir.xend*self.line.m], linestyle='-', marker=' ',\
            color='darkgreen', label='incoming')
            return fig,ax
        fig,ax = plt.subplots(1)
        x_int, y_int = self.logspir.return_cart_coords(int_theta)
        theta_range = self.logspir.return_theta_range()
        theta_range = np.linspace(theta_range[0],theta_range[1],100)
        ax.plot(*self.logspir.return_cart_coords(theta_range),\
         linestyle='-', marker=' ', color='black', label='mirror')
        ax.plot([0, x_int],[self.line.y0, y_int], linestyle='-', marker=' ',\
        color='darkgreen', label='incoming')
        x_back = x_int*2
        vx = 1
        vy = self.line.m
        #print('incoming', vx, vy)
        #show the normal as well for fixing kot
        nx, ny = self.return_normal_vec(int_theta)
        ny /= nx
        nx = 1
        ax.plot([x_int, x_int-0.5], [y_int, y_int -0.5*ny],\
        linestyle='-', marker=' ')
        rx, ry = self.return_reflect_dir(theta=int_theta, vx=vx, vy=vy)
        #print('reflected', rx, ry)
        y_back = y_int + (x_back-x_int)*ry
        ax.plot([x_int, x_back], [y_int, y_back], color='red', label='reflected')
        return fig, ax

