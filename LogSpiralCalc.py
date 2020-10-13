import numpy as np
import matplotlib.pyplot as plt

class LogSpiral:
    """Allow to calculates parameters for a Log Spiral with two parameters
    """
    def __init__(self, xstart: float, xend: float, psi: float):
        """initializes the logarithmic mirror

        Args:
            xstart (float): where the logspiral cuts the x axis
            xend (float): at which x value the log spiral shall end, important for the linear approx
            psi (float): the angle at which the spiral cuts the xaxis
        """
        self.xstart = xstart
        self.theta_start = 0
        self.xend = xend
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
        return x,y

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
        for k in range(max_iterations):
            x_n = x_0 - self.function(x_0)/self.derivative(x_0)

            if abs(x_0-x_n) < precision:
                return x_n
            x_0 = x_n
        return None

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
        self.theta_end = self.return_precise_theta(precision=10**-7)
        xend,yend = self.return_cart_coords(self.theta_end)
        if plot:
            fig,ax = plt.subplots(1)
            ax.plot([xstart, xend], [ystart, yend], linestyle='-', marker=' ')
            theta_range = np.linspace(0, self.theta_end, 10)
            print(self.return_cart_coords(theta_range))
            ax.plot(*self.return_cart_coords(theta_range), linestyle='-', marker=' ')
            xlim = ax.get_xlim()
            ax.set_xlim(0,xlim[1])
            ax.set_aspect('equal')
        m = (yend-ystart)/(xend-xstart)
        y0 = yend - m*xend
        if plot:
            y = lambda x: x*m +y0
            ax.plot([0,xend],[y(0),y(xend)])
        return m,y0

class Line2D:
    def __init__(self, m, y0):
        self.m = m
        self.y0 = y0

class Intersection:
    """class to determine the intersection of a line and
    a logarithmic spiral
    """
    def __init__(self, logspir, line):
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
        k = self.logspir.k
        m = self.line.m
        xstart = self.logspir.xstart
        y0 = self.line.y0
        function = lambda theta: np.exp(k*theta)*(np.sin(theta) - m*np.cos(theta))-y0/xstart
        derivative = lambda theta: np.exp(k*theta)*\
        (np.cos(theta)*(1-k*m)+np.sin(theta)*(k+m))

        theta0 = self.return_prelim_theta_int()
        fig,ax = plt.subplots(1)
        theta_plot = np.linspace(0,theta0*2,100)
        ax.plot(theta_plot, function(theta_plot))
        for k in range(max_iterations):
            thetan = theta0 - function(theta0)/derivative(theta0)
            if abs(thetan-theta0) < precision:
                return thetan
            theta0 = thetan
        return None
