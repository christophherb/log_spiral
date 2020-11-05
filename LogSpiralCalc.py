import numpy as np
import matplotlib.pyplot as plt

class LogSpiral:
    """Allow to calculates parameters for a Log Spiral with two parameters
    """
    def __init__(self, xstart: float, xend: float, psi: float, theta_range: float):
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
        self.psi = psi
        self.psi_rad = np.pi/180*psi
        self.k = 1/np.tan(self.psi_rad)
        #function and derivative are used to find the angle theta at which the x_value of the spiral equals xend
        self.function = lambda theta: np.cos(theta)*self.xstart*np.exp(self.k*theta)-self.xend
        self.derivative = lambda theta: self.xstart*np.exp(self.k*theta)*(np.cos(theta)*self.k-np.sin(theta))
        self.theta_end = self.return_precise_theta()
        self.branches = int(theta_range/180*np.pi/self.theta_end) + 1

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
        """returns a close approx to the theta_end given a specific xend

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

    def return_equiv_vector(self, plot=False):
        """returns the line approximating the spiral

        Args:
            plot (bool, optional): Whether to plot the resulting equivalent line. Defaults to False.

        Returns:
            tuple: (x0, y0, xdir, ydir) 
        """
        xstart = self.xstart
        ystart = 0
        if self.theta_end == None:
            self.update_theta_end()
        xend, yend = self.return_cart_coords(self.theta_end)
        m = (yend-ystart)/(xend-xstart)
        #y0 = yend - m*xend
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
        return (xstart, ystart, 1, m)

    def plot_log_spir(self):
        """plots all log spirals in the branch and returns the figure and axis

        Returns:
            tuple: (fig, ax)
        """
        theta_rot = self.return_precise_theta()
        fig, ax = plt.subplots(1)
        theta_range = np.linspace(0, theta_rot, 101)
        for branch in range(self.branches):
            sin = np.sin(theta_rot*branch)
            cos = np.cos(theta_rot*branch)
            x, y = self.return_cart_coords(theta_range)
            x, y = cos*x - sin*y, sin*x + cos*y
            ax.plot(x, y, linestyle='-', marker=' ', color='black')
            ax.plot([0,x[-1]], [0, y[-1]], linestyle='-', marker=' ')
        return fig, ax

class Line2D:
    """quick class for 2D line objects
    """
    def __init__(self, x0, y0, xdir, ydir):
        """init line with m and y0 according to X = X0 + t*X_dir where large X is 2D vector

        Args:
            x0 (float): x value of point on the line
            y0 (float): y value of point on the line
            xdir (float): x direction of line
            ydir (float): y direction of line
        """
        self.x0 = x0
        self.y0 = y0
        self.xdir = xdir
        self.ydir = ydir

    def return_coords(self, t: float):
        """return x and y for a given t

        Args:
            t (float): t value at which to evaluate X0 + X_dir*t

        Returns:
            tuple: (x, y)
        """
        return self.x0 + self.xdir*t, self.y0 + self.ydir*t

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
        rot_x0 = cos*self.x0 - sin*self.y0
        rot_y0 = sin*self.x0 + cos*self.y0
        rot_xdir = cos*self.xdir - sin*self.ydir
        rot_ydir = sin*self.xdir + sin*self.ydir
        #rot_m = (sin + cos * m)/(cos - sin * m)
        #rot_y0 = self.y0*(cos + sin*rot_m)
        return rot_x0, rot_y0, rot_xdir, rot_ydir

    def rotate_line(self, theta: float):
        """rotates the line around the origin by an angle theta

        Args:
            theta (float): angle of rotation (deg)
        """
        self.x0, self.y0, self.xdir, self.ydir = self.return_rotated_line(theta)


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
        log_x0, log_y0, log_xdir, log_ydir = self.logspir.return_equiv_vector()
        line_x0, line_y0, line_xdir, line_ydir = self.line.x0, self.line.y0, self.line.xdir, self.line.ydir
        #find solutions to the equation line0 + linedir*t == log0 + logdir*lam and finally check if the x value makes sense
        try:
            t = (log_ydir*(line_x0-log_x0) + log_xdir*(log_y0-line_y0))/\
                (log_xdir*line_ydir - line_xdir*log_ydir)
            x_intersection = line_x0 + line_xdir*t
        except ZeroDivisionError:
            try:
                lam = (line_xdir*(line_y0-log_y0) + line_ydir*(log_x0-line_x0))/\
                    (log_ydir*line_xdir - log_xdir*line_ydir)
                x_intersection = log_x0 + log_xdir*lam
            except ZeroDivisionError:
                x_intersection = None
        if x_intersection and (self.logspir.xstart < x_intersection < self.logspir.xend):
            return x_intersection
        return None

    def return_prelim_theta_int(self):
        """returns the theta calculated from the line approximation
        to be used in a newton refinement

        Returns:
            float: theta_intersection (rad)
        """
        x_intersection = self.return_prelim_x()
        if x_intersection:
            y_intersection = self.line.ydir*(x_intersection - self.line.x0)/self.line.xdir + self.line.y0
            r_intersection = (x_intersection**2 + y_intersection**2)**0.5
            theta_prelim = np.log(r_intersection/self.logspir.xstart)\
                /self.logspir.k
            return theta_prelim
        return None

    def return_precise_thetaint(self, max_iterations=10, precision=10**-5,p=True):
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
        #transform linevector to m*x + y0
        m = self.line.ydir/self.line.xdir
        xstart = self.logspir.xstart
        y0 = self.line.y0+(0-self.line.x0)/self.line.xdir*self.line.ydir#intersection with y-axis

        def function(theta):#line intersects spiral where function(theta) = 0 
            return np.exp(k*theta)*(np.sin(theta) - m*np.cos(theta))-y0/xstart

        def derivative(theta):#theta derivative of the above function for newton iterations
            return np.exp(k*theta)*\
                    (np.cos(theta)*(1-k*m)+np.sin(theta)*(k+m))

        theta0 = self.return_prelim_theta_int()
        if theta0 is None:
            return None
        if p:
            _,ax = plt.subplots(1)
            theta_plot = np.linspace(0,self.logspir.theta_end*2,100000)
            ax.plot(theta_plot, function(theta_plot), linestyle='-', marker=' ', color='red')
            ax.plot([theta0], [function(theta0)], marker='s', color='orange')
            ax.plot([theta_plot[0], theta_plot[-1]],[0,0], marker=' ', color='green')
            ax.plot(theta_plot, derivative(theta0)*(theta_plot-theta0)+function(theta0), color='blue')
        for _ in range(max_iterations):
            thetan = theta0 - function(theta0)/derivative(theta0)
            if p:
                ax.plot(theta_plot, derivative(thetan)*(theta_plot-thetan)+function(thetan), linestyle='-', marker=' ', color='blue')
                ax.plot([thetan], [function(thetan)], marker='s', color='orange')
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
        if theta is None:
            return None
        return self.logspir.return_cart_coords(theta)

    def return_normal_vec(self, theta):
        """returns the normal vector of the spiral at given theta,
        relevant for the calculation of the reflection

        Args:
            theta (float): angle of intersection (deg)

        Returns:
            tuple: nx,ny x and y component of the normal vector
        """
        #print('theta', theta*180/np.pi)
        k = self.logspir.k
        norm_prefac = 1/(1 + k**2)**0.5
        nx = (np.cos(theta)+k*np.sin(theta))*norm_prefac
        ny = (np.sin(theta)-k*np.cos(theta))*norm_prefac
        return nx, ny

    def return_reflect_dir(self, theta, xdir, ydir):
        """returns the direction of a ray reflected at theta

        Args:
            theta (float): angle of intersection
            xdir (float): x component of incoming velocity vector
            ydir (float): y component of incoming velocity vector

        Returns:
            tuple: rx, ry; x and y component of reflected beam
        """
        nx, ny = self.return_normal_vec(theta)
        vtimesn = nx*xdir + ny*ydir
        #print('leftscalar', vtimesn/((nx*nx+ny*ny)*(vx*vx+vy*vy))**0.5)
        rx, ry = xdir-2*vtimesn*nx, ydir-2*vtimesn*ny
        #print('rightscalar', (rx*nx+ry*ny)/((nx*nx+ny*ny)*(rx*rx+ry*ry)**0.5))
        #print('after', rx**2 + ry**2)
        return rx, ry

    def return_all_branches(self):
        """tries to intersect all possible spirals with the neutron path

        Returns:
            list: list of 
        """
        logspir = self.logspir
        theta_rot = logspir.theta_end
        intersects = []
        fix_line = Line2D(self.line.x0, self.line.y0, self.line.xdir, self.line.ydir)
        hit = False
        for branch in range(logspir.branches):
            sin = np.sin(theta_rot*branch)
            cos = np.cos(theta_rot*branch)
            self.line.x0, self.line.y0, self.line.xdir, self.line.ydir \
                = fix_line.return_rotated_line(-branch*theta_rot*180/np.pi)
            print('rotated line', self.line.x0, self.line.y0, self.line.xdir, self.line.ydir)
            try:
                x_int, y_int = self.return_scatter_coords()
                x_int, y_int = cos*x_int - sin*y_int, sin*x_int + cos*y_int
                intersects.append((x_int, y_int))
                hit = True
            except TypeError:
                print('no can do ', branch)
                if hit==True:
                    break
        self.line = fix_line
        return intersects

    def plot_intersection(self):
        """plots the situation and returns the fig,ax
        """
        int_theta = self.return_precise_thetaint()
        if int_theta is None:
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

