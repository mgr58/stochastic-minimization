"""
Particles are moving on a circular path with incommensurable frequencies.
If they start at time t = 0 at the same point, when and where on the path they
all meet again?

The method uses global optimization to find a minimum of the penalty function
that measures the angular spread of the particles.
"""

import numpy as np
from math import pi


def penalty(ttime, *parms):
    """
        The function weighting the angular separation of particles
        rotation on a circular path with different frequencies.
    """
    tmin = parms[-2]
    tmax = parms[-1]
    omega = parms[:-2]

    if ttime > tmax:
        pen = 100. + ttime - tmax
    elif ttime < tmin:
        pen = 100. + tmin - ttime
    else:
        pen = 0.
        for j in range(len(omega)-1):
            pen = pen - np.cos(((omega[j]-omega[j+1])*ttime)%(2*pi))

    return pen


from matplotlib.pyplot import show, figure


def show_graph(xmin, xmax, npt, fun, *pars):
    """
        Plot a graph of a function fun(x, *pars)
        for xmin < x < xmax)
    """
    step = (xmax - xmin)/(npt - 1)
    xxs = xmin + step*np.linspace(0, npt-1, npt)
    funvals = np.zeros(npt, dtype=np.float64)
    for j in range(npt):
        funvals[j] = fun(xxs[j], *pars)

    fig = figure()

    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(xxs, funvals)
    ax1.grid()

    show(block=False)

if __name__ == '__main__':

    from scipy import optimize
    import time

    PARAMS = (np.sqrt(2.), np.sqrt(3.), np.sqrt(7.), np.sqrt(11.))

    METHODS = ('anneal', 'basinhopping')
    TMIN = 10.
    TMAX = 500.
    PARPACK = PARAMS + (TMIN, TMAX)

    NPOINT = 400

    """
        Plot the graph of the penalty function
    """

    show_graph(TMIN, TMAX, NPOINT, penalty, *PARPACK)

    GUESS = (TMIN + TMAX)/2    # Initial guess.
    np.random.seed(1)

    for method in METHODS:

        print "\nmethod: %s" % method
        start_time = time.time()

        if method == 'anneal':
            res = optimize.anneal(penalty, GUESS, args=PARPACK,
                                  schedule='boltzmann',
                                  T0=10., full_output=True, maxiter=600,
                                  dwell=250, disp=True)
            timem = res[0]
            penmin = res[1]
        elif method == 'basinhopping':
            minimizer_kwargs = {"args": PARPACK}
            ret = optimize.basinhopping(penalty, 15., stepsize=30.,
                                        minimizer_kwargs=minimizer_kwargs,
                                        niter=500)
            timem = ret.x
            penmin = ret.fun
        else:
            print "Unknow algorithm requested"

        print "Elapsed time: %s seconds" % (time.time() - start_time)

        print("found global minimum: time = %.4f, penalty = %.6f"
              % (timem, penmin))

        for i in range(len(PARAMS)):
            angle = (timem*PARAMS[i])%(2*pi)
            print(" planet %2d:  angle rad = %10.4f  deg = %10.4f"
                  % (i, angle, angle/pi*180))

    show()      # release the graph
