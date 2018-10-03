import os
import sys
import argparse as ag
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import discoUtil as du

def get_r_bh( R, phif, om=None ):

    phi = 0.5*( phif[0:-1] + phif[1:] )        #Avg phifs
    x = R*np.cos(phi)
    y = R*np.sin(phi)

    #a = 1.0
    #q = 1.0
    #mu = q/(1.+q)
    #bh_p = [a*mu, 0]  #Cartesion
    #bh_m = [-a*(1-mu), 0]
    bh_p = [ 0.5, 0]
    bh_m = [-0.5, 0]

    rbh_p = np.sqrt( ( bh_p[0] - x )**2 + ( bh_p[1] - y )**2 )
    rbh_m = np.sqrt( ( bh_m[0] - x )**2 + ( bh_m[1] - y )**2 )

    return rbh_p, rbh_m


def get_theta( r, r_bh, q=1.0, which='plus' ):

    a = 1.0
    q = 1.0
    mu = q/(1.+q)

    if which=='plus':
        a05 = a*mu
    elif which=='minus':
        a05 = a*(1-mu)
    else:
        print "Get theta 'which' must be 'plus' or 'minus' "

    costh = ( r*r + r_bh*r_bh - a05*a05 )/( 2*r*r_bh )

    return np.arccos( costh )


def bhPlot(fig, ax, rjph, piph, r, prim, label, pars, vmin=None, vmax=None,
            noGhost=False, colorbar=True, xlabel=None, ylabel=None, log=False,
            rmax=None, planets=None):


    #Always in binary rotating frame
    w = 1.0
    rho = prim[:,0]
    vr  = prim[:,2]
    vp  = (prim[:,3] - w)*r
    #vz  = prim[:,4]

    q_r = rho*vr
    q_p = rho*vp
    q_mag = np.sqrt( q_r*q_r + q_p*q_p )

    phi_max = pars['Phi_Max']

    N = 50
    dphi = phi_max/N
    slice_tot = np.zeros(N)
    slice_ct  = np.zeros(N)

    try:
        #cmap = mpl.cm.inferno
        #cmap = mpl.cm.RdYlBu_r
        cmap = mpl.cm.coolwarm_r
    except:
        cmap = mpl.cm.afmhot

    if vmin is None:
        vmin = q_r.min()
    if vmax is None:
        vmax = q_r.max()

    Rs = np.unique(r)
    if noGhost:
        Rs = Rs[:-2]

    if rmax is None or rmax <=0.0:
        lim_float = True
    else:
        lim_float = False

    if log:
        norm = mpl.colors.SymLogNorm( vmin=vmin, vmax=vmax, linthresh=0.1, linscale=0.5 )
    else:
        norm = mpl.colors.Normalize( vmin, vmax )

    for i, R in enumerate(Rs):
        ind = r==R
        imax = np.argmax(piph[ind])

        apiph = np.roll(piph[ind], -imax-1)
        #aq = np.roll(q[ind], -imax-1)
        #rho_a = np.roll( rho[ind], -imax-1 )
        aq_r  = np.roll( q_r[ind] , -imax-1 )
        aq_p  = np.roll( q_p[ind] , -imax-1 )

        phif = np.empty(len(apiph)+1)
        phif[1:] = apiph[:]
        phif[0] = apiph[-1] - phi_max

        rf = rjph[i:i+2]

        x = rf[:,None] * np.cos(phif)[None,:]
        y = rf[:,None] * np.sin(phif)[None,:]
        #Move to BH_p frame
        x -= 0.5

    #=================== Do rotation into bh frame ==========================

        mu = 0.5
        which = 'plus'
        rbh_p, rbh_m = get_r_bh( R, phif )
        if which=='plus':
            r_bh = rbh_p
        elif which=='minus':
            r_bh = rbh_m
        theta = get_theta( R, r_bh, q=1.0, which=which )
        bq = aq_r*np.cos( theta ) - aq_p*np.sin( theta ) #Rotate polar coords

    #=======================================================================

        C = ax[0].pcolormesh(x, y, bq[None,:],
                cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)

        if lim_float and rf.max() > rmax:
            rmax = rf.max()

    #=================== Make azimutal plot ================================
        # TODO: a list or array with N bins such that dphi = 2pi/N
        #       and a list/array of number of elements included
        phi = 0.5*( phif[0:-1] + phif[1:] )        #Avg phifs
        x = R*np.cos(phi) - 0.5
        y = R*np.sin(phi)

        r_cut = 0.2
        cut = ( r_bh<r_cut )
        phi_bh = np.arctan2( y[cut], x[cut] )
        #figure out which bin to average with
        ns = np.floor_divide( phi_bh, dphi )
        #print phif[mask]
        #print phi[mask], ns
        for i, n in enumerate( ns ):
            #add to that bin + add to counter
            n = int(n)
            slice_tot[n] += bq[cut][i]
            slice_ct[n]  += 1
    #at end, average and plot
    slice_ct[slice_ct==0] = 1  #slice_tot here should still be 0
    slice_avg = slice_tot/slice_ct

    eps = 0.2
    Ns = dphi*np.arange(N)

    #TODO: Make this do the same thing but with a bar chart instead of lines
    #      since each point represents a bin
    points = np.array([ Ns, slice_avg ]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection( segments, cmap=cmap, norm=norm )
    lc.set_array( slice_avg )
    lc.set_linewidth( 4 )
    ax[1].add_collection(lc)
    ax[1].set_xlim(0.0, 2*np.pi )
    ax[1].set_ylim( np.min(slice_avg)-eps, np.max(slice_avg)+eps )
    ax[1].set_xlabel( r'$\phi\,(radians)$', fontsize=18 )
    ax[1].set_ylabel( r'$\overline{ \Sigma\,v_r }$' , fontsize=24)


    if planets is not None:
        rpl = planets[:,3]
        ppl = planets[:,4]
        xpl = rpl * np.cos(ppl)
        ypl = rpl * np.sin(ppl)
        if which=='plus':
            xpl -= 0.5
        elif which=='minus':
            xpl += 0.5
        ax[0].plot(xpl, ypl, color='black', ls='', marker='o', mew=0, ms=5)

    ax[0].set_aspect('equal')
    ax[0].set_xlim(-rmax, rmax)
    ax[0].set_ylim(-rmax, rmax)
    ax[0].plot( -0.5, 0.0, marker='x', ms=7, mew=1.5, color='grey' )

    if colorbar:
        cb = fig.colorbar(C, ax=ax[1] )
        cb.set_label(label, fontsize=24)

    if xlabel == None:
        xlabel = r'$x$'
    if ylabel == None:
        ylabel = r'$y$'

    ax[0].set_xlabel(xlabel, fontsize=18)
    ax[0].set_ylabel(ylabel, fontsize=18)


def accrete_profile( fig, ax, r, prim ):

    w = 1.0  #Always in binary rotating frame
    rho = prim[:,0]
    vr  = prim[:,2]
    vp  = (prim[:,3] - w)*r




def plotCheckpoint(file, log=None, noGhost=False, om=None, bounds=None,
                    rmax=None, planets=False, k=None):

    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = du.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]
    pars = du.loadPars(file)

    if k is None:
        k = int(zkph.shape[0]/2-1)

    zind = (z>zkph[k]) * (z<zkph[k+1])
    r = r[zind]
    phi = phi[zind]
    prim = prim[zind,:]
    primPhi0 = primPhi0[k,:]
    piph = piph[zind]


    if planets:
        planetDat = dat[4]
    else:
        planetDat = None

    if om is not None:
        phi1 = phi - om*t
        piph1 = piph - om*t
        if planetDat is not None:
            planetDat[:,4] -= om*t
    else:
        phi1 = phi
        piph1 = piph


    varnames, vartex, num_c, num_n = du.getVarNames(file)
    #nq = prim.shape[1]
    nq = num_c + num_n

    Zs = np.unique(z)
    z_eq = Zs[len(Zs)/2]
    eq_ind = (z==z_eq)
    title = "DISCO t = {0:.1f}".format(t)
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

#    if vars is None:
#        vars = range(nq)
#    if logvars is None:
#        logvars = []

#============ Make Plot--for now just \Sigma*v_r ===============

    fig, ax = plt.subplots( 1, 2, figsize=(18,10) )

    label = r"$\Sigma\,v_r$"
    vecname = 'accrt'

    print("   Plotting...")

    if bounds is not None:
        vmin, vmax = bounds[q]
    else:
        vmin, vmax = None, None


    bhPlot( fig, ax, rjph, piph1, r, prim, label, pars, log=log,
            vmin=vmin, vmax=vmax, rmax=rmax, planets=planetDat )
    #accrete_profile( fig, ax[1], r, piph1, prim )

    fig.suptitle(title, fontsize=24)

    if log:
        plotname = "plot_bh_{0:s}_log_{1:s}.png".format( name, vecname )
    else:
        plotname = "plot_bh_{0:s}_lin_{1:s}.png".format( name, vecname )

    print("   Saving {0:s}...".format(plotname))
    fig.savefig(plotname)
    plt.close(fig)

#==========================================================================

def calcBounds(files):

    f = files[0]
    t, r, phi, z, prim, dat = du.loadCheckpoint(f)

    num_q = prim.shape[1]
    bounds = [[prim[:,q].min(), prim[:,q].max()] for q in range(num_q)]
    bounds = np.array(bounds)

    for f in files[1:]:
        t, r, phi, z, prim, dat = du.loadCheckpoint(f)
        for q in range(num_q):
            bounds[q,0] = min(bounds[q,0], prim[:,q].min())
            bounds[q,1] = max(bounds[q,1], prim[:,q].max())

    return bounds

def writeBoundsFile(filename, names, bounds):

    lines = ["{0:s} {1:.12g} {2:.12g}\n".format(names[q], bounds[q,0], bounds[q,1])
                for q in range(len(names))]

    f = open(filename, "w")
    for line in lines:
        f.write(line)

    f.close()

def readBoundsFile(filename, num_q):

    bounds = np.loadtxt(filename, usecols=[1,2])

    bounds = bounds[:num_q]

    return bounds

def getBounds(use_bounds, names, files):

    num_q = len(names)

    if use_bounds is not None:
        if use_bounds is True:
            bounds = calcBounds(files)
        else:
            if os.path.isfile(use_bounds):
                bounds = readBoundsFile(use_bounds, num_q)
            else:
                bounds = calcBounds(files)
                writeBoundsFile(use_bounds, names, bounds)
    else:
        bounds = None

    return bounds

if __name__ == "__main__":

    parser = ag.ArgumentParser(description="Create 2D plots of Disco variables.")
    parser.add_argument('checkpoints', nargs='+',
                            help="Checkpoint (.h5) files to plot.")
#    parser.add_argument('-v', '--vars', nargs='+', type=int,
#                            help="Variables to plot.")
#    parser.add_argument('-l', '--logvars', nargs='+', type=int,
#                            help="Variables to plot logscale.")
    parser.add_argument( '-l', '--logscale', action='store_true',
                            help='Set color scale to log')
    parser.add_argument('-p', '--planets', action='store_true',
                            help="Plot planets.")
    parser.add_argument('-b', '--bounds', nargs='?', const=True,
                            help="Use global max/min for bounds. Optional argument BOUNDS is a file. If it exists, it will be read for parameter bounds. If it does not exist the global max/min will be calculated and saved to the file.")
    parser.add_argument('-r', '--rmax', type=float,
                            help="Set plot limits to RMAX.")
    parser.add_argument('-o', '--omega', type=float,
                            help="Rotate frame at rate OMEGA.")
    parser.add_argument('--noghost', action='store_true',
                            help="Do not plot ghost zones.")

    args = parser.parse_args()

#    vars = args.vars
#    logvars = args.logvars
    log = args.logscale
    om  = args.omega
    rmax = args.rmax
    use_bounds = args.bounds
    planets = args.planets
    noghost = args.noghost

    files = args.checkpoints

    names, texnames, num_c, num_n = du.getVarNames(files[0])

    bounds = getBounds(use_bounds, names, files)

    if om is None:
        om = 1.0
    if rmax is None:
        rmax = 0.3

    for f in files:
        plotCheckpoint(f, log=log, bounds=bounds, om=om,
                        rmax=rmax, noGhost=noghost, planets=planets)
