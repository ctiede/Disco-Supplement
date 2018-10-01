import sys
import h5py as h5
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import discoUtil as du

time_avg = []

def plot_density_avg( R ):

    fig, ax = plt.subplots( 1 )

    #Average each annulus over number of checkpoints
    Sig = sum( time_avg )/float( len( time_avg ) )
    np.savetxt('sig_time_azi_avg.txt', Sig )

    #ax.plot( R, Sig, color='b' )
    #ax.set_xlabel(r'$r$')
    #ax.set_ylabel(r'$\Sigma$')
    #ax.set_xlim([0,6])

    #plotname = 'plot_R_avg.png'
    #print("   Saving {0:s}...".format(plotname))
    #fig.savefig(plotname)


def density_azi_avg( file ):

    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = du.loadCheckpoint(file)
    pars = du.loadPars(file)
    phi_max = pars['Phi_Max']

    R   = np.unique(r)
    Sig = np.zeros( len(R) )
    piph = dat[3]

    for i in range( len(R) ):
        #Take azimuthal average for each annulus
        ind = (r==R[i])
        dphi = get_anulus_dphi( piph[ind], phi_max )
        #Sig[i] = np.sum( prim[np.where(ind),0] )/np.float( len(r[ind]) )
        rho = prim[np.where(ind),0]
        Sig[i] = np.sum( rho*dphi/phi_max )

    #append to list of 1d-arrays holding azimuthal averaged data
    time_avg.append(Sig)

    return R

def get_anulus_dphi( phi_ind, phi_max ):

    imax = np.argmax( phi_ind )
    apiph = np.roll( phi_ind, -imax-1 )
    phif = np.empty( len(apiph)+1 )
    phif[:-1] = apiph[:]
    phif[-1]  = apiph[0] + phi_max

    return np.diff( phif )



if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    for f in files:
        R = density_azi_avg( f )

    plot_density_avg( R )
    np.savetxt('R.txt', R )
