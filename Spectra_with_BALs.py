from astropy.io import fits
import numpy as np
import fitsio
from astropy.table import Table,vstack,join,Column
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib.pyplot as plt
import scipy as sp

plt.rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
plt.rc("text", usetex=True)

g_fontsize_legend = 9
g_fontsize_xlabel = 15
g_fontsize_ylabel = 15
g_fontsize_xyticklabels = 16

g_ylabel_flux = "Flux [10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]"
g_xlabel = "$\lambda$ [\AA]"

# Data from https://people.ast.cam.ac.uk/~rfc/atomdat.html & Harris et al. 2016
l_lya = 1215.67
l_nv = 1238.82
l_siv = 1393.76
l_civ = 1548.21
l_mgii = 2803

bal_tr = 'truth-16-1187_bal.fits'
bal_rr = 'rr_zbest-16-1187.fits'
mas_rr = 'rr_zbest-16-1187_mod.fits'
bal_sp = 'spectra-16-1187_bal.fits'
mas_sp = 'spectra-16-1187_mod.fits'

bal_tr1 = Table.read(bal_tr, hdu = 1)
bal_rr = Table.read(bal_rr, hdu = 1)
mas_rr = Table.read(mas_rr, hdu = 1)
bal_tr1.sort('TARGETID')
bal_rr.sort('TARGETID')
mas_rr.sort('TARGETID')

z_tr_b = bal_tr1['TRUEZ']
z_in_b = bal_tr1['Z_INPUT']
z_rr_b = bal_rr['Z']
z_rr_m = mas_rr['Z']

waveb_bal = fitsio.read(bal_sp,ext=2)
waver_bal = fitsio.read(bal_sp,ext=7)
fluxb_bal = fitsio.read(bal_sp,ext=3)
fluxr_bal = fitsio.read(bal_sp,ext=8)
waveb_mas = fitsio.read(mas_sp,ext=2)
waver_mas = fitsio.read(mas_sp,ext=7)
fluxb_mas = fitsio.read(mas_sp,ext=3)
fluxr_mas = fitsio.read(mas_sp,ext=8)

def main():
    #print(len(z_in_l), z_rr_l.shape, z_tr_l.shape, flux_b.shape)
    g_xlim = [3500,7800]
    g_ylim = [-1,10]
    
    for i in range(896):
        xcoords_civ = [(1.+z_tr_b[i])*l_civ,(1.+z_in_b[i])*l_civ,(1.+z_rr_b[i])*l_civ,(1.+z_rr_m[i])*l_civ]
        xcoords_siv = [(1.+z_tr_b[i])*l_siv,(1.+z_in_b[i])*l_siv,(1.+z_rr_b[i])*l_siv,(1.+z_rr_m[i])*l_siv]
        xcoords_lya = [(1.+z_tr_b[i])*l_lya, (1.+z_in_b[i])*l_lya,(1.+z_rr_b[i])*l_lya,(1.+z_rr_m[i])*l_lya]
        xcoords_nv = [(1.+z_tr_b[i])*l_nv,(1.+z_in_b[i])*l_nv,(1.+z_rr_b[i])*l_nv,(1.+z_rr_m[i])*l_nv]
        xcoords_mgii = [(1.+z_tr_b[i])*l_mgii,(1.+z_in_b[i])*l_mgii,(1.+z_rr_b[i])*l_mgii,(1.+z_rr_m[i])*l_mgii]
        
        colors = ['r','k','g', 'y']
        """

        #print(i)
        #print(z_tr_b[i], z_in_b[i], z_rr_b[i], z_rr_m[i])
        print(i, (1+z_tr_b[i])*l_lya, (1+z_in_b[i])*l_lya, (1+z_rr_b[i])*l_lya, (1+z_rr_m[i])*l_lya)
        """
        fig = plt.figure(figsize=(16, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(g_xlim)
        ax.set_ylim(g_ylim)
        ysmoothed_b = gaussian_filter1d(fluxb_bal[i,:], sigma=3)
        ysmoothed_r = gaussian_filter1d(fluxr_bal[i,:], sigma=3)
        plt.plot(waveb_bal,ysmoothed_b,color='darkblue',label='Blue channel')
        plt.plot(waver_bal,ysmoothed_r,color='blue',label='Red channel')
        #plt.plot(wavez,ysmoothed_z,color='deepskyblue',label='NIR channel')
        plt.ylabel(g_ylabel_flux, fontsize=g_fontsize_ylabel)
        plt.xlabel(g_xlabel, fontsize=g_fontsize_xlabel)

        for xc,c in zip(xcoords_civ,colors):
            plt.axvline(x=xc, c=c, linewidth = 1)    
        for xc,c in zip(xcoords_siv,colors):
            plt.axvline(x=xc, c=c, linewidth = 1)  
        for xc,c in zip(xcoords_lya,colors):
            plt.axvline(x=xc, c=c, linewidth = 1)  
        for xc,c in zip(xcoords_nv,colors):
            plt.axvline(x=xc, c=c, linewidth = 1)  
        for xc,c in zip(xcoords_mgii,colors):
            plt.axvline(x=xc, c=c, linewidth = 1)    

        plt.savefig("spectra_"+str(i)+".png")
        #plt.show()
        plt.close(fig)
        
if __name__ == "__main__":
    main()
    
