from astropy.io import fits
import numpy as np
import fitsio
from astropy.table import Table,vstack,join,Column
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib.pyplot as plt
import scipy as sp

plt.rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
plt.rc("text", usetex=True)

g_fontsize_legend = 10
g_fontsize_xlabel = 15
g_fontsize_ylabel = 10
g_fontsize_ylabel_2 = 15
g_fontsize_xyticklabels = 16

g_ylabel_flux = "Flux [10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]"
g_ylabel_ivar = "$\sigma^{-2}$"
g_xlabel = "$\lambda$ [\AA]"
g_xlim = [3500,5500]

spec_1 = "spectra_bal.fits"
spec_2 = "spectra_mod.fits"
waveb1 = fitsio.read(spec_1,ext=2)
waver1 = fitsio.read(spec_1,ext=7)
flux_b1 = fitsio.read(spec_1,ext=3)
ivar_b1 = fitsio.read(spec_1,ext=4)
flux_r1 = fitsio.read(spec_1,ext=8)
ivar_r1 = fitsio.read(spec_1,ext=9)
waveb2 = fitsio.read(spec_2,ext=2)
waver2 = fitsio.read(spec_2,ext=7)
flux_b2 = fitsio.read(spec_2,ext=3)
ivar_b2 = fitsio.read(spec_2,ext=4)
flux_r2 = fitsio.read(spec_2,ext=8)
ivar_r2 = fitsio.read(spec_2,ext=9)

truth = "truth-file.fits"
tru_2 = Table.read(truth, hdu = 2)
tru_3 = Table.read(truth, hdu = 3)
tarID_q = tru_2['TARGETID'].data #quasar
tarID_b = tru_3['TARGETID'].data #bal
z = tru_3['Z'].data
AI_CIV =tru_3['AI_CIV'].data
NCIV = tru_3['NCIV_450'].data 
vmin = tru_3['VMIN_CIV_450'].data
vmax = tru_3['VMAX_CIV_450'].data
baltem = tru_3['BAL_TEMPLATEID'].data

#fbal = 'bal_templates_v3.0.fits'
#tru_1 = Table.read(fbal, hdu = 1)
#temp = tru_1['TEMP'].data

def main():
    
    #i = 133
    k = 46
    """
    w0 = 944.6
    array = []
    for j in range(2474):
        val = w0 + (j)*(0.3)
        array.append(val)
    wave_out = np.array(array)
    print(wave_out)
    wave = (1+z[int(i)])*wave_out
    """
    
    fig, axarr = plt.subplots(2, sharex=True)
    axarr[0].set_title("Mock DESI Y1 at $z = $ 1.808")#+ str(z[int(k)]))
    axarr[0].plot(waveb1, ivar_b1[k,:],color='red')
    #axarr[0].plot(waver1, ivar_r1[k,:],color='red')
    axarr[0].plot(waveb2, ivar_b2[k,:],color='darkblue')
    #axarr[0].plot(waver2, ivar_r2[k,:],color='blue')
    #axarr[0].vlines(4330, -2, 20, color='grey', linestyles='--')
    #axarr[0].vlines(4450, -2, 20, color='grey', linestyles='--')
    #axarr[0].vlines(6345, -2, 20, color='grey', linestyles='--')
    #axarr[0].vlines(6370, -2, 20, color='grey', linestyles='--') 
    axarr[0].set_ylabel(g_ylabel_ivar, fontsize=g_fontsize_ylabel_2)
    axarr[0].set_ylim(-1,5)
    #axarr[1].plot(waveb1, flux_b1[int(i),:],color='darkblue')
    ysmoothed_b = gaussian_filter1d(flux_b1[k,:], sigma=3)
    ysmoothed_r = gaussian_filter1d(flux_r1[k,:], sigma=3)
    #axarr[1].plot(waver1,flux_r1[int(i),:],color='darkblue')
    #axarr[1].vlines(4330, -2, 20, color='grey', linestyles='--')
    #axarr[1].vlines(4450, -2, 20, color='grey', linestyles='--')
    #axarr[1].vlines(6345, -2, 20, color='grey', linestyles='--')
    #axarr[1].vlines(6370, -2, 20, color='grey', linestyles='--')
    axarr[1].plot(waveb1,ysmoothed_b,color='darkblue')
    axarr[1].plot(waver1,ysmoothed_r,color='blue')
    #axarr[1].plot(waveb2, flux_b2[int(i),:],color='blue')
    #axarr[1].plot(waver2, flux_r2[int(i),:],color='darkblue')
    #axarr[1].set_xlabel(g_xlabel, fontsize=g_fontsize_xlabel)
    axarr[1].set_ylabel(g_ylabel_flux, fontsize=g_fontsize_ylabel)
    axarr[1].set_ylim(0,10)
    axarr[1].set_xlabel(g_xlabel, fontsize=g_fontsize_xlabel)
    axarr[1].set_xlim(g_xlim)
    legend_lines = []
    legend_labels = []
    legend_lines.append(plt.Line2D((0, 0), (0, 0), color='red', linestyle='-', linewidth=2.0, markeredgecolor='red'))
    legend_lines.append(plt.Line2D((0, 0), (0, 0), color='darkblue', linestyle='-', linewidth=2.0, markeredgecolor='darkblue'))
    #legend_lines.append(plt.Line2D((0, 0), (0, 0), color='blue', linestyle='-', linewidth=2.0, markeredgecolor='blue'))
    legend_labels.append("Original spectrum - B band")
    legend_labels.append("Masked BAL - B band")
    #legend_labels.append("Masked BAL - Red channel")
    """
    legend_lines.append(plt.Line2D((0, 0), (0, 0), color='darkred', linestyle='-', linewidth=1.0, markeredgecolor='darkred'))
    legend_lines.append(plt.Line2D((0, 0), (0, 0), color='darkblue', linestyle='-', linewidth=1.0, markeredgecolor='darkblue'))
    legend_labels.append("BAL-QSO bands R")
    legend_labels.append("Masked BAL bands R")
    """
    axarr[0].legend(legend_lines, legend_labels, loc="upper left", numpoints=1, ncol=1, fontsize=g_fontsize_legend, handlelength=1.0,frameon=False)
    plt.savefig("ivar_flux_indx_"+str(k)+".pdf")
    plt.show(fig)
    plt.close(fig)

if __name__ == "__main__":
    main()
    
