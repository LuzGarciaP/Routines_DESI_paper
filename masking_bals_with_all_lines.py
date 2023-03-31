import numpy as np
import fitsio
from astropy.io import fits
from astropy.table import Table,vstack,join,Column
import scipy as sp

fbal = "truth-files.fits"
spectra =fits.open("spectra-file.fits")

tru_1 = Table.read(fbal, hdu = 1)
tru_2 = Table.read(fbal, hdu = 2)
tru_3 = Table.read(fbal, hdu = 3)

tarID_q = tru_2['TARGETID'].data #quasar
tarID_b = tru_3['TARGETID'].data #bal
z = tru_3['Z'].data
NCIV = tru_3['NCIV_450'].data ## N_CIV is the actual number of BAL troughs in each pixel
ai =tru_3['AI_CIV'].data
vmin = tru_3['VMIN_CIV_450'].data
vmax = tru_3['VMAX_CIV_450'].data
baltem = tru_3['BAL_TEMPLATEID'].data

# Implement the masks on the BAL features:
wave_b=spectra['B_WAVELENGTH'].data
wave_r=spectra['R_WAVELENGTH'].data
wave_z=spectra['Z_WAVELENGTH'].data
ivar_b=spectra['B_IVAR'].data 
ivar_r=spectra['R_IVAR'].data
ivar_z=spectra['Z_IVAR'].data

lmin_B = wave_b[0]
lmax_B = wave_b[2379]
lmin_R = wave_r[0]
lmax_R = wave_r[2115]
lmin_Z = wave_z[0]
lmax_Z = wave_z[2115]

#for i,targetid in enumerate(tarID_b[0:58]):
for i,targetid in enumerate(tarID_b):  
    indx = sp.where(targetid==tarID_q)
    #print(int(indx))
    #ls = sp.constants.c*10**-3 ##Speed of light in km/s
    ls = 3E5 #Speed of light in km/s
        
    lLya = 1216 #\AA
    lNV = 1241 #\AA
    lSiIV=1394 #\AA
    lCIV=1549 #\AA
                  
    # Calculate mask width for each absorption line (ref and obs frames)
    # Maximum number of BAL troughs that can be generated = 27

    # Wavelengths associated to Lya
    lMin_lya = lLya*(1-vmin/ls)
    lMin_obs_lya = (1.+z[i])*lMin_lya
    lMax_lya = lLya*(1-vmax/ls)
    lMax_obs_lya = (1.+z[i])*lMax_lya

    # Wavelengths associated to NV
    lMin_NV = lNV*(1-vmin/ls)
    lMin_obs_NV = (1.+z[i])*lMin_NV
    lMax_NV = lNV*(1-vmax/ls)
    lMax_obs_NV = (1.+z[i])*lMax_NV

    # Wavelengths associated to SiIV
    lMin_SiIV = lSiIV*(1-vmin/ls)
    lMin_obs_SiIV = (1.+z[i])*lMin_SiIV
    lMax_SiIV = lSiIV*(1-vmax/ls)
    lMax_obs_SiIV = (1.+z[i])*lMax_SiIV

    # Wavelengths associated to CIV
    lMin = lCIV*(1-vmin/ls)
    lMin_obs = (1.+z[i])*lMin
    lMax = lCIV*(1-vmax/ls)
    lMax_obs = (1.+z[i])*lMax

    num_bal = int(NCIV[i])
    #print(num_bal)
       
    for j in range(num_bal):
        lmin_tmp = float(lMin_obs[i,j]) #CIV
        lmax_tmp = float(lMax_obs[i,j])
        lmin_tmp1 = float(lMin_obs_lya[i,j]) #Lya
        lmax_tmp1 = float(lMax_obs_lya[i,j])
        lmin_tmp2 = float(lMin_obs_NV[i,j]) #NV
        lmax_tmp2 = float(lMax_obs_NV[i,j])
        lmin_tmp3 = float(lMin_obs_SiIV[i,j]) #SiIV
        lmax_tmp3 = float(lMax_obs_SiIV[i,j])
        # Create the wavelength mask and set it in the corresponding ivar 
        #Mask only inside the region of filter B:
        if(lmax_tmp < lmin_R).any():
            w = np.where((wave_b > lmin_tmp) & (wave_b < lmax_tmp))
            #print('Mask only goes in B filter')
            ivar_b[indx,w] = 0

        if(lmax_tmp1 < lmin_R).any():
            w = np.where((wave_b > lmin_tmp1) & (wave_b < lmax_tmp1))
            ivar_b[indx,w] = 0
        if(lmax_tmp2 < lmin_R).any():
            w = np.where((wave_b > lmin_tmp2) & (wave_b < lmax_tmp2))
            ivar_b[indx,w] = 0    
        if(lmax_tmp3 < lmin_R).any():
            w = np.where((wave_b > lmin_tmp3) & (wave_b < lmax_tmp3))
            ivar_b[indx,w] = 0    
            
            
        #Mask both filters B and R:    
        elif(lmax_tmp < lmax_B).any():
            w1 = np.where((wave_b > lmin_tmp) & (wave_b < lmax_tmp))
            w2 = np.where((wave_r > lmin_tmp) & (wave_r < lmax_tmp))
            ivar_b[indx,w1] = 0
            ivar_r[indx,w2] = 0

        elif(lmax_tmp1 < lmax_B).any():
            w1 = np.where((wave_b > lmin_tmp1) & (wave_b < lmax_tmp1))
            w2 = np.where((wave_r > lmin_tmp1) & (wave_r < lmax_tmp1))
            ivar_b[indx,w1] = 0
            ivar_r[indx,w2] = 0            
        elif(lmax_tmp2 < lmax_B).any():
            w1 = np.where((wave_b > lmin_tmp2) & (wave_b < lmax_tmp2))
            w2 = np.where((wave_r > lmin_tmp2) & (wave_r < lmax_tmp2))
            ivar_b[indx,w1] = 0
            ivar_r[indx,w2] = 0  
        elif(lmax_tmp3 < lmax_B).any():
            w1 = np.where((wave_b > lmin_tmp3) & (wave_b < lmax_tmp3))
            w2 = np.where((wave_r > lmin_tmp3) & (wave_r < lmax_tmp3))
            ivar_b[indx,w1] = 0
            ivar_r[indx,w2] = 0  
            
        #Region around R
        elif(lmax_tmp < lmin_Z).any():
            # Mask falls only in the R region:
            if(lmin_tmp > lmax_B):
                #print('Mask only goes in R filter')
                w = np.where((wave_r > lmin_tmp) & (wave_r < lmax_tmp))
                ivar_r[indx,w] = 0
            # Mask falls both in B and R regions:     
            else:
                #print('Mask both B and R filters')
                w1 = np.where((wave_b > lmin_tmp) & (wave_b < lmax_tmp))
                w2 = np.where((wave_r > lmin_tmp) & (wave_r < lmax_tmp))
                ivar_b[indx,w1] = 0                
                ivar_r[indx,w2] = 0
                
        elif(lmax_tmp1 < lmin_Z).any():
            # Mask falls only in the R region:
            if(lmin_tmp1 > lmax_B):
                #print('Mask only goes in R filter')
                w = np.where((wave_r > lmin_tmp1) & (wave_r < lmax_tmp1))
                ivar_r[indx,w] = 0                                
            # Mask falls both in B and R regions:     
            else:
                #print('Mask both B and R filters')
                w1 = np.where((wave_b > lmin_tmp1) & (wave_b < lmax_tmp1))
                w2 = np.where((wave_r > lmin_tmp1) & (wave_r < lmax_tmp1))
                ivar_b[indx,w1] = 0                
                ivar_r[indx,w2] = 0

        elif(lmax_tmp2 < lmin_Z).any():
            # Mask falls only in the R region:
            if(lmin_tmp2 > lmax_B):
                #print('Mask only goes in R filter')
                w = np.where((wave_r > lmin_tmp2) & (wave_r < lmax_tmp2))
                ivar_r[indx,w] = 0                                
            # Mask falls both in B and R regions:     
            else:
                #print('Mask both B and R filters')
                w1 = np.where((wave_b > lmin_tmp2) & (wave_b < lmax_tmp2))
                w2 = np.where((wave_r > lmin_tmp2) & (wave_r < lmax_tmp2))
                ivar_b[indx,w1] = 0                
                ivar_r[indx,w2] = 0                

        elif(lmax_tmp3 < lmin_Z).any():
            # Mask falls only in the R region:
            if(lmin_tmp3 > lmax_B):
                #print('Mask only goes in R filter')
                w = np.where((wave_r > lmin_tmp3) & (wave_r < lmax_tmp3))
                ivar_r[indx,w] = 0                                
            # Mask falls both in B and R regions:     
            else:
                #print('Mask both B and R filters')
                w1 = np.where((wave_b > lmin_tmp3) & (wave_b < lmax_tmp3))
                w2 = np.where((wave_r > lmin_tmp3) & (wave_r < lmax_tmp3))
                ivar_b[indx,w1] = 0                
                ivar_r[indx,w2] = 0                
                
        #Mask both filters R and Z:         
        elif(lmax_tmp < lmax_R).any():
            #print('Mask both R and Z filters')
            w1 = np.where((wave_r > lmin_tmp) & (wave_r < lmax_tmp))
            w2 = np.where((wave_z > lmin_tmp) & (wave_z < lmax_tmp))
            ivar_r[indx,w1] = 0
            ivar_z[indx,w2] = 0
            
        #Region around Z
        elif(lmax_tmp < lmax_Z).any():
            # Mask falls only in the Z region:
            if(lmin_tmp > lmax_R):
                #print('Mask only goes in Z filter')
                w = np.where((wave_z > lmin_tmp) & (wave_z < lmax_tmp))
                ivar_z[indx,w] = 0
            # Mask falls both in R and Z regions:     
            else:
                #print('Mask both R and Z filters')
                w1 = np.where((wave_r > lmin_tmp) & (wave_r < lmax_tmp))
                w2 = np.where((wave_z > lmin_tmp) & (wave_z < lmax_tmp))
                ivar_r[indx,w1] = 0
                ivar_z[indx,w2] = 0
                
#print(wave_b.shape,ivar_b.shape)
#print(wave_r.shape,ivar_r.shape)
#print(wave_z.shape,ivar_z.shape)
print(lmin_B,lmax_B,lmin_R,lmax_R,lmin_Z,lmax_Z)
#spectra.writeto("spectra-16-0_mod.fits")
spectra.close()        

