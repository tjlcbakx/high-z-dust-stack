import numpy as np
from astropy.io import fits
# import LineStacker
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord 
from astropy.wcs import WCS
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
import matplotlib.pyplot as plt


fontSize = 14
fontBlack = {'family' : 'sans-serif',
#        'serif' = 'Helvetica Neue'
'weight' : 'normal',
'size'   : fontSize,
}

#class for easier book keeping
class Galaxy:
    def __init__(self, ra=-1, dec=-1, z=-1, image='', weight=-1, ra_pix=-1, dec_pix=-1):
        self.ra = ra #RA of source
        self.dec = dec #Dec of source
        self.z = z  #redshift
        self.image = image #name of the fits file
        self.weight = weight #weight for the stack
        self.ra_pix = ra_pix #RA pixel coordinate of the source
        self.dec_pix = dec_pix #Dec pixel coordinate of the source
        self.name=image.replace('.fits','') #name of the source

    #function to retrieve the coordinates in pixel
    def get_pix_coords(self): 
        open_data = fits.open(self.image)
        head = open_data[0].header
        w = WCS(head)
        open_data.close()
        sky_coord = SkyCoord(self.ra, self.dec, frame='icrs')
        self.ra_pix, self.dec_pix = sky_coord.to_pixel(w)
        self.ra_pix, self.dec_pix = int(self.ra_pix), int(self.dec_pix)

    #function that retrieves the stamp to stack
    def get_to_stack(self, stamp_size=91):
        open_data = fits.open(self.image)
        data = open_data[0].data
        #data is transposed to have ra, dec, stokes, frequency
        if len(data.shape)==4:
            data = np.transpose(data, (3, 2, 0, 1))
        else:
            data=np.reshape(data, (data.shape[0], data.shape[1],1,1))
        self.data=data
        open_data.close()
        half_stamp = int(stamp_size/2)
        assert data.shape[2:]==(1,1), f'wrong shape, {data.shape}, {data.shape[2:]}'
        self.to_stack = data[self.ra_pix-half_stamp:self.ra_pix+half_stamp+1, self.dec_pix-half_stamp:self.dec_pix+half_stamp+1, 0 ,0]
        return self.to_stack



stackFiles = np.array([
'all.txt',
'158um.txt',
'88um.txt',
'lowBeta.txt',
'highBeta.txt',
'lowMass.txt',
'highMass.txt',
'lowRedshift.txt',
'highRedshift.txt'])



stackFilesNames = np.array([
'All',
r'158$\mu$m',
r'88$\mu$m',
r'$\beta_{\rm UV} < -2$',
r'$\beta_{\rm UV} > -2$',
r'$M_{\rm \ast} < 10^{9}$ M$_{\odot}$',
r'$M_{\rm \ast} > 10^{9}$ M$_{\odot}$',
r'$z < 9.5$',
r'$z > 9.5$'])

import matplotlib.patheffects as patheffects
import matplotlib.cm as cm

for q in ['reprojectCont','dustMassFits','dustToStellarMassFits']:
    for abc in range(len(stackFiles)):
        #the file containing the coordinates etc. of each target
        coord_file=stackFiles[abc]
        coord_file_info=pd.read_csv(coord_file, sep=',', header=0)
        #will contain all Galaxy class objects
        all_galaxies=[]
        #the size of the stamp used for stacking, in pixel
        stamp_size=121
        #the numpy array as will be filled and before averaging, size of number_of_galaxies*stamp_size*stamp_size
        to_stack=np.zeros((len(coord_file_info), stamp_size, stamp_size))
        #do_sigma2_w is set to True if we want to set the weights to 1/sigma^2
        do_sigma2_w=True
        method='mean' #stacking method, can be 'mean' or 'median' (note that median don't use weights)
        #parsing the coordinate file one row at a time
        for row in range(len(coord_file_info)):
            this_row = coord_file_info.iloc[row]
            #extracting RA and Dec for each source
            ra, dec = this_row['ra'], this_row['dec'] 
            #making it into a astropy SkyCoord
            coord = SkyCoord( ra=ra,  dec=dec,  unit=(u.hourangle, u.deg), frame='icrs') 
            #making a Galaxy class object
            this_galaxy=Galaxy(ra=coord.ra, dec=coord.dec, z=this_row['redshift'], image='../'+q+'/'+this_row['#fitsname'], weight=this_row['weighting']) 
            #getting the pixel coordinate for the target
            this_galaxy.get_pix_coords()
            #getting the stack stamp for the target and filling the numpy array that will be averaged later on
            to_stack[row,:,:]=this_galaxy.get_to_stack(stamp_size=stamp_size)
            # if q == 'dustToStellarMassFits':
            #     if this_row['#fitsname'] == 'J2140+0241_cont.fits':
            #         to_stack[row,:,:]= 1e200 * this_galaxy.get_to_stack(stamp_size=stamp_size)
            #if do_sigma2_w, compute the sigma2 and using as weight
            if do_sigma2_w:
                #the lower 5 flux of the image, to avoid artefacts
                low5=np.nanpercentile(this_galaxy.data,0) 
                #the higher 95 flux of the image, to avoid artefacts/bright sources
                high5=np.nanpercentile(this_galaxy.data,100)
                #to std conists of all pixels in the image with 5%>flux<95%
                to_std=this_galaxy.data[this_galaxy.data>low5]
                to_std=to_std[to_std<high5]
                #weights are 1/sigma2
                sigma_weight=1/np.std(to_std)**2
                #setting the weight to the object
                this_galaxy.weight=sigma_weight
            #adding the galaxy object to the list of all galaxies
            all_galaxies.append(this_galaxy)
        #depending on stacking method chosen do the stack
        if method=='mean':
            #using weights
            stack=np.average(to_stack, axis=0, weights=[gal.weight for gal in all_galaxies])
        elif method=='median':
            stack=np.average(to_stack, axis=0)
        #plotting the stack
        if q == 'reprojectCont':
            texttoshow = r'$S_{\nu}$ stack'+ '\n' 
        elif q == 'dustMassFits':
            texttoshow = r'$M_{\rm dust}$ stack'+ '\n' 
        elif q == 'dustToStellarMassFits':
            texttoshow = r'$M_{\rm dust} / M_{\ast}$ stack'+ '\n' 
        texttoshow  = texttoshow + stackFilesNames[abc]
        fig=plt.figure(figsize=(4.5,4.5))
        ax=fig.add_subplot(111)
        val1 = -1*np.percentile(stack.T[np.isnan(stack.T) == False],15.865)#
        val2 = np.percentile(stack.T[np.isnan(stack.T) == False],84.135)#
        val3 = np.nanstd(stack.T)
        # val3 = np.nanstd(data[0].data[int(datalen/2 - 25):int(datalen/2 + 25),int(datalen/2 - 25):int(datalen/2 + 25)])
        stdtemp = min(val1,val2,val3)
        variableName = ax.imshow(stack.T, origin='lower', cmap='RdGy_r',vmin=-3*stdtemp,vmax=3*stdtemp)
        ax.contour(stack.T, levels=stdtemp*np.array([-4,-3,-2,2,3,4]), origin='lower', colors='w',linewidths=1)
        ax.set_xticks([20, 30,40,50,60,70,80,90,100])
        ax.set_xticklabels([r'$-$4',r'$-$3',r'$-$2',r'$-$1','0','1','2','3','4'],fontsize=14)
        ax.set_xlim(20,100)
        ax.set_xlabel('RA [arcsec]',fontsize=16)
        ax.set_yticks([20, 30,40,50,60,70,80,90,100])
        ax.set_yticklabels([r'$-$4',r'$-$3',r'$-$2',r'$-$1','0','1','2','3','4'],fontsize=14)
        ax.set_ylim(20,100)
        ax.set_ylabel('Dec [arcsec]',fontsize=16)
        ax.plot([60,60],[50,40],lw=3,color='k',zorder=30)
        ax.plot([70,80],[60,60],lw=3,color='k',zorder=30)
        ax.text(25,90 ,texttoshow,
        horizontalalignment='left',
        verticalalignment='center',
        rotation=0,
        size=24,
        color='w',
        path_effects=[patheffects.withStroke(linewidth=5, foreground='k', capstyle="round")],
        fontdict = fontBlack)
        # ax.set_colorbar(label='Flux density [µJy]')
        # fig.colorbar(variableName,label='Flux density [µJy]')
        fig.tight_layout()
        val1 = -1*np.percentile(stack.T[np.isnan(stack.T) == False],15.865)#
        val2 = np.percentile(stack.T[np.isnan(stack.T) == False],84.135)#
        val3 = np.nanstd(stack.T)
        # val3 = np.nanstd(data[0].data[int(datalen/2 - 25):int(datalen/2 + 25),int(datalen/2 - 25):int(datalen/2 + 25)])
        stdtemp = min(val1,val2,val3)
        print('stack'+q+'_'+coord_file[:-4] + ' ' + str(stdtemp*3) + ' ' + str(stack.T[60,60]) + ' ' + str(stack.T[60,60]/stdtemp))
        fig.savefig(f'stack_'+q+'_'+coord_file[:-4]+'.fits'+'.png')
        fig.savefig(f'stack_'+q+'_'+coord_file[:-4]+'.fits'+'.pdf',transparent=True)
        # plt.show()
        plt.close(fig)
        #writting the stack to a fits file
        hdu = fits.PrimaryHDU(data=stack)
        hdu.writeto('stack'+q+'_'+coord_file[:-4]+'.fits', overwrite=True)
        # hdu.close()


# plt.show()
# 







#if set plot_all to True just plot all the individual stamp for each galaxy, useful to identify bright sources in the field etc.
plot_all=False
if plot_all:
    for gal in all_galaxies:
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.imshow(gal.to_stack.T, origin='lower', cmap='RdGy_r')
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(f'{gal.name}.png')
        plt.close(fig)