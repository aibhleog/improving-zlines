'''
This is a script that can be implemented into your .bashrc and used
to plot a redshifted spectrum and show relevant lines of interest.
--> NOTE that this outputs a PDF file, but feel free to change to PNG

As Taylor mostly works in the NIR and IR, the bandpasses plotted only
cover that wavelength space, although this code can easily be adapted
to fit your needs.

One possible way to do this would be, based upon the redshift and 
resulting spectral range covered, you could only plot the bandpasses
that would show up there.  That could cut down on the number of
unnecessary filters shown in your legend.

Filters for most every telescope/instrument pair can be found using
the SVO Filter Profile Service:
http://svo2.cab.inta-csic.es/svo/theory/fps3/

Planned updates:  I would like to add a quiescent spectrum eventually.
                  
Credit: 	Taylor Hutchison
        aibhleog@tamu.edu
        Texas A&M University
'''

_author_ = 'Taylor Hutchison'


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
import igm_absorption as igm # 	another script written by Taylor Hutchison
                 #	which adds in IGM absorption for the higher redshifts
import os
import add_lines
import filters as filt

import warnings # it's the plt.tight_layout() yelling, we're ignoring it for now
warnings.filterwarnings('ignore',category=UserWarning)


# if you want to run this script in your terminal via an alias, edit this variable
# with the path to this directory -- then you can run the script as a command from
# any location on your computer and it will work perfectly!

#path = '/path/to/file/' # include the '/' at the end!
path = '' # as long as it's in the same folder

def bandpass_zlines(redshift,which='jwst'):
    '''
    Returns a plotted galaxy model spectrum at the specified redshift.
    Underlaid are bandpasses for spectroscopy & photometry from various 
    telescopes.

    INPUTS:
    > redshift ----- the numerical redshift of interest
    > which (opt) -- "jwst" or "mosfire" to specify which bands to plot

    RETURNS:
    > creates & opens a plot showing the redshifted spectrum & bands
    
    '''

    # Cloudy model provided by Taylor Hutchison
    # Single stellar population from BPASS with age=10Myr & IMF that goes to 100 Msolar, 
    # Z(stellar)=Z(nebular)=0.2 Zsolar, ionization parameter log_10(U)=-2.1, n_H=300 cm^(-3)	
    z,zneb,u = 0.2,0.2,-2.1
    model = 'age7z%szneb%su%s_100.con'%(z,zneb,u)
    data = '' # if you place it somewhere else, you can put the path here
    con = np.loadtxt(data+model,usecols=[0,6])
    wave,vfv = con[:,0],con[:,1]
    #print(wave[0]*1e4,wave[-1]*1e4)
    nu = 2.998e+14 / wave
    sed_0 = vfv / nu

    # plotting the redshifted spectrum
    plt.figure(figsize=(16.5,7.5))
    gs00 = gridspec.GridSpec(2,1,height_ratios=[1,0.6],hspace=0.39)
    
    ax = plt.subplot(gs00[0]) # for the bandpasses, spectrum, lines, etc
    ax2 = plt.subplot(gs00[1]) # for the jwst filters
    
    # setting up the axes
    ax.spines['bottom'].set_visible(True)
    ax.spines[['left','right','top']].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.axes.get_yaxis().set_visible(False)
    
    # setting up the bottom area
    #ax2.axis('off')
    ax2.spines['top'].set_visible(True)
    ax2.spines[['left','right','bottom']].set_visible(False)
    ax2.get_xaxis().tick_top()
    ax2.axes.get_yaxis().set_visible(False)
    


    # --- adding filters --- #
    # ---------------------- #
    # plotting the filters that are in the range
    
    filts, filt_total = filt.is_it_in_range(redshift)
    for telinst in filts.keys(): # going through telescope/instruments
        bands = filts[telinst]
        if len(bands) > 0:
            for band in bands:
                f = filt.get_filter(telinst,band)
                f_info = filt.get_filter_info(telinst,band)
                
                scale = (1e-14/max(f.throughput.values)) # just to scale them all the same
                ax.fill_between(f.wave/1e4,f.throughput*scale,0,alpha=0.3,zorder=0,\
                    label=f_info['name'],facecolor=f_info['color'],edgecolor=f_info['edgecolor'])
                ax.plot(f.wave/1e4,f.throughput*scale,alpha=0.8,color=f_info['color'])
    
    # adding jwst filters!
    inst_count = 0
    filts, filt_total = filt.is_it_in_range(redshift,jwst=True)
    for telinst in filts.keys(): # going through JWST/instruments
        bands = filts[telinst]
        if len(bands) > 0:
            band_count = 0
            for band in bands:
                f = filt.get_filter(telinst,band,jwst=True)
                f_info = filt.get_filter_info(telinst,band,jwst=True)
                
                scale = 1.75*inst_count + band_count*0.2
                ax2.plot(f.wave/1e4,f.throughput+scale,alpha=0.4,color=f_info['color'],lw=9)
                
                # adding the label
                if band_count == 0:
                    txt = ax2.text(0.668*(1+redshift)*0.99,f.throughput[0]+scale+0.2,telinst.upper(),
                            ha='right', fontsize=13, color=f_info['color'])
                    txt.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground='k')])
                    ax2.axhline(f.throughput[0]+scale,ls=':',color='grey',lw=1,zorder=0)

                band_count += 1
            inst_count += 1
    
    

    # applying IGM absorption depending upon z
    sed = sed_0 * igm.igm_absorption(wave*1e4*(1+redshift),redshift)
    ax.plot(wave*(1+redshift),sed,color='k',lw=2.)

    # adding emission line notations
    add_lines.modify_yval('halpha',-0.4)
    add_lines.modify_yval('nii1',-0.3)
    add_lines.rename_line('halpha',r'H$\alpha$')
    # add_lines.rename_line('nii1','[NII]$\lambda$6550')
    add_lines.draw_lines(ax,2e-14,redshift,xlim=['lya','halpha']) # adding in the line labels

    ax.set_yscale('log')
    #ax.set_yticklabels([])
    ax.tick_params(labelsize=16)
    ax.set_xlabel(f'observed wavelength for $z=\,${redshift} [microns]',labelpad=5,fontsize=16)
    ax.set_xlim(0.08*(1+redshift),0.668*(1+redshift))
    ax.set_ylim(1e-15,3.5e-13)
    
    ax2.tick_params(labelsize=16)
    ax2.set_xlim(0.08*(1+redshift),0.668*(1+redshift))
    ax2.set_ylim(7,-0.5)
    
    leg = ax.legend(frameon=False,loc=9,fontsize=11,ncol=filt_total,bbox_to_anchor=(0.5,1.2)) # was 1.12
    for lh in leg.legend_handles: 
        lh.set_alpha(1)
    

    plt.tight_layout()
    # plt.show() # uncomment this if you want to edit the script
    plt.savefig(path+'figure.pdf')
    plt.close('all')

    # opening image from the terminal, comment out if you're in editing mode
    # os.system(f'gnome-open {path}figure.pdf') # if linux
    os.system(f'open {path}figure.pdf') # if mac


# reads in input for scripted version
if __name__ == "__main__":
    import sys
    try:
        if sys.argv[1] == 'help':
            print('''
Given a redshift, this command will show a spectrum that
has been redshifted with relevant lines marked. In addition,
it will show relevant NIR and IR bandpasses.

Use the following notation:   zlines [redshift]
''')	
        else: bandpass_zlines(float(sys.argv[1]))
    except IndexError: # if you don't list a redshift
        z = 7.5032 # my favorite redshift, see Hutchison et al. 2019
        print('Redshift not specified, set to z=7.5032',end='\n\n')
        bandpass_zlines(z)





