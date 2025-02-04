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


# if you want to run this script in your terminal via an alias, edit this variable
# with the path to this directory -- then you can run the script as a command from
# any location on your computer and it will work perfectly!

#path = '/path/to/file/' # include the '/' at the end!
path = '' # as long as it's in the same folder

def bandpass_zlines(redshift,which='spec'):
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
    if which == 'phot': phot = True
    else: phot = False

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
    plt.figure(figsize=(16.5,7))
    ax = plt.gca()
    ax.set_zorder(10)
    ax.set_frame_on(False)
    
    # setting up the axes
    ax.spines['bottom'].set_visible(True)
    ax.spines[['left','right','top']].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.axes.get_yaxis().set_visible(False)

    xlims = [0.08*(1+redshift),0.668*(1+redshift)]
    ax.set_xlim(xlims)
    ax.set_ylim(7.4,-0.6)


    # adding model spectrum
    # -----------------------
    ax2 = ax.twinx()
    ax2.set_zorder(0)
    ax2.spines['bottom'].set_visible(True)
    ax2.spines[['left','right','top']].set_visible(False)
    ax2.get_xaxis().tick_bottom()
    ax2.axes.get_yaxis().set_visible(False)
    # applying IGM absorption depending upon z
    sed = sed_0 * igm.igm_absorption(wave*1e4*(1+redshift),redshift)
    ax2.plot(wave*(1+redshift),sed,color='k',lw=2.,alpha=0.6)
    ax2.set_yscale('log')
    ax2.set_ylim(6e-16,2.5e-13)


    # --- adding filters --- #
    # ---------------------- #
    # plotting the JWST filters that are in the range
    inst_count = 0
    filts, filt_total = filt.is_it_in_range(redshift,jwst=True,phot=phot)
    for telinst in filts.keys(): # going through JWST/instruments
        bands = filts[telinst]
        if len(bands) > 0:
            band_count = 0
            for band in bands:
                f = filt.get_filter(telinst,band,jwst=True,phot=phot)
                f_info = filt.get_filter_info(telinst,band,jwst=True,phot=phot)

                # cutting the JWST throughputs until they're 1/2 or more
                f = f.query(f'throughput > {np.nanmax(f.throughput)/2}').copy()
                f.reset_index(inplace=True,drop=True)
                f['throughput'] = np.ones(len(f)) # for when we're plotting just horz. lines

                flip = 1
                if band_count % 2 == 1: flip = -1
                scale = 1.75*inst_count + flip*0.3 #band_count*0.2
                ax.plot(f.wave/1e4,f.throughput+scale,alpha=0.4,color=f_info['color'],lw=9)
                
                # adding the label
                if band_count == 0:
                    txt = ax.text(0.668*(1+redshift)*0.99,
                            f.throughput[0]+scale-0.25,telinst.upper(),zorder=10,
                            ha='right', fontsize=13, color='#a6d3a9',weight='demibold')
                    txt.set_path_effects([PathEffects.withStroke(linewidth=3.5, foreground='k')])
                    ax.axhline(f.throughput[0]+scale,ls=':',color='grey',
                               alpha=0.6,lw=1,zorder=0)

                # adding filter name
                # have to do some noodling to adjust for edge cases & redshifts
                if band == 'lrs' or band == 'mrs': flip *= -6
                y_scale = 1.1 + 1.75*inst_count + flip*0.1*-1
                x_scale = np.nanmedian(f.wave/1e4)-(0.015*(redshift+1))
                if x_scale < xlims[0]: x_scale = xlims[0]
                elif x_scale > xlims[1]: x_scale = xlims[1]-0.2
                
                txt = ax.text(x_scale,y_scale,band.upper(),
                        fontsize=10, color=f_info['color'],alpha=0.2)
                txt.set_path_effects([PathEffects.withStroke(linewidth=1.5, 
                                        foreground=f_info['color'],alpha=0.2)])

                band_count += 1
            inst_count += 1
    

    
    add_lines.modify_yval('lya',0)
    add_lines.modify_yval('hbeta',0.2)
    add_lines.modify_yval('halpha',0.5)
    add_lines.draw_lines(ax,4,redshift,xlim=['lya','halpha']) # adding in the line labels
    
    ax.tick_params(labelsize=16)
    ax.set_xlabel(f'observed wavelength for $z=\,${redshift} [microns]',labelpad=5,fontsize=16)
    

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

Use the following notation:   zlines [redshift] [(opt) spec/phot]
''')	
        try:
            len(sys.argv[2])
            bandpass_zlines(float(sys.argv[1]),str(sys.argv[2]))
        except:
            bandpass_zlines(float(sys.argv[1]))
    except IndexError: # if you don't list a redshift
        z = 7.5032 # my favorite redshift, see Hutchison et al. 2019
        print('Redshift not specified, set to z=7.5032',end='\n\n')
        bandpass_zlines(z)





