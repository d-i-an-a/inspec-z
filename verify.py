from spectrum import Spectrum
from catalog_io import read_1d, save_file
from astropy.table import Table
import os
import yaml

class Verify(Spectrum):
 
    def __init__(self, fcat, config):
        import matplotlib
        import numpy as np

        Spectrum.__init__(self)

        # Plotting parameters        
        matplotlib.rcParams['backend'] = 'TkAgg'
        matplotlib.rcParams['figure.autolayout'] = False
        # spectrum smooth factor
        self.box_size = config['box_size']
        self.sm_temp = self.box_size
        # padding around 1d spec plot
        self.padding = config['padding']

        # Input files
        self.fcat = fcat
        self.config = config
        self.path =self.config['1dpath']

        # Emission lines for plot
        emlines_file = self.config['emlines_file']
        emlines = Table.read(emlines_file, format=self.config['emlines_file_type'])
    
        if self.config['emlines_show'].lower() == 'primary':
            emlines = emlines[emlines['primary']==1]

        if self.config['wavelength_type'].lower() == 'air':
            self.emlines_wave = emlines['wave_air']
        elif self.config['wavelength_type'].lower() == 'vacuum':
            self.emlines_wave = emlines['wave_vac']
        else:
            raise ValueError(f"{self.config['wavelength_type']} not supported, use \'vacuum\' or \'air\'.")
        self.emlines_name = emlines['Ion']

        # output files
        self.outfilename = os.path.join(self.config['out_path'], self.config['out_filename'])
        self.buffer_name = os.path.join(self.config['buffer_path'], self.config['buffer_name'])
        # save buffer config file to check consistency when resuming
        buffer_config_name = '.'.join(self.buffer_name.split('.')[:-1])+'.yaml'
        with open(buffer_config_name, 'w') as yaml_file:
            yaml.dump(config, yaml_file, default_flow_style=False)

        if os.path.exists(self.buffer_name):

            # if buffer exists, recover last entry
            self.buffer = Table.read( self.buffer_name, format=self.config['buffer_format'])
            self.buffer.add_index('idx')

            # find last entry
            self.idx = self.buffer['idx'][-1]
            self._get_buffer()        
        else:
            self.idx = int(0)
            
            # if buffer does not exist, start from the top
            self.buffer = Table(names=('idx', 'inspecz_id', 'id_', 
                                        'z', 'z_temp', 'z_phot',
                                        'flag', 'flag_temp', 
                                        'verified', 'verified_temp'),
                                dtype=('i', np.int64, 'U100',
                                        'f', 'f', 'f',
                                        'i', 'i',
                                        'i', 'i' ))

            self.id_ = self.fcat[self.config['ID']][self.idx]
            self.inspecz_id = self.fcat['inspecz_id'][self.idx]
            self.z = self.fcat[self.config['zspec']][self.idx]
            self.z_temp = self.z
            try:
                self.zphot = self.fcat[self.config['zphot']+'_zphot'][self.idx]
            except:
                self.zphot = -99

            self.flag = int(self.fcat[self.config['new_flag_name']][self.idx])
            self.flag_temp = self.flag
            self.verified = int(self.fcat['verified'][self.idx])
            self.verified_temp = self.verified
            
            self.buffer.add_row((self.idx, self.inspecz_id, str(self.id_), 
                                self.z, self.z_temp, self.zphot,
                                self.flag, self.flag_temp, 
                                self.verified, self.verified_temp))                    

            self.buffer.write(self.buffer_name, format=self.config['buffer_format'], overwrite=True)

            self.buffer.add_index('idx')
            
        # read 1d spectrum
        try:
            self.fsp = self.fcat[self.config['1dname']][self.idx].strip()

            self.wave, flux, self.mask, self.error = read_1d(os.path.join(self.path, self.fsp.strip()), self.config['1dsp_profile'])
            self.flux = flux * float(self.config['flux_calibration'])

            self.set_1d(self.wave, self.flux)
        except:
            raise ValueError(f'File {self.fcat[self.config["1dname"]][self.idx]} not found.')

    
    def _get_buffer(self):

        current_idx = self.buffer.loc_indices[self.idx]
        
        self.id_ = self.buffer['id_'][current_idx]
        self.inspecz_id = self.buffer['inspecz_id'][self.idx]
        self.z = self.buffer['z'][current_idx]
        self.z_temp = self.buffer['z_temp'][current_idx]
        try:
            self.zphot = self.buffer['z_phot'][current_idx]
        except:
            self.zphot = -99
        self.flag = self.buffer['flag'][current_idx]
        self.flag_temp = self.buffer['flag_temp'][current_idx]
        self.verified = int(self.buffer['verified'][current_idx])
        self.verified_temp = int(self.buffer['verified_temp'][current_idx])
        return

    def _update_buffer(self):
        try:
            self.buffer.loc[self.idx] = [self.idx, self.inspecz_id, self.id_, 
                                        self.z, self.z_temp, self.zphot, 
                                        self.flag, self.flag_temp, 
                                        self.verified, self.verified_temp]
            
        except:        
            self.buffer.add_row((self.idx, self.inspecz_id, str(self.id_), 
                            self.z, self.z_temp, self.zphot,
                            self.flag, self.flag_temp, 
                            self.verified, self.verified_temp)) 
        
        self.buffer.write(self.buffer_name, format=self.config['buffer_format'], overwrite=True)

    def _update_text(self):
        # update info box
        textstr = '\n'.join((
            r'$%i$/$%i$' % (self.idx+1, len(self.fcat)),
            r'$zspec=%.4f$' % (self.z, ),
            r'$zspec_{new}=%.4f$' % (self.z_temp, ),
            '',
            r'$zphot=%.4f$' % (self.zphot, ),
            '',
            r'$zflag=%i$' % (self.flag, ),
            r'$zflag_{new}=%i$' % (self.flag_temp, ),
            ))
        self.text_info.set_text(textstr)
        pass

    def _update_flag(self, f):
        self._set_flag(f)
        self._update_text()
        pass

    def _update_1d(self,update_axes=True):
        from astropy.convolution import convolve, Box1DKernel
        from astropy.stats import sigma_clip
        import numpy as np
        from emission_lines import plot_lines
        import matplotlib.pyplot as plt

        for l in self.spec1d_plot:
            l.set_xdata(self.wave)
            l.set_ydata(self.flux)

        # Smoothed spectrum
        kernel = Box1DKernel(self.sm_temp)
        self.smoothed_flux = convolve(self.flux, kernel)
        
        for l in self.spec1d_smoothed:
            l.set_xdata(self.wave)
            l.set_ydata(self.smoothed_flux)
       
        if update_axes == True:
            self.ax.set_xlim([np.min(self.sp1d_wave), np.max(self.sp1d_wave)])

            # set axis limits ignoring extreme values 
            self.flux_sigma_clipped = sigma_clip(self.smoothed_flux, sigma=self.config['sigma_clip'])
            self.ax.set_ylim([np.min(self.flux_sigma_clipped)*(1-self.padding), 
                              np.max(self.flux_sigma_clipped)*(1+self.padding)])
        
            self._update_emlines(self.z_temp)
        
        # update title
        self.ax.set_title(self.id_)

        self._update_text()
        
        # Update legend
        # get old labels
        handles, labels = self.ax.get_legend_handles_labels()
        # update smoothing parameter
        new_labels = []
        for label in labels:
            if 'box_size' in label:
                label = f'box_size {self.sm_temp}'
            new_labels.append(label)
        # Draw new legend
        self.ax.legend(handles,new_labels, loc=3)

        # Update plot
        plt.draw()
        

    def _update_emlines(self, val):
        import numpy as np
        import matplotlib.pyplot as plt

        if self.config['emlines_show'].lower() != 'none':
            # get current limits from figure
            xlims = self.ax.get_xlim()
            ylims = self.ax.get_ylim()

            old_x = self.lines[0].get_xdata()
            z_past = np.abs(self.emlines_wave[0] - old_x[0])/self.emlines_wave[0]

            for line, text in zip(self.lines, self.texts):
                old_x = line.get_xdata()
                (x_text, y_text) = text.get_position()

                new_x  = x_text*(1+val)/(1+z_past)
                new_y = ylims[1]*(1-self.padding)

                line.set_xdata([new_x for x in old_x])
                text.set_position((new_x, new_y))

                if new_x < xlims[0] or new_x > xlims[1]:
                    text.set_visible(False)
                else:
                    text.set_visible(True)

            self.z_temp = val
            self._update_text()
            self.zbox.set_val(self.z_temp)

            plt.draw()
        pass

    def _reset_z(self, event):
        
        self.z_temp = self.z
        self.sredshift.set_val(self.z_temp)
        self._update_emlines(self.z_temp)
        self._update_buffer()
        pass

    def _set_value_z(self, val):
        
        self.z_temp = round(float(val),6)
        self.sredshift.set_val(self.z_temp)
        self.zbox.set_val(self.z_temp)
        self._update_emlines(self.z_temp)
        pass

    def _set_value_box_size(self, val):
        try:
            self.sm_temp = int(val)
            self.sbox.set_val(self.sm_temp)
            self._update_1d(update_axes=False)
        except:
            pass
       
    def _save(self, event):
        # save progress in buffer
        self._update_buffer()
        pass

    def _save2(self, event):
        # set verified flag, save to buffer, and advance to next spectrum
        self._set_verified('Y')
        self._update_buffer()
        
        if self.idx+1 < len(self.fcat):
            self.idx += 1
            self.get_1d()
        else:
            print('Reached the end of the catalog.')

    def _save_buffer(self, event):
        self._update_buffer()

        # save buffer to final file, overwrites file
        key_dict = { 'updated_zspec': 'z_temp',
                 'updated_flag': 'flag_temp',
                 'verified': 'verified_temp'}
        save_file(self.fcat, self.buffer, key_dict, 'verified_temp', self.outfilename, self.config['out_format'])
        pass
    #######################################################
    
    def get_1d(self):
        
        try:
            self._get_buffer()
        except:
            self.id_ = self.fcat[self.config['ID']][self.idx]
            self.inspecz_id = self.fcat['inspecz_id'][self.idx]
            self.z = self.fcat[self.config['zspec']][self.idx]
            self.z_temp = self.z
            try:
                self.zphot = self.fcat[self.config['zphot']+'_zphot'][self.idx]
            except:
                self.zphot = -99
            self.flag = int(self.fcat[self.config['new_flag_name']][self.idx])
            self.flag_temp = self.flag
            self.verified = int(self.fcat['verified'][self.idx])
            self.verified_temp = self.verified
            self._update_buffer()

        try:
            # read spectrum
            self.fsp = self.fcat[self.config['1dname']][self.idx]

            self.wave, flux, self.error, self.mask = read_1d(os.path.join(self.path, self.fsp), self.config['1dsp_profile'])
            self.flux = flux * float(self.config['flux_calibration'])
            self.set_1d(self.wave, self.flux)
        except:
            print(f'File {self.fcat[self.config["1dname"][self.idx]]} not found.')


        # update window
        self.sredshift.set_val(self.z_temp)
        self.radio.set_active(self.flag_temp)
        self.radio2.set_active(self.verified_temp)

        self._update_1d()
        

    def next(self, event):
        if self.idx+1 < len(self.fcat):
            self.idx += 1
            self.get_1d()
        else:
            print('Reached the end of the catalog.')
            pass        


    def prev(self, event):
        if self.idx-1 >= 0:
            self.idx -= 1
            self.get_1d()        
        else:
            print('No more previous sources.')

            pass

    def _on_ylim_change(self, *args):  
        self._update_emlines(self.z_temp)
    #######################################################


    def verify(self, show=False):
        """
        plot interactive window with:
        observed frame spectrum
        emission lines at redshift of source
        button for verification 0/1
        """
        
        import numpy as np
        import matplotlib.pyplot as plt
        from emission_lines import plot_lines
        from matplotlib.widgets import RadioButtons, CheckButtons, Slider, Button
        from matplotlib.widgets import TextBox
        import matplotlib.gridspec as gridspec
        from astropy.convolution import convolve, Box1DKernel
        from astropy.stats import sigma_clip

        fig = plt.figure(figsize=(12,6), constrained_layout=False)

        fig.subplots_adjust(left=0.065, bottom=0.2, top=0.95, right=0.98)
        gs = fig.add_gridspec(ncols=6, nrows=1, figure=fig)
        self.ax = fig.add_subplot(gs[0, :5])

        # Raw spectrum
        # TODO: add mask
        # TODO: add error
        self.spec1d_plot = self.ax.plot(self.sp1d_wave, self.sp1d_flux, '', 
                                        color='gray', alpha=0.5, label='raw')

        # Smoothed spectrum
        kernel = Box1DKernel(self.box_size)
        self.smoothed_flux = convolve(self.sp1d_flux, kernel)
        self.spec1d_smoothed = self.ax.plot(self.sp1d_wave, self.smoothed_flux, '', 
                                        color='k', alpha=1, label=f'box_size {self.sm_temp}')
        self.flux_5s_clipped = sigma_clip(self.smoothed_flux, sigma=self.config['sigma_clip'])

        # Mask
        # FUTURE: highlight masked regions

        # Emission lines
        self.lines, self.texts = plot_lines(ax=self.ax, 
                    redshift=self.z_temp, 
                    wave = self.emlines_wave, 
                    name = self.emlines_name, 
                    min_wave = np.min(self.sp1d_wave), 
                    max_wave = np.max(self.sp1d_wave), 
                    bottom = np.min(self.flux_5s_clipped), 
                    top = np.max(self.flux_5s_clipped)*(1-self.padding), 
                    offset = None, 
                    style = 'emission'
                    )

        plt.legend(loc=3)
        plt.title(self.id_)

        # redshift slider
        ax_z = plt.axes([0.15, 0.075, 0.6, 0.03], facecolor=None)
        self.sredshift = Slider(ax_z, 'redshift', 
                                self.config['slider_zmin'], self.config['slider_zmax'], 
                                valinit = self.z,
                                valstep=5e-6, dragging=True,
                                valfmt='%1.4f',
                                capstyle='round')
        self.sredshift.vline.set_ls('')              
        self.sredshift.vline.set_marker(None)
        self.sredshift.on_changed(self._update_emlines)


        #
        # Info box
        #
        self.textbox_ax = plt.axes([0.85, 0.95, 0.1, 0.2])
        self.textbox_ax.set_axis_off()
        textstr = '\n'.join((
                    r'$%i$/$%i$' % (self.idx+1, len(self.fcat)),
                    r'$zspec=%.4f$' % (self.z, ),
                    r'$zspec_{new}=%.4f$' % (self.z_temp, ),
                    '',
                    r'$zphot=%.4f$' % (self.zphot, ),
                    '',
                    r'$zflag=%i$' % (self.flag, ),
                    r'$zflag_{new}=%i$' % (self.flag_temp, ),
                    ))
        props = dict(boxstyle='round', facecolor='gray', alpha=0.15)
        self.text_info = self.textbox_ax.text(0,0, textstr, 
                            verticalalignment='top', bbox=props)
        
        # Interactivity
        #
        # Radio button for quality flag
        # [left, bottom, width, height]
        rax = plt.axes([0.85, 0.45, 0.1, 0.2], facecolor=None)
        rax.text(0,1.05,"Quality Flag")
        # FUTURE: dymanic allocation from configuration file
        self.radio = RadioButtons(rax, ('0', '1', '2', '3', '4', '5'), active=self.flag_temp)
        self.radio.on_clicked(self._update_flag)

        # Radio button for verification
        rax2 = plt.axes([0.85, 0.3, 0.1, 0.1], facecolor=None)
        rax2.text(0,1.05,"Verified")
        self.radio2 = RadioButtons(rax2, ('N', 'Y'), active=self.verified_temp)
        self.radio2.on_clicked(self._set_verified)


        # reset zspec button
        axreset = plt.axes([0.85, 0.25, 0.1, 0.05])
        breset = Button(axreset, 'Reset z')
        breset.on_clicked(self._reset_z)

        # save to buffer button
        axsave = plt.axes([0.85, 0.2, 0.1, 0.05])
        bsave = Button(axsave, 'Save to buffer')
        bsave.on_clicked(self._save)
        
        # one step verify, save to buffer, next spectrum button
        axsave = plt.axes([0.64, 0.0125, 0.2, 0.05])
        bsave2 = Button(axsave, 'Save to buffer, verified, next')
        bsave2.on_clicked(self._save2)

        # interactive redshift value box
        ax_zbox = plt.axes([0.82, 0.13, 0.075, 0.05], facecolor=None)
        self.zbox = TextBox(ax_zbox, "z:")
        self.zbox.set_val(round(self.z_temp, 6))
        self.zbox.on_submit(self._set_value_z)

        # box smooth value control
        ax_smbox = plt.axes([0.945, 0.13, 0.03, 0.05], facecolor=None)
        self.sbox = TextBox(ax_smbox, "smooth:")
        self.sbox.set_val(self.sm_temp)
        self.sbox.on_submit(self._set_value_box_size)
        
        # previous spectrum 
        axprev = plt.axes([0.82, 0.07, 0.075, 0.05])
        bprev = Button(axprev, 'Previous')
        bprev.on_clicked(self.prev)

        # next spectrum
        axnext = plt.axes([0.9, 0.07, 0.075, 0.05])
        bnext = Button(axnext, 'Next')
        bnext.on_clicked(self.next)

        # save buffer to catalog
        axsave = plt.axes([0.85, 0.0125, 0.1, 0.05])
        bsave_file = Button(axsave, 'Save progress')
        bsave_file.on_clicked(self._save_buffer)

        # Update emission lines when zooming in
        self.ax.callbacks.connect('ylim_changed', self._on_ylim_change)
        
        self.ax.set_xlabel('wavelength')
        self.ax.set_ylabel('flux')

        self.ax.set_xlim([np.min(self.sp1d_wave), np.max(self.sp1d_wave)])
        self.ax.set_ylim([np.min(self.flux_5s_clipped)*(1-self.padding), 
                          np.max(self.flux_5s_clipped)*(1+self.padding)])
        
        #plt.draw()
        #if show: plt.show()
        plt.show()
        return 


if __name__ == "__main__":

    """
    Test catalog
    """    
    import os
    import yaml

    import matplotlib.pyplot as plt
    from catalog_io import get_catalog


    # configuration
    yaml_file = open("config_files/example.yaml", 'r')
    config = yaml.load(yaml_file, Loader=yaml.FullLoader)

    fcat = get_catalog(config)

    s1 = Verify(fcat, config)
    s1.verify(show=True)
