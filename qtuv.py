#!/usr/bin/env python
"""
qtuv.py
========

Qt based GUI for viewing interferometric visibility data.

(c) 2013 Danny Price
dprice@cfa.harvard.edu

"""



# Imports
import sys, os
from optparse import OptionParser
from datetime import datetime

# Python metadata

__version__  = 'v1.0'
__author__   = 'Danny Price'
__email__    = 'dprice@cfa.harvard.edu'
__license__  = 'GNU GPL'
__modified__ = datetime.fromtimestamp(os.path.getmtime(os.path.abspath( __file__ )))
progname     = 'InterFits viewer'

try:
    import lib.qt_compat as qt_compat
    QtGui = qt_compat.import_module("QtGui")
    QtCore = qt_compat.QtCore
    
    USES_PYSIDE = qt_compat.is_pyside()
    
    #import PyQt4
    #from PyQt4 import QtGui, QtCore
except:
    print "Error: cannot load PySide or PyQt4. Please check your install."
    exit()
    
try:    
    import numpy as np
except:
    print "Error: cannot load Numpy. Please check your install."
    exit()

try:    
    import pyfits as pf
except:
    print "Error: cannot load PyFITS. Please check your install."
    exit()

import matplotlib

if matplotlib.__version__ == '0.99.3':
    print "Error: your matplotlib version is too old to run this. Please upgrade."
    exit()
else:
    matplotlib.use('Qt4Agg')
    if USES_PYSIDE:
        matplotlib.rcParams['backend.qt4']='PySide'
    else:
        matplotlib.rcParams['backend.qt4']='PyQt4'
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure
try:
    import matplotlib.pylab as plt
except:
    print "Error: cannot load Pylab. Check your matplotlib install."
    exit()

from interfits import *

class InterFitsGui(QtGui.QWidget):
    """ Qt GUI class
    
    A Qt4 Widget that uses matplotlib to display various information
    
    Parameters
    ----------
    filename: str
        Name of file to open. Defaults to blank, in which case no file is opened.
    """
    def __init__(self, filename=''):
        super(InterFitsGui, self).__init__()
        
        # Initialize user interface
        self.filename = filename
        self.current_row = 0
        self.initUI(width=1200, height=800)
        
       

    def initUI(self, width=1200, height=800):
        """ Initialize the User Interface 
        
        Parameters
        ----------
        width: int
            width of the UI, in pixels. Defaults to 1200px
        height: int
            height of the UI, in pixels. Defaults to 808px
        """
        
        self.setWindowTitle('InterFits %s'%__version__) 
        
        self.main_frame = QtGui.QWidget()   
        self.setWindowIcon(QtGui.QIcon('lib/icon.gif'))  
        #self.gen_gui = generateGui()
        
        # Create buttons/widgets
        self.but_open     = QtGui.QPushButton("Open")
        self.but_open.clicked.connect(self.onButOpen)
        self.but_plot     = QtGui.QPushButton("Plot")
        self.but_plot.clicked.connect(self.updatePlot)
        
        self.lab_info     = QtGui.QLabel(" ")
        
        self.axes_select, self.laxes_select = self.createSpinner("Axis", self.updateAxes, 0, 1, 1)
        self.spin_ref_ant, self.lspin_ref_ant = self.createSpinner("Ant 1", self.updateAxes, 1, 2, 1)
        self.spin_ref_ant2, self.lspin_ref_ant2 = self.createSpinner("Ant 2", self.updateAxes, 1, 2, 1)
        
        self.plot_select = QtGui.QComboBox(self)
        self.plot_select.addItem("Single Baseline")
        self.plot_select.addItem("Single Baseline: dual pol")
        self.plot_select.addItem("Multi baseline: Bandpass")
        self.plot_select.addItem("Multi baseline: Amplitude")
        self.plot_select.addItem("Multi baseline: Phase")
        self.plot_select.addItem("UV coverage")
        self.plot_select.activated.connect(self.updateSpinners)
        
        self.current_plot = ""
        
        self.axes_select = QtGui.QComboBox(self)
        for v in ['Stokes I','Stokes Q','Stokes U','Stokes V']:
            self.axes_select.addItem(v)
        
        # Create plots
        self.sp_fig, self.sp_ax = self.createBlankPlot()
        self.sp_canvas   = FigureCanvas(self.sp_fig)
        self.mpl_toolbar = NavigationToolbar(self.sp_canvas, self.main_frame)
        
        # Widget layout
        layout = QtGui.QVBoxLayout()
        h_layout = QtGui.QHBoxLayout()
        h_layout.addWidget(self.plot_select)
        h_layout.addStretch(1)
        h_layout.addWidget(self.laxes_select)
        h_layout.addWidget(self.axes_select)
        h_layout.addWidget(self.lspin_ref_ant)
        h_layout.addWidget(self.spin_ref_ant)
        h_layout.addWidget(self.lspin_ref_ant2)
        h_layout.addWidget(self.spin_ref_ant2)
        h_layout.addWidget(self.but_plot)
        layout.addLayout(h_layout)
        h_layout = QtGui.QHBoxLayout()
        h_layout.addWidget(self.sp_canvas)
    
        layout.addLayout(h_layout)
        h_layout = QtGui.QHBoxLayout()
        h_layout.addStretch(1)

        layout.addLayout(h_layout)
        layout.addWidget(self.mpl_toolbar)
        
        bbox = QtGui.QHBoxLayout()
        bbox.addWidget(self.lab_info)
        bbox.addStretch(1)
        bbox.addWidget(self.but_open)
        layout.addLayout(bbox)

        self.setLayout(layout)    
        #textEdit = QtGui.QTextEdit()
        #self.setCentralWidget(textEdit)
        #self.setCentralWidget(sp_canvas)
        
        # Load file if command line argument is passed
        if self.filename != '':
            try:
                self.uv = InterFits(self.filename)
                #self.openSdFits(self.filename)
                self.onFileOpen()
                self.plot_single_baseline(1,1)
                self.updateSpinners()
            except:
                print "Error: cannot open %s"%self.filename
                #raise
        
        self.setGeometry(300, 300, width, height)   
        self.show()
        
        def on_click(event):
            """Enlarge or restore the selected axis."""
            ax = event.inaxes
            
            if ax is None:
                # Occurs when a region not in an axis is clicked...
                return
                
            if self.current_plot == 'single':
                if event.button is 1:
                    if not self.ax_zoomed:
                        # Change over to a single baseline plot
                        try:
                            self.ax_zoomed = True
                            ax.set_position([0.1, 0.05, 0.85, 0.80])
                            ax.set_xlabel("Frequency")
                            ax.set_ylabel("Time")
                            
                            for axis in self.sp_fig.axes:
                                if axis is not ax:
                                    axis.set_visible(False)
                                
                        except ValueError:
                            raise
                            self.sp_fig.canvas.mpl_disconnect(self.fig_connect)
                    
                elif event.button is 3:
                    if self.ax_zoomed:
                        self.ax_zoomed = False
                        #self.sp_fig.canvas.mpl_disconnect(self.fig_connect)
                        self.updatePlot()
                        
                else:
                    # No need to re-draw the canvas if it's not a left or right click
                    return
                
            elif self.current_plot == 'multi':
                if ax is None:
                    # Occurs when a region not in an axis is clicked...
                    return
                if event.button is 1:
                    if not self.ax_zoomed:
                        # Change over to a single baseline plot
                        try:
                            ant1, ant2 = ax.get_title().split(" ")
                        except:
                            ant1 = int(ax.get_title())  
                            ant2 = ant1                          
                        try:
                            self.spin_ref_ant.setValue(int(ant1))
                            self.spin_ref_ant2.setValue(int(ant2))
                            self.plot_select.setCurrentIndex(0)
                            self.current_plot = 'single'
                            
                            self.updatePlot()
                        except:
                            raise
                            self.sp_fig.canvas.mpl_disconnect(self.fig_connect)
                
                elif event.button is 3:
                    if not self.ax_zoomed:
                        ax.set_position([0.1, 0.1, 0.85, 0.85])
                        # TODO: fix labelling of zoom plots
                        ax.set_xlabel("Frequency")
                        ax.set_ylabel("Time")
                        self.orig_position = ax.get_position()
                        for axis in event.canvas.figure.axes:
                           # Hide all the other axes...
                           if axis is not ax:
                               axis.set_visible(False)
                        self.ax_zoomed=True
                    else:
                        self.updatePlot()
                    
                else:
                    # No need to re-draw the canvas if it's not a left or right click
                    return
                    
            event.canvas.draw()
                
        
        self.fig_connect = self.sp_fig.canvas.mpl_connect('button_press_event', on_click)
    
    def createSpinner(self, label, action, smin=0, smax=1, step=1):
        """ Create a QtGui.QSpinner with sensible values """
        spinner    = QtGui.QSpinBox()
        spinner.setRange(smin, smax)
        spinner.setSingleStep(step)
        spinner.valueChanged.connect(self.updateAxes)
        label = QtGui.QLabel(label)
        return spinner, label
        
    def createBlankPlot(self):
          """ Creates a single pylab plot for displaying a spectrum """

          fig = plt.figure(figsize=(8,6),dpi=80)
          fig.set_facecolor('#ededed')
          
          # Format plot
          ax = plt.subplot(111)
        
          fig.canvas.draw()
      
          return fig, ax
    
    def plot_single_baseline(self, ant1, ant2, axis=0):
        """ Plot single baseline 
        
        ant1: int
            antenna number of first antenna
        ant2: int
            antenna number of second antenna
        """
        self.current_plot = 'single'
        
        bls = self.uv.d_uv_data['BASELINE']
        
        if ant1 > ant2: bl_id = 256*ant2 + ant1
        else: bl_id = 256*ant1 + ant2
        
        x_data    = self.uv.d_uv_data['DATA'][bls == bl_id,0,0,0,:,axis]  # Baselines, freq and stokes
        x         = x_data[:,:,0] + 1j * x_data[:,:,1]
        
        fig = self.sp_fig
        self.ax_zoomed = False
        
        figtitle = '%s %s: %s -- %s\n'%(self.uv.telescope, self.uv.instrument, self.uv.source, self.uv.date_obs)
        figtitle += "Baseline %i %i"%(ant1, ant2)
        fig.suptitle(figtitle, fontsize=18)
        
        ax = plt.subplot(211)
        x_pow     = np.abs(x)
        #print x_pow.shape
        #x_avg    = np.average(x_pow, axis=0)
        x_max     = np.max(x_pow, axis=0)  
        x_med     = np.median(x_pow, axis=0)
        x_min     = np.min(x_pow, axis=0)  
        
        ax.plot(x_med, label='median')
        ax.plot(x_min, label='min')
        ax.plot(x_max, label='max')
        plt.minorticks_on()
        plt.xlabel("Frequency")
        plt.ylabel("Amplitude")
        plt.legend()
        
        
        ax = plt.subplot(223)
        plt.imshow(x_pow)
        plt.title("Amplitude")
        plt.xlabel("Frequency channel")
        plt.ylabel("Time")
        ax.set_aspect(x.shape[1] / x.shape[0] * 3. / 4)
        plt.colorbar(orientation='horizontal')
        
        ax = plt.subplot(224)
        plt.imshow(np.angle(x))
        plt.title("Phase")
        plt.xlabel("Frequency channel")
        plt.ylabel("Time")
        ax.set_aspect(x.shape[1] / x.shape[0] * 3. / 4)
        cbar = plt.colorbar(orientation='horizontal')
        cbar.set_ticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi-0.05])
        cbar.set_ticklabels(["$-\pi$","$-\pi/2$",0,"$\pi/2$","$\pi$"])
        
        plt.subplots_adjust(left=0.05, right=0.98, top=0.9, bottom=0.05, wspace=0.25, hspace=0.3)
        
        return fig, ax

    def plot_single_baseline_dual_pol(self, ant1, ant2, axis=0):
        """ Plot single baseline 
        
        ant1: int
            antenna number of first antenna
        ant2: int
            antenna number of second antenna
        """
        self.current_plot = 'single'
        
        bls = self.uv.d_uv_data['BASELINE']
        
        if ant1 > ant2: bl_id = 256*ant2 + ant1
        else: bl_id = 256*ant1 + ant2
        
        x_data    = self.uv.d_uv_data['DATA'][bls == bl_id,0,0,0,:,0]  # Baselines, freq and stokes
        x         = x_data[:,:,0] + 1j * x_data[:,:,1]

        y_data    = self.uv.d_uv_data['DATA'][bls == bl_id,0,0,0,:,1]  # Baselines, freq and stokes
        y         = y_data[:,:,0] + 1j * y_data[:,:,1]
        
        fig = self.sp_fig
        self.ax_zoomed = False
        
        figtitle = '%s %s: %s -- %s\n'%(self.uv.telescope, self.uv.instrument, self.uv.source, self.uv.date_obs)
        figtitle += "Baseline %i %i"%(ant1, ant2)
        fig.suptitle(figtitle, fontsize=18)
        
        ax = plt.subplot(221)
        x_pow     = np.abs(x)
        x_max     = np.max(x_pow, axis=0)  
        x_med     = np.median(x_pow, axis=0)
        x_min     = np.min(x_pow, axis=0)  
        
        ax.plot(x_med, label='median')
        ax.plot(x_min, label='min')
        ax.plot(x_max, label='max')
        plt.minorticks_on()
        plt.xlabel("Frequency")
        plt.ylabel("Amplitude")
        plt.title(self.uv.stokes_axis[0])
        plt.legend()

        ax = plt.subplot(222)
        y_pow     = np.abs(y)
        y_max     = np.max(y_pow, axis=0)  
        y_med     = np.median(y_pow, axis=0)
        y_min     = np.min(y_pow, axis=0)  
        
        ax.plot(y_med, label='median')
        ax.plot(y_min, label='min')
        ax.plot(y_max, label='max')
        plt.minorticks_on()
        plt.xlabel("Frequency")
        plt.ylabel("Amplitude")
        plt.title(self.uv.stokes_axis[1])
        plt.legend()        
        
        ax = plt.subplot(245)
        plt.imshow(x_pow)
        plt.title("Amplitude")
        plt.xlabel("Frequency channel")
        plt.ylabel("Time")
        ax.set_aspect(x.shape[1] / x.shape[0] * 3. / 4)
        plt.colorbar(orientation='horizontal')
        
        ax = plt.subplot(246)
        plt.imshow(np.angle(x))
        plt.title("Phase")
        plt.xlabel("Frequency channel")
        plt.ylabel("Time")
        ax.set_aspect(x.shape[1] / x.shape[0] * 3. / 4)
        cbar = plt.colorbar(orientation='horizontal')
        cbar.set_ticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi-0.05])
        cbar.set_ticklabels(["$-\pi$","$-\pi/2$",0,"$\pi/2$","$\pi$"])
        
        ax = plt.subplot(247)
        plt.imshow(y_pow)
        plt.title("Amplitude")
        plt.xlabel("Frequency channel")
        plt.ylabel("Time")
        ax.set_aspect(x.shape[1] / x.shape[0] * 3. / 4)
        plt.colorbar(orientation='horizontal')
        
        ax = plt.subplot(248)
        plt.imshow(np.angle(y))
        plt.title("Phase")
        plt.xlabel("Frequency channel")
        plt.ylabel("Time")
        ax.set_aspect(x.shape[1] / x.shape[0] * 3. / 4)
        cbar = plt.colorbar(orientation='horizontal')
        cbar.set_ticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi-0.05])
        cbar.set_ticklabels(["$-\pi$","$-\pi/2$",0,"$\pi/2$","$\pi$"])
        
        plt.subplots_adjust(left=0.05, right=0.98, top=0.9, bottom=0.05, wspace=0.25, hspace=0.3)
        
        return fig, ax
    
    def plot_visibilities(self, ref_ant=1, axis=0, n_rows=4, n_cols=8, plot_type='phase'):
        """ Plot visibilities
        
        ref_ant: int
            reference antenna
        axis: int
            which axis (stokes parameter) to view
        plot_type: str
            Plot type. Can be amplitude 'amp', phase 'phase'.
        """
        self.current_plot = 'multi'
        self.ax_zoomed = False
        
        bls = self.uv.d_uv_data['BASELINE']
        
        
        # Extract the relevant baselines using a truth array
        bls = bls.tolist()
        bl_ids = self.uv.get_baseline_ids(ref_ant)
        bl_truths = np.array([(b in bl_ids) for b in bls])
        
        x_data    = self.uv.d_uv_data['DATA'][bl_truths,0,0,0,:,axis]  # Baselines, freq and stokes
        x_cplx    = x_data[:,:,0] + 1j * x_data[:,:,1]
        
        
        # Plot the figure
        fig = self.sp_fig
        
        figtitle = '%s %s: %s -- %s'%(self.uv.telescope, self.uv.instrument, self.uv.source, self.uv.date_obs)
        for i in range(n_rows):
            for j in range(n_cols):
                ax = fig.add_subplot(n_rows, n_cols, i*n_cols + j +1)
                ax.set_title("%s %s"%(ref_ant, i*n_cols + j +1))
                #ax.set_title("%s %s"%(i, j))
                x = x_cplx[i*n_cols+j::self.uv.n_ant]

                #img.set_interpolation('nearest') # Bicubic etc
                #img.set_cmap('jet')               # summer, hot, spectral, YlGnBu
                
                                
                if plot_type == 'phase':
                    img = ax.imshow(np.angle(x), vmin=-np.pi, vmax=np.pi)
                elif plot_type == 'amp':
                    img = ax.imshow(np.abs(x))
                else:
                    print "Error: plot_type %s not understood"%plot_type
                    raise

                if i == n_rows-1:
                    ax.set_xlabel('Freq')
                else: 
                    ax.set_xlabel('')
                if j == 0:
                    ax.set_ylabel('Time')
                else:
                    ax.set_ylabel('')
                ax.set_aspect(x.shape[1] / x.shape[0])
        
        if plot_type == 'phase':
            # Add phase colorbar
            cax = fig.add_axes([0.925,0.08,0.015,0.8])
            cbar = fig.colorbar(img, cax=cax)
            #cbar.set_label('Phase [rad]')
            cbar.set_ticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
            cbar.set_ticklabels(["$-\pi$","$-\pi/2$",0,"$\pi/2$","$\pi$"])
                    
            plt.subplots_adjust(left=0.05, right=0.9, top=0.9, bottom=0.05)
        else:
            plt.tight_layout()
        
        return fig, ax
    
    def plot_autocorrs(self, axis=0, n_rows=4, n_cols=8):
        """ Plot autocorrelations for all antennas
        """
        self.current_plot = 'multi'
        self.ax_zoomed = False
        
        bls = self.uv.d_uv_data['BASELINE']

        # Extract the relevant baselines using a truth array
        bls = bls.tolist()
        bl_ids = [256*i + i for i in range(1,self.uv.n_ant+1)]
        bl_truths = np.array([(b in bl_ids) for b in bls])
        
        #print self.uv.d_uv_data['DATA'].shape
        #x_data    = self.d_uv_data['DATA'][bl_truths,0,0,:,0,axis]  # Baselines, freq and stokes
        x_data    = self.uv.d_uv_data['DATA'][bl_truths, 0,0,0,:,axis,:]
        x_cplx    = x_data[:,:,0] + 1j * x_data[:,:,1]
        
        #print x_cplx.shape
        
        # Plot the figure
        #print self.uv.n_ant
        fig = self.sp_fig
        figtitle = '%s %s: %s -- %s'%(self.uv.telescope, self.uv.instrument, self.uv.source, self.uv.date_obs)
        for i in range(n_rows):
            for j in range(n_cols):
                ax = fig.add_subplot(n_rows, n_cols, i*n_cols + j +1)
                ax.set_title("%s"%( i*n_cols + j +1))
                #ax.set_title("%s %s"%(i, j))
                
                x = x_cplx[i*n_cols+j::self.uv.n_ant]
                
                x_pow     = np.abs(x)
                #x_avg     = np.average(x_pow, axis=0)
                x_max     = np.max(x_pow, axis=0)
                x_med     = np.median(x_pow, axis=0)
                x_min     = np.min(x_pow, axis=0)
                
                ax.plot(x_med)
                ax.plot(x_min)
                ax.plot(x_max)
                
                if i == n_rows-1:
                    ax.set_xlabel('Freq')
                if j == 0:
                    ax.set_ylabel('Amplitude')
                
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                plt.tick_params(axis='both', which='major', labelsize=10)
                plt.tick_params(axis='both', which='minor', labelsize=8)
                plt.xticks(rotation=30)
                #ax.set_aspect(x.shape[1] / x.shape[0])
        
        plt.subplots_adjust(left=0.05, right=0.98, top=0.9, bottom=0.05, wspace=0.25, hspace=0.3)
        return fig, ax
        
    
    def plot_uv(self):
        """ Plot UV points. Just because its easy. """
        self.current_plot = 'uv'
        self.ax_zoomed = False
        
        uu = self.uv.d_uv_data['UU']*1e6
        vv = self.uv.d_uv_data['VV']*1e6
        
        pmax, pmin = np.max([uu, vv])*1.1, np.min([uu,vv])*1.1
        
        fig = self.sp_fig
        plt.subplot(111, aspect='equal')
        ax = plt.plot(uu, vv, 'bo')
       
        plt.xlabel("UU [$\\mu s$]")
        plt.ylabel("VV [$\\mu s$]")
        plt.xlim(pmin, pmax)
        plt.ylim(pmin, pmax)
        return fig, ax
        
        
    def updatePlot(self):
        """ Update plot """
        self.sp_ax.clear()
        self.sp_fig.clear()
        
        self.updateSpinners()

        axis = self.axes_select.currentIndex()
        ref_ant = self.spin_ref_ant.value()
        ref_ant2 = self.spin_ref_ant2.value()
        
        if self.plot_select.currentIndex() == 0:
            self.plot_single_baseline(ref_ant, ref_ant2, axis=axis)
            
        if self.plot_select.currentIndex() == 1:
            self.plot_single_baseline_dual_pol(ref_ant, ref_ant2, axis=axis)
        
        if self.plot_select.currentIndex() == 2:
            self.plot_autocorrs(axis=axis)
        
        if self.plot_select.currentIndex() == 3:
            self.plot_visibilities(plot_type='amp', axis=axis, ref_ant=ref_ant)
        
        if self.plot_select.currentIndex() == 4:
            self.plot_visibilities(plot_type='phase', axis=axis, ref_ant=ref_ant)

        if self.plot_select.currentIndex() == 5:
            self.plot_uv()            
        self.sp_canvas.draw()
        
    def updateSpinners(self):
        """ Update spinner values """
        
        self.spin_ref_ant.setRange(1, self.uv.n_ant)
        self.spin_ref_ant2.setRange(1, self.uv.n_ant)
        
        axis = self.axes_select.currentIndex()
        ref_ant = self.spin_ref_ant.value()
        ref_ant2 = self.spin_ref_ant2.value()
        
        if self.current_plot == 'single':
            self.spin_ref_ant2.setEnabled(True)
        else:
            self.spin_ref_ant2.setDisabled(True)
        
        if self.current_plot == 'multi':
            self.spin_ref_ant.setEnabled(False)
        else:
            self.spin_ref_ant.setEnabled(True)
        
                      
    def updateAxes(self):
        """ Spinner action: update data axes """
        pass
    
    def onFileOpen(self):
        """ Do this whenever a new file is opened """
        # Recreate combobox whenever file is loaded
        c_ind = self.axes_select.currentIndex()
        for i in range(self.axes_select.count()):
            self.axes_select.removeItem(0)            
        for v in self.uv.stokes_axis:
            self.axes_select.addItem(v)
        if self.axes_select.count() <= c_ind:
            self.axes_select.setCurrentIndex(c_ind)
            
    def onButOpen(self):
        """ Button action: Open station file """
        self.file_dialog    = QtGui.QFileDialog()
        fileparts = self.file_dialog.getOpenFileName()
        
        # Compatibility check
        if isinstance(fileparts, str) or isinstance(fileparts, unicode):
            filename = fileparts
        
        elif not USES_PYSIDE:
            if isinstance(fileparts, PyQt4.QtCore.QString):
                filename = str(fileparts)
            else:
                filename = fileparts
        else:
            filename = fileparts[0]
        
        self.uv = InterFits(filename)

        self.onFileOpen()
        self.updatePlot()
    
    def onButStats(self):
        """ Button action: Open statistics """
        pass


def main():
    
    # Basic option parsing 
    p = OptionParser()
    p.set_usage('fits_pattern_viewer.py [filename] [options]')
    p.set_description(__doc__)
    (options, args) = p.parse_args()

    print "Starting %s..."%progname
    global main_gui
    app = QtGui.QApplication(sys.argv)
    
    try:
        filename = args[0]
        main_gui = InterFitsGui(filename)
    except:
        main_gui = InterFitsGui()
    app.exec_()
    sys.exit()    

if __name__ == '__main__':
    main()
