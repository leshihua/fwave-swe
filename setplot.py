
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import os
import math

import matplotlib.pyplot as plt
from clawpack.visclaw import colormaps

import clawpack.clawutil.data as clawdata
import clawpack.geoclaw.surge.plot as surge
import clawpack.geoclaw.geoplot as geoplot

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    amrdata = clawdata.ClawInputData(2)
    amrdata.read(os.path.join(plotdata.outdir,'claw.data'))

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def set_drytol(current_data):
        # The drytol parameter is used in masking land and water and
        # affects what color map is used for cells with small water depth h.
        # The cell will be plotted as dry if h < drytol.
        # The best value to use often depends on the application and can
        # be set here (measured in meters):
        current_data.user['drytol'] = 1.e-2

        # # Square island
        # current_data.user['island_x'] = [-52e3, -52e3, 52e3, 52e3, -52e3]
        # current_data.user['island_y'] = [-52e3, 52e3, 52e3, -52e3, -52e3]

        # Strip island
        current_data.user['island_x'] = [-52e3, -52e3, 52e3, 52e3, -52e3]
        current_data.user['island_y'] = [-102e3, 102e3, 102e3, -102e3, -102e3]

    plotdata.beforeframe = set_drytol

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    # surface_range = 1e-4
    # speed_max = 1.e-5

    surface_range = 1e-2
    speed_max = 1.e-2

    xlimits = [amrdata.lower[0], amrdata.upper[0]]
    ylimits = [amrdata.lower[1], amrdata.upper[1]]
    surface_limits = [-surface_range,surface_range]
    speed_limits = [0.0, speed_max]

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
        gaugenos='all', format_string='ko', add_labels=True)

    def after_axes(cd):
        surge.days_figure_title(cd)
        addgauges(cd)
        plt.plot(cd.user['island_x'], cd.user['island_y'], 'k--')

    def after_axes_slice(cd, y_bounds=[-1, 1]):
        surge.days_figure_title(cd)
        for x in cd.user['island_x']:
            plt.plot([x, x], y_bounds, 'k--')

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    # Surface
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)

    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = after_axes

    surge.add_surface_elevation(plotaxes,bounds=surface_limits)
    surge.add_land(plotaxes)
    # plotaxes.plotitem_dict['surface'].amr_celledges_show = [1,1,1,1]
    # plotaxes.plotitem_dict['land'].amr_celledges_show = [1,1,1,1]

    # Speed
    plotfigure = plotdata.new_plotfigure(name='Speed', figno=1)

    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Currents'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = after_axes

    surge.add_speed(plotaxes,bounds=speed_limits)
    surge.add_land(plotaxes)
    # plotaxes.plotitem_dict['speed'].amr_celledges_show = [1,1,1,1]
    # plotaxes.plotitem_dict['land'].amr_celledges_show = [1,1,1,1]

    # ========================================================================
    #  Water Velocity Components - Entire Gulf
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Velocity Components - Entire Domain',  
                                         figno=2)
    plotfigure.kwargs['figsize'] = (16,6)
    plotfigure.show = False

    # X-Component
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(121)"
    plotaxes.title = 'Velocity, X-Component'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = after_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surge.water_u
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.pcolor_cmin = -speed_limits[1]
    plotitem.pcolor_cmax = speed_limits[1]
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0, 0, 0]

    surge.add_land(plotaxes)

    # Y-Component
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(122)"
    plotaxes.title = 'Velocity, Y-Component'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = after_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = surge.water_v
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.pcolor_cmin = -speed_limits[1]
    plotitem.pcolor_cmax = speed_limits[1]
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0, 0, 0]
    
    surge.add_land(plotaxes)

    # ========================================================================
    #  Water Momenta Components
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Momentum Components - Entire Domain',  
                                         figno=3)
    plotfigure.kwargs['figsize'] = (16,6)
    plotfigure.show = False

    # X-Component
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(121)"
    plotaxes.title = 'Momentum, X-Component'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = after_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 1
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.pcolor_cmin = -speed_limits[1] * 100.0
    plotitem.pcolor_cmax = speed_limits[1] * 100.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0, 0, 0]

    surge.add_land(plotaxes)

    # Y-Component
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(122)"
    plotaxes.title = 'Momentum, Y-Component'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = after_axes

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    plotitem.pcolor_cmin = -speed_limits[1] * 100.0
    plotitem.pcolor_cmax = speed_limits[1] * 100.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0, 0, 0]
    
    surge.add_land(plotaxes)

    # ============
    #  Topography
    # ============
    plotfigure = plotdata.new_plotfigure(name="Topography")
    plotfigure.show = False
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Topography"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.topo
    cmap = colormaps.make_colormap({-1:[0.3,0.2,0.1],
                                       -0.00001:[0.95,0.9,0.7],
                                       0.00001:[.5,.7,0],
                                       1:[.2,.5,.2]})
    plotitem.pcolor_cmap = cmap
    plotitem.pcolor_cmin = -100.0
    plotitem.pcolor_cmax = 0.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0, 0, 0]

    # =====================
    #  Slice through y = 0
    # =====================
    slice_index = 3
    def slice_surface(cd):
        return cd.x[:, slice_index], cd.q[3, :, slice_index]
    
    def slice_momentum(cd):
        return cd.x[:, slice_index], cd.q[1, :, slice_index]

    def slice_velocity(cd):
        u = cd.q[1, :, slice_index] / cd.q[0, :, slice_index]
        return cd.x[:, slice_index], u

    plotfigure = plotdata.new_plotfigure(name="Slice eta, y=0")
    plotfigure.show = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Surface Slice y = 0"
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = surface_range
    plotaxes.afteraxes = lambda cd: after_axes_slice(cd, surface_limits)

    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = slice_surface
    plotitem.plot_var = 3
    plotitem.plotstyle = 'k-o'
    plotitem.show = True

    plotfigure = plotdata.new_plotfigure(name="Slice u, y=0")
    plotfigure.show = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Velocity-x Slice y = 0"
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-speed_limits[1], speed_limits[1]]
    plotaxes.afteraxes = lambda cd: after_axes_slice(cd, [-speed_limits[1], speed_limits[1]])

    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = slice_velocity
    plotitem.plot_var = 1
    plotitem.plotstyle = 'k-o'
    plotitem.show = True

    plotfigure = plotdata.new_plotfigure(name="Slice hu, y=0")
    plotfigure.show = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "hu Slice y = 0"
    plotaxes.xlimits = xlimits
    momentum_limits = [-speed_limits[1] * 1000.0, speed_limits[1] * 1000.0] 
    plotaxes.ylimits = momentum_limits
    plotaxes.afteraxes = lambda cd: after_axes_slice(cd, momentum_limits)

    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = slice_momentum
    plotitem.plot_var = 1
    plotitem.plotstyle = 'k-o'
    plotitem.show = True

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Gauge Surface', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    def gauge_afteraxes(cd,label=None):
        surge.gauge_afteraxes(cd)
        plt.ylabel(label)

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo

    # Surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.afteraxes = lambda cd: gauge_afteraxes(cd,label='Surface Height (m)')
    plotaxes.title = 'Surface'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'
    try:
        plotaxes.xlimits = [amrdata.t0,amrdata.tfinal]
    except:
        pass

    # Momenta
    plotfigure = plotdata.new_plotfigure(name='Gauge X-Momentum', figno=301, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.afteraxes = lambda cd: gauge_afteraxes(cd,label=r'X-Momentum ($m/s^2$)')
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'X Momentum'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = True
    plotitem.plot_var = 1
    try:
        plotaxes.xlimits = [amrdata.t0,amrdata.tfinal]
    except:
        pass

    plotfigure = plotdata.new_plotfigure(name='Gauge Y-Momentum', figno=302, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.afteraxes = lambda cd: gauge_afteraxes(cd,label=r'Y-Momentum ($m/s^2$)')
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Y Momentum'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = True
    plotitem.plot_var = 2
    try:
        plotaxes.xlimits = [amrdata.t0,amrdata.tfinal]
    except:
        pass

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = range(0, 36, 4)# list of frames to print
    plotdata.print_gaugenos = [1, 2, 3]  # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.format = 'ascii'                # Format of output
    # plotdata.format = 'netcdf'             

    return plotdata

    
