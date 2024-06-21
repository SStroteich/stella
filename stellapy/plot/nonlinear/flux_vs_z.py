"""

#===============================================================================
#                                 Plot phi2(z)                                 #
#===============================================================================

Based on <folder> create a <research> and create a figure with nspec subplots:
        - ax * nspec: flux_s(z)

Arguments
--------- 
    folder : pathlib.Path directory where the command has been executed
    x_quantity : {z, pol, tor}
    y_quantity : {phi2, phi, phi_real, phi_imag}
    geometry : {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0} overlayed on the background 
    folderIsExperiment : {False, True} 
    interpolate : {int, False} where <int> (e.g. 20) is the interpolation step
    Choose between {qflux, pflux, vflux}: plot_flux_spectra --qflux (or --pflux or --vflux) 

Hanne Thienpondt 
16/12/2022

"""

#!/usr/bin/python3   
import pathlib
import sys, os
import numpy as np
import matplotlib as mpl  
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.plot.utils.style.get_styleForLinesAndMarkers import get_styleForLinesAndMarkers 
from stellapy.data.input.read_inputFile import read_nonlinearFullFluxSurfaceFromInputFile
from stellapy.plot.utils.labels.get_timeFrameString import get_timeFrameString
from stellapy.plot.utils.species.recognize_species import recognize_species
from stellapy.utils.files.get_firstInputFile import get_firstInputFile
from stellapy.plot.utils.labels.standardLabels import standardLabels  
from stellapy.plot.utils.style.create_figure import create_figure 
from stellapy.utils.decorators.exit_program import exit_program
from stellapy.simulations.Research import create_research 
from stellapy.plot.utils.style import Axis, Legend, Plot
from stellapy.utils.commandprompt.bash import Bash

#===============================================================================
#                                 Plot phi2(z)                                 #
#===============================================================================

def plot_flux_vs_z(
        folder, 
        # Quantities to be plotted
        x_quantity="z", 
        y_quantity="qflux", 
        specie=None, 
        # Overlay a geometric quantity in gray
        geometry="bmag", 
        # Plotting options
        log=False,
        color=None, 
        interpolate=20, 
        normalize_to_one=False,
        folderIsExperiment=False):        
    
    # Create a <research> based on the given <folder>
    try: research = create_research(folders=folder, folderIsExperiment=folderIsExperiment)
    except: research = create_research(folders=folder, folderIsExperiment=True)
    simulation = research.experiments[0].simulations[0]
    nspec = simulation.input.nspec if specie==None else 1
    species = range(simulation.input.nspec) if specie==None else [specie]

    nk = 1

    # Create a figure
    if y_quantity=="qflux": title = "Parallel mode structure of the heat flux spectra"
    elif y_quantity=="pflux": title = "Parallel mode structure of the particle flux spectra"
    elif y_quantity=="vflux": title = "Parallel mode structure of the momentum flux spectra"
    else: exit_program("<y_quantity> = "+y_quantity+" is not a valid choice.", plot_flux_vs_z, sys._getframe().f_lineno) 
    fig = plt.figure(figsize=(18, 9)); axes = []
    grid_specifications = gridspec.GridSpec(nk, nspec)
    if nspec==1: grid_specifications.update(top=0.93, left=0.08, right=0.97, bottom=0.1, wspace=0.2, hspace=0.3)
    if nspec==2: grid_specifications.update(top=0.93, left=0.08, right=0.97, bottom=0.1, wspace=0.2, hspace=0.3)
    if nspec==3: grid_specifications.update(top=0.92, left=0.08, right=0.97, bottom=0.1, wspace=0.3, hspace=0.3)
    for i in range(nspec*nk): axes.append(plt.subplot(grid_specifications[i]))
    update_figure_style(fig, axes)  
    fig.suptitle(title)


    # Add the data to the plot
    for i, specie in enumerate(species):
        subplot_flux_vs_z(axes[i], research, x_quantity, y_quantity, specie=specie, log=log, color=color, geometry=geometry, normalize_to_one=normalize_to_one, interpolate=interpolate)
    
    # Appearance
    if normalize_to_one: ax.set_ylim([0,1]) 
    ax.yaxis.labelpad = 15
    
    # Show the figure   
    mpl.rcParams["savefig.directory"] = folder
    if len(ax.get_lines())!=0: plt.show() 
    return

#-----------------------------
def subplot_flux_vs_z(
        ax, research, 
        # Quantities to be plotted
        x_quantity="z", 
        y_quantity="qflux", 
        specie=0, 
        # Overlay a geometric quantity in gray
        geometry="bmag", 
        # Plotting options
        log=False,
        color=None, 
        interpolate=20, 
        normalize_to_one=False):
    
    # Automate the axis limits and legend
    plot = Plot(); plot.update_legend(loc="upper right");  legend = Legend(ax, plot) 
    axis = Axis(ax, plot, xbot_pos=0, ytop_neg=0, ybot_pos=0, overshoot_y=1.02, logy=log) 
    plot.process_plottingVariables(research); maximum = np.nan; minimum = np.nan;
    tstarts = [np.nan,np.nan]; tends = [np.nan,np.nan] 
    
    # Check whether <specie> is a valid choice 
    if specie >= research.experiments[0].simulations[0].dim.species:
        dim_species = research.experiments[0].simulations[0].dim.species
        exit_reason = "Error: specie = "+str(specie)+" was selected but only "+str(dim_species)+" species are present in the simulation.\n"
        exit_reason += "Please choose specie from {"+", ".join([str(s) for s in range(dim_species)])+"} instead.\n"
        exit_program(exit_reason, subplot_flux_vs_wavenumber, sys._getframe().f_lineno)

    # Iterate over the experiments and simulations
    for experiment in research.experiments:
        for simulation in experiment.simulations:
            
            # Style of the data
            style = get_styleForLinesAndMarkers(plot, legend, research, experiment, simulation)
            style["color"] = color if color!=None else style["color"]
            
            # Read potential(z)
            x = get_field_line_coordinate(x_quantity, simulation)
            y, tstarts, tends = get_quantity_versus_z(y_quantity, simulation, specie, tstarts, tends)
            
            # Interpolate the y-data and normalize the data to have max(y)=1 
            if interpolate: x, y = interpolate_data(x, y, interpolate) 
            if normalize_to_one: y = y/np.max(np.abs(y))

            # Plot potential(z)
            ax.plot(x, y, **style)
            
            # Keep track of the axis limits
            axis.update_axisLimits(x, y)
            maximum = np.nanmax([maximum, np.nanmax(y)])
            minimum = np.nanmin([minimum, np.nanmin(y)])
    
    # Axis labels
    ax.set_xlabel(standardLabels["normalized"][x_quantity]) 

    specie_label = recognize_species(research, specie) 
    extra = "$\\sum_{k_x}$" if x_quantity=="ky" else "$\\sum_{k_y}$"
    ylabel = standardLabels["normalized"][y_quantity].replace("{s}",specie_label)
    ylabel = extra+"$\\langle$" + ylabel + "$\\rangle_{z, t="+get_timeFrameString(tstart, tend)+"}$"
    ax.set_ylabel(ylabel)         
    # Automatically set the axis limits and legends 
    legend.add_legend()
    axis.rescale_axis()
            
    # Plot geometric quantity
    plot_geometric_quantity(ax, x_quantity, geometry, simulation, maximum, minimum, interpolate, plot)
    return


#-----------------------------
def interpolate_data(x, y, interpolate):
    xnew = np.linspace(x[0], x[-1], num=len(x)*interpolate, endpoint=True) 
    f = interp1d(x, y, kind='cubic')
    y = f(xnew)
    x = xnew
    return x, y

#-----------------------------
def add_extra_legend(ax, color, label, plot):
    fake_line = [Line2D([0], [0], color=color, linewidth=3, linestyle="-")]
    extra_legend = ax.legend(fake_line, [label], labelspacing=0.0, shadow=True, \
        loc="upper left", prop={'size':plot.fontsize},handlelength=plot.handlelength)  
    ax.add_artist(extra_legend) 
    return 

#-----------------------------
def plot_geometric_quantity(ax, x_quantity, geometry, simulation, maximum, minimum, interpolate, plot):
    if geometry=="bmag":
        x = get_field_line_coordinate(x_quantity, simulation)
        bmag = simulation.geometry.bmag[:,0]
        bmag = (bmag-np.min(bmag))/np.max(bmag-np.min(bmag))*(maximum-minimum)+minimum
        if interpolate: x, bmag = interpolate_data(x, bmag, interpolate) 
        ax.plot(x, bmag, color="gray")
        add_extra_legend(ax, "gray", "$B$", plot) 
    elif geometry in ["alpha", "zed", "gradpar", "gds2", "gds21", "gds22", "gds23","gds24", "cvdrift", "gbdrift0", "bmag_psi0"]: 
        x = get_field_line_coordinate(x_quantity, simulation)
        geometry_data = getattr(simulation.geometry, geometry)
        geometry_data = geometry_data/np.max(np.abs(geometry_data))*maximum
        if interpolate: x, geometry_data = interpolate_data(x, geometry_data, interpolate) 
        ax.plot(x, geometry_data, color="gray") 
        add_extra_legend(ax, "gray", geometry+'/max('+geometry+')', plot) 
    elif geometry!=None and geometry!=False:
        print("WARNING: <geometry> = "+geometry+" is not recognized.") 
    return 
         
#-----------------------------
def get_quantity_versus_z(y_quantity, simulation, specie, tstarts, tends):   
        
    # Get the quantity versus (t,z)
    if y_quantity=="qflux": 
        vec_quantity = simulation.fluxes.qflux_vs_tsz.qflux[:,specie,:]; t 
        vec_time = simulation.fluxes.qflux_vs_tsz.t
    if y_quantity=="pflux":
        vec_quantity = simulation.fluxes.pflux_vs_tsz.pflux[:,specie,:]
        vec_time = simulation.fluxes.pflux_vs_tsz.t
    if y_quantity=="vflux":
        vec_quantity = simulation.fluxes.vflux_vs_tsz.vflux[:,specie,:]
        vec_time = simulation.fluxes.vflux_vs_tsz.t
        
    # Average over the saturated time frame
    tstart, tend = simulation.time.tstart, simulation.time.tend
    vec_quantity = np.mean(vec_quantity[(vec_time>=tstart)&(vec_time<=tend),:], axis=0)  
        
    # Keep track of the time frames
    tstarts[0] = np.nanmin([tstarts[0], tstart]) 
    tstarts[1] = np.nanmax([tstarts[1], tstart]) 
    tends[0] = np.nanmin([tends[0], tend]) 
    tends[1] = np.nanmax([tends[1], tend]) 
    return vec_quantity, tstarts, tends

#-----------------------------
def get_field_line_coordinate(x_quantity, simulation):
    " Get the quantity on the z-axis"
    if x_quantity == "z":     vec_z = simulation.vec.z 
    if x_quantity == "zeta":  vec_z = simulation.vec.zeta
    if x_quantity == "pol":   vec_z = simulation.vec.pol
    if x_quantity == "tor":   vec_z = simulation.vec.tor
    return vec_z

#===============================================================================
#                             RUN AS BASH COMMAND                              #
#===============================================================================
 
if __name__ == "__main__":
    
    # Create a bash-like interface
    # Toggles are defined through (name, value, shortoption, longoption, explanation)
    # Options are defined through (name, datatype, shortoption, longoption, explanation)
    bash = Bash(plot_flux_vs_z, __doc__)   
    
    # Check which script to launch      
    input_file = get_firstInputFile(pathlib.Path(os.getcwd())) 
    nonlinear, full_flux_surface = read_nonlinearFullFluxSurfaceFromInputFile(input_file)
    if not nonlinear: os.system("python3 $STELLAPY/plot/linear/potential_vs_z.py "+" ".join(sys.argv[1:])); sys.exit()

    # Data 
    bash.add_option('geometry', 'str', 'g', '', 'Plot geometric quantity from {bmag, gradpar, gds2, gds21, gds22, gds23, gds24, cvdrift, gbdrift0, bmag_psi0}.')  
    
    # Data toggles for the y-quantities
    bash.add_toggleheader("y-quantity")
    bash.add_toggle('y_quantity', 'qflux', '', '', 'Plot the kx- and ky-spectra of the heat flux.') 
    bash.add_toggle('y_quantity', 'pflux', '', '', 'Plot the kx- and ky-spectra of the particle flux.') 
    bash.add_toggle('y_quantity', 'vflux', '', '', 'Plot the kx- and ky-spectra of the momentum flux.') 
    bash.add_togglespace()

    # Species
    bash.add_toggleheader("specie")
    bash.add_option('specie', 'int', 's', '', 'Select the species as "-s 0" for ions and "-s 1" for electrons.') 
    bash.add_toggle('specie', 0, '', 'ions', 'Plot the fluxes for the ions.') 
    bash.add_toggle('specie', 1, '', 'electrons', 'Plot the fluxes for the electrons (assuming s=1 for electrons).') 
    bash.add_toggle('specie', 2, '', 'impurity', 'Plot the fluxes for the impurities (assuming s=2 for impurities).')
    bash.add_togglespace() 
    
    # Plotting options 
    bash.add_toggle('normalize_to_one', True, 'n', 'normalize_to_one', 'Normalize the data to a maximum of one.')  
    bash.add_toggle('log', True, 'l', 'log', 'Use logaritmic scales.')  
    
    # Research options  
    bash.add_toggleheader("other")    
    bash.add_option('folderIsExperiment', True, 'f', 'folderIsExperiment', 'Each folder is an experiment.')  
    
    args = bash.get_arguments() 

    # Get the arguments and execute the script
    plot_flux_vs_z(**args)    





    




