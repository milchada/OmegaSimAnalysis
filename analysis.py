
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 
import yt 
from yt import derived_field
import sys
from yt.units import G
import random

@derived_field(name="t_ff", units="s")
def _t_ff(field, data) :
    densities = data[('Density')]
    radii = data[('radius')]
    enclosed_density = [np.mean(densities[radii<radii[ind]]) for ind in range(len(radii))]
    return np.sqrt(3*np.pi/(32*G*enclosed_density)) #this is the cell density, need the mean density within r
    """
   YT should have a function to compute enclosed mass/density. Find it.
    """

@derived_field(name="tcool_tff", units="dimensionless")
def _tctff(field, data):
    ratio = data["Cooling_Time"]/data[("t_ff")]
    return ratio

ts = yt.load("/home/fas/nagai/etl28/group_scratch/Yuan/DD????/stest_????")

# Initialization
i = 0

# Loop over timesteps
def plot_timesteps():
    i = 0
    for ds in ts[1::5]:
        i = i + 1
        #   Only do a couple steps for debugging
        # if i > 2:
        #    break
        print i

        # Select a spherical region with radius of 2000 kpc centered in the simulation box. 
        # ds.add_field('t_ff',function=_t_ff, units=ds.unit_system['time'],force_override=True)
        sp = ds.sphere(ds.domain_center, (2000, "kpc"))
        # Coarsening of the grid

        # Phase diagram plot
        nbins = 2000/3.8 #total radius divided by our resolution
        profile = yt.create_profile( data_source=sp, bin_fields=["radius","tcool_tff"], fields=["cell_mass"], units=dict(radius="kpc",cell_mass="Msun"), logs=dict(radius=True,tcool_tff=True),weight_field=None, n_bins=nbins, extrema={'radius':(1,2000),'tcool_tff':(1e-5, 1e5)})
        plot = yt.PhasePlot.from_profile(profile)
        plot.set_cmap("cell_mass", "YlOrRd")
        prefix = "%s" % (profile.ds)
        plot.save(prefix+'_density_tcool_tff.png')

        profile = yt.create_profile( data_source=sp, bin_fields=["radius","cooling_time"], fields=["cell_mass"], units=dict(radius="kpc",cell_mass="Msun", cooling_time='Gyr'), logs=dict(radius=True,cooling_time=True),weight_field=None, n_bins=nbins, extrema={'radius':(1,2000),'cooling_time':(1e-7, 1e5)})
        plot = yt.PhasePlot.from_profile(profile)
        plot.set_cmap("cell_mass", "YlOrRd")
        prefix = "%s" % (profile.ds)
        plot.save(prefix+'_density_tcool.png')

        profile = yt.create_profile( data_source=sp, bin_fields=["radius","t_ff"], fields=["cell_mass"], units=dict(radius="kpc",cell_mass="Msun", t_ff='Gyr'), logs=dict(radius=True,t_ff=True),weight_field=None, n_bins=nbins, extrema={'radius':(1,2000),'t_ff':(1e-3, 1e2)})
        plot = yt.PhasePlot.from_profile(profile)
        plot.set_cmap("cell_mass", "YlOrRd")
        prefix = "%s" % (profile.ds)
        plot.save(prefix+'_density_tff.png')

        profile = yt.create_profile( data_source=sp, bin_fields=["radius","temperature"], fields=["cell_mass"], units=dict(radius="kpc",cell_mass="Msun"), logs=dict(radius=True,temperature=True),weight_field=None, n_bins=nbins, extrema={'radius':(1,2000),'temperature':(1e2, 1e11)})
        plot = yt.PhasePlot.from_profile(profile)
        plot.set_cmap("cell_mass", "YlOrRd")
        prefix = "%s" % (profile.ds)
        plot.save(prefix+'_density_temp.png')


def plot_min_precip_ratio():
    min_tcool_tff = []
    timestep = []
    for ds in ts[1::100]:
        # Select a spherical region representing the cluster center - how did Yuan do this in her paper?
        sp = ds.sphere(ds.domain_center, (20, "kpc"))
        timestep.append(ds.current_time.in_units("Gyr"))
        min_tcool_tff.append(min(sp["tcool_tff"]))
        print len(sp['tcool_tff'])

    plt.plot(timestep, min_tcool_tff)
    plt.xlabel('Time (Gyr)')
    plt.ylabel(r'Minimum $t_{cool}/t_{ff}$ (s)')
    plt.savefig('cooling_ratio_evolution.png')


def temp_cooling_ratio(precip_threshold=10, nmin=10):
    tcool_tff = []
    temp = []
    for ds in ts[1::5]:
        sp = ds.sphere(ds.domain_center,(100,"kpc"))
        cooling_ratio = sp['tcool_tff']
        temperature = sp['temperature']
        timestep = ds.current_time.value.min()
        nparticles =  len(cooling_ratio[cooling_ratio<precip_threshold])
        if nparticles>0:
            min_particles = min(nmin,nparticles)
            random.shuffle(cooling_ratio[cooling_ratio<precip_threshold])
            random.shuffle(temperature[cooling_ratio<precip_threshold])
            plt.scatter(temperature[:min_particles], cooling_ratio[:min_particles], lw=0, s=0.5, cmap=plt.cm.Reds)
            tcool_tff.append(cooling_ratio[:min_particles])
            temp.append(temperature[:min_particles])
            print "Finished for timestep ", timestep
    plt.savefig('temperature_coolingratio.png')
    return tcool_tff, temp
