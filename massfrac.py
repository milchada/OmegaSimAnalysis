#omega plots

import numpy as np
import matplotlib.pylab as plt
from matplotlib.pyplot import cm 
import glob
import os
import mass_headers as m
import gastherm as gt

sciencedir = '/lustre/scratch/client/fas/nagai/projects/L100_0'
homedir = '/home/fas/nagai/uc24'
os.chdir(sciencedir)
# runs = ['/CSF', '/AGN/CCA/recentering', '/AGN/CCA/largescalethermal', '/AGN/Bondi/recentering', '/AGN/Kick/fiducial']#'/CCA/nocool_lowTcrit','/Bondi/fiducial','/Bondi/eff_0.1','/Bondi/eff_0.3','/Bondi/sf_den_9'
runs = ['/CSF','/AGN/CCA/recentering/fiducial','/AGN/CCA/largethermal/fiducial', '/AGN/CCA/largethermal/rfb10', '/AGN/CCA/largethermal/rfb50','/AGN/CCA/largethermal/slope1','/AGN/CCA/largethermal/slope3']
runname = ['No AGN', r'$r_{FB}$=3, $\alpha$=2', r'$r_{FB}$=25, $\alpha$=2', r'$r_{FB}$=10, $\alpha$=2', r'$r_{FB}$=50, $\alpha$=2', r'$r_{FB}$=25, $\alpha$=1', r'$r_{FB}$=25, $\alpha$=3']
#but if you could just print r_K in kpc we could convert r_fb to kpc and that's SUPER USEFUL

properties = ['gas_frac', 'star_frac', 'm_cum']
propname = [r'$f_{gas}(<r)$', r'$f_*$(<r)', r'M(<r)($10^{14}M_\odot$)']

samplefile = glob.glob(sciencedir+runs[3]+'/profile_analysis/profiles/0.6257/*profile_radius_bin_a*')[0]
radii = np.genfromtxt(samplefile)[:,0]	

def read_run(snapshot, prop='gas'):
	gasfile = glob.glob(snapshot+'/*profile_'+prop+'_a*')[0]
	gas = np.genfromtxt(gasfile)
	return gas

def mass_fractions(snapshot):
	mass = read_run(snapshot, prop='mass')
	gas_index = m.profile_columns['gas_M_cum']
	dm_index = m.profile_columns['dark_M_cum']
	star_index = m.profile_columns['star_M_cum']
	total_mass = mass[:,gas_index]+mass[:,dm_index]+mass[:,star_index]
	gas_frac = np.divide(mass[:,gas_index], total_mass)
	star_frac = np.divide(mass[:,star_index], total_mass)
	return gas_frac, star_frac, total_mass

def closest_snap_run(a, rundir=runs[2]):
	rundir += '/profile_analysis/profiles/*'
	# print rundir
	snaps = sorted(glob.glob(os.getcwd()+rundir))
	# print snaps
	a_s = np.asarray([float(snap.split('/')[-1]) for snap in snaps])
	closest_snap = snaps[np.where(abs(a_s - a) == min(abs(a_s - a)))[0]]
	return closest_snap

def applyPlotStyle(ax, ylabel, title, snapshot):
	ax.set_xlabel('R (kpc)')
	ax.set_title(title)
	ax.set_xscale('log')
	ax.set_xlim(10,1000)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels, loc='best')
	ax.ylabel(ylabel)
	fine_radii, delta_r = gt.delta_r(snapshot)
	r500 = gt.Rdelta(delta_r, 500, fine_radii)
	r2500 = gt.Rdelta(delta_r, 2500, fine_radii)
	ymin, ymax = ax.get_ylim()
	ax.vlines(r500, ymin, ymax, linestyle='dashed', label=r'$R_{500c}$')
	ax.vlines(r2500, ymin, ymax, linestyle='dashed', label=r'$R_{2500c}$')

def compare_runs_a(a=None):
	fig1, ax1 = plt.subplots()
	fig2, ax2 = plt.subplots()
	colors = cm.rainbow(np.linspace(0,1, len(runs)))
	for run, color in zip(runs, colors):
		print run.split('/')
		last_snapshot = closest_snap_run(a, run)
		gas_frac, mass_frac = mass_fractions(last_snapshot)
		label = runname[runs==run]
		ax1.plot(radii, gas_frac, color=color,label=label)
		ax1.hlines(0.15, 10,1000,linestyle='dashed',label=r'$f_{b,cosmic}$')
		ax2.plot(radii, star_frac, color=color, label=label)
	title = 'a = '+str(a)+' '+property
	ylabel_gas = propname[0]+' ('+units[0]+')'
	applyPlotStyle(ax1, ylabel_gas, title)
	ylabel_star = propname[1]+' ('+units[1]+')'
	applyPlotStyle(ax2, ylabel_star, title)
	fig1.savefig(homedir+'fgas_'+str(a)+'.png')
	fig2.savefig(homedir+'fstar_'+str(a)+'.png')

def plot_evolution(property):
	for run in runs:
		rundir = os.getcwd()+'/'+run+'/profile_analysis/profiles/*'
		snapshots = sorted(glob.glob(rundir)) #this gives a list of paths to each snapshot
		fig1, ax1 = plt.subplots()
		colors = cm.rainbow(np.linspace(0,1, len(snapshots)))
		plt.rcParams['legend.fontsize']= 'small'
		for snapshot, color in zip(snapshots,colors):
			snapname = snapshot.split('/')[-1]
			print run, snapname
			try:
				if property in properties:
					index = properties.index(property)
					profiles = mass_fractions(snapshot)
					profile = profiles[index]
				else:
					profile = gt.compute_profile(property, snapshot)
				lines=ax1.plot(radii, profile, c=color,label=snapname[0:4])
			except (KeyError, ValueError, IndexError) as e:
				print 'oops!'
				continue
		ax1.set_title(property+' in '+run)
		if (property != 'gas_frac') and (property != 'star_frac'):
			ax1.set_yscale('log')
		ax1.set_xscale('log')
		ax1.set_xlim(10,1000)
		handles, labels = ax1.get_legend_handles_labels()
		ax1.legend(handles, labels)
		ax1.set_ylabel(propname[properties==property]+units[properties.index(property)])
		ax1.set_xlabel('R (kpc)')
		# fig1.colorbar(lines)
		label = runname[runs==run]
		fig1.savefig(homedir+'/gas_'+property+'_evoln_'+label+'.png')

def core_evolution(property):
	colors = cm.rainbow(np.linspace(0,1, len(runs)))
	plt.rcParams['legend.fontsize']= 'small'
	plt.clf()
	for run, color in zip(runs, colors):
		rundir = os.getcwd()+run+'/profile_analysis/profiles/*'
		snapshots = sorted(glob.glob(rundir)) #this gives a list of paths to each snapshot
		profile = []
		snaps =[]
		for snapshot in snapshots:	
			try:
				if property == 'gas_frac':
					snap_profile, star_frac = mass_fractions(snapshot)
				else:
					gas_frac, snap_profile = mass_fractions(snapshot)
				profile.append(snap_profile[0])
				snaps.append(snapshot.split('/')[-1])
			except (KeyError, ValueError, IndexError) as e:
				# print run, snapshot, 'oops!'
				continue
		label = runname[runs==run]
		plt.plot(snaps, profile, c=color, label = label)
		print run, ' done!'
	plt.ylabel(property)
	plt.legend(loc='best')
	plt.title('Evolution of Central Bin')
	plt.xlabel('Scale Factor')
	plt.savefig(homedir+'/'+property+'_in_core_allruns.png')
	#add merger lines here?


if __name__ == '__main__':
	a_latest = float(input("Snapshot at which to compare runs: "))
	compare_runs_a(property, a=a_latest)
	# for property in properties:
	# 	core_evolution(property)
