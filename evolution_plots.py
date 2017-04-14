#omega plots

import numpy as np
import matplotlib.pylab as plt
from matplotlib.pyplot import cm 
import glob
import os
import mass_headers as m

sciencedir = '/lustre/scratch/client/fas/nagai/projects/L100_0'
homedir = '/home/fas/nagai/uc24'
os.chdir(sciencedir)
# runs = ['/CSF', '/AGN/CCA/recentering', '/AGN/CCA/largescalethermal', '/AGN/Bondi/recentering', '/AGN/Kick/fiducial']#'/CCA/nocool_lowTcrit','/Bondi/fiducial','/Bondi/eff_0.1','/Bondi/eff_0.3','/Bondi/sf_den_9'
runs = ['/CSF','/AGN/CCA/recentering/fiducial','/AGN/CCA/largethermal/fiducial', '/AGN/CCA/largethermal/rfb10', '/AGN/CCA/largethermal/rfb50','/AGN/CCA/largethermal/slope1','/AGN/CCA/largethermal/slope3']
runname = ['No AGN', r'$r_{FB}$=3, $\alpha$=2', r'$r_{FB}$=25, $\alpha$=2', r'$r_{FB}$=10, $\alpha$=2', r'$r_{FB}$=50, $\alpha$=2', r'$r_{FB}$=25, $\alpha$=1', r'$r_{FB}$=25, $\alpha$=3']
#but if you could just print r_K in kpc we could convert r_fb to kpc and that's SUPER USEFUL

keys = ['gas_temperature(mass-weighted) [K','gas_entropy(mass-weighted) [keV cm ^ 2', 'gas_density(emission-weighted) [g cm^{-3}', 'gas X-ray emission [erg /s']
properties = ['temperature','entropy', 'density', 'gas_frac', 'star_frac']
propname = ['T', 'K', r'$\rho$', r'$f_{gas}(<r)$', r'$f_*$(<r)']
units = ['K', r'keV cm$^{-2}$', r'g cm$^{-3}$','','']

samplefile = glob.glob(sciencedir+runs[3]+'/profile_analysis/profiles/0.6257/*profile_radius_bin_a*')[0]
radii = np.genfromtxt(samplefile)[:,0]	

def read_run(snapshot, prop='gas'):
	gasfile = glob.glob(snapshot+'/*profile_'+prop+'_a*')[0]
	gas = np.genfromtxt(gasfile)
	return gas

def compute_profile(property, snap):
	gas = read_run(snap)

	with open(homedir+'/gas_colnames.txt','r') as file:
	    colnames=file.readlines()[1]
	gas_colnames = colnames.split('] ')
	
	propkey = [key for key in keys if property in key][0]
	# print gas_colnames.index(propkey), gas.keys()
	return gas[:,gas_colnames.index(propkey)]

def mass_fractions(snapshot):
	mass = read_run(snapshot, prop='mass')
	gas_index = m.profile_columns['gas_M_cum']
	dm_index = m.profile_columns['dark_M_cum']
	star_index = m.profile_columns['star_M_cum']
	total_mass = mass[:,gas_index]+mass[:,dm_index]+mass[:,star_index]
	gas_frac = np.divide(mass[:,gas_index], total_mass)
	star_frac = np.divide(mass[:,star_index], total_mass)
	return gas_frac, star_frac

def closest_snap_run(a, rundir=runs[2]):
	rundir += '/profile_analysis/profiles/*'
	# print rundir
	snaps = sorted(glob.glob(os.getcwd()+rundir))
	# print snaps
	a_s = np.asarray([float(snap.split('/')[-1]) for snap in snaps])
	closest_snap = snaps[np.where(abs(a_s - a) == min(abs(a_s - a)))[0]]
	return closest_snap

def compare_runs_a( property,a=None):
	fig1 = plt.figure()
	fig2 = plt.figure()
	ax1 = fig1.add_subplot(111)
	ax2 = fig2.add_subplot(111)
	# print last_fid_snapshot_cca, last_fid_snapshot_bondi
	colors = cm.rainbow(np.linspace(0,1, len(runs)))
	for run, color in zip(runs, colors):
		print run.split('/')
		if a== None:
			last_snapshot = sorted(glob.glob(os.getcwd()+run+'/profile_analysis/profiles/*'))[-1]
		else: 
			last_snapshot = closest_snap_run(a, run)
		if property == 'gas_frac':
			profile, mass_frac = mass_fractions(last_snapshot)
		elif property == 'star_frac':
			gas_frac, profile = mass_fractions(last_snapshot)
		else:
			profile = compute_profile(property,last_snapshot)
		label = runname[runs==run]
		ax1.plot(radii, profile,color=color,label=label)
	if a != None:
		title = 'a = '+str(a)+' '+property
		ax1.set_title(title)
	else:
		ax1.set_title('a = 1 '+property)
	if (property != 'gas_frac') and (property != 'star_frac'):
		ax1.set_yscale('log')
	ax1.set_xscale('log')
	ax1.set_xlim(10,1000)
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels, loc='best')
	ax1.set_ylabel(propname[properties==property]+units[properties.index(property)])
	ax1.set_xlabel('R (kpc)')
	fig1.savefig(homedir+property+'.png')

def plot_evolution():
	for property in properties:
		for run in runs:
			rundir = os.getcwd()+'/'+run+'/profile_analysis/profiles/*'
			snapshots = sorted(glob.glob(rundir)) #this gives a list of paths to each snapshot
			fig1 = plt.figure()
			ax1 = fig1.add_subplot(111)
			colors = cm.rainbow(np.linspace(0,1, len(snapshots)))
			plt.rcParams['legend.fontsize']= 'small'
			for snapshot, color in zip(snapshots,colors):
				snapname = snapshot.split('/')[-1]
				print run, snapname
				try:
					if property == 'gas_frac':
						profile, star_frac = mass_fractions(snapshot)
					elif property == 'star_frac':
						gas_frac, profile = mass_fractions(snapshot)
					else:
						profile = compute_profile(property, snapshot)
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
				elif property == 'star_frac':
					gas_frac, snap_profile = mass_fractions(snapshot)
				else:
					snap_profile = compute_profile(property,snapshot)
				profile.append(snap_profile[0])
				snaps.append(snapshot.split('/')[-1])
			except (KeyError, ValueError, IndexError) as e:
				# print run, snapshot, 'oops!'
				continue
		label = runname[runs==run]
		plt.plot(snaps, profile, c=color, label = label)
		print run, ' done!'
	plt.ylabel(property)
	if (property != 'gas_frac') and (property != 'star_frac'):
		plt.yscale('log')
	plt.legend(loc='best')
	plt.title('Evolution of Central Bin')
	plt.xlabel('Scale Factor')
	plt.savefig(homedir+'/'+property+'_in_core_allruns.png')


if __name__ == '__main__':
	a_latest = float(input("Snapshot at which to compare runs: "))
	for property in properties:
		compare_runs_a(property, a=a_latest)
		core_evolution(property)
	plot_evolution()
