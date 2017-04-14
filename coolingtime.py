import numpy as np
import matplotlib.pylab as plt

def cooling_time(snapshot):
	emission_profile = compute_profile('emission', snapshot) #erg/s
	temp_profile = compute_profile('temperature', snapshot)
	internal_energy = temp_profile*1.5*1.38e-16 #K * erg/K
	return np.divide(internal_energy,emission_profile) #s

def freefall_time(snapshot):
	density = compute_profile('density', snapshot)
	return 66430*1/np.sqrt(density) #s

def tcool_tff():

	# fiducial_snapshots = sorted(glob.glob(os.getcwd()+'/fiducial/profile_analysis/profiles/*'))
	# fiducial_snapnums = np.asarray([float(elt[-6:-1]) for elt in fiducial_snapshots])
	
	for run in runs:
		rundir = os.getcwd()+'/'+run+'/profile_analysis/profiles/*'
		snapshots = sorted(glob.glob(rundir)) #this gives a list of paths to each snapshot
		fig = plt.figure()
		ax = fig.add_subplot(111)
		colors = cm.rainbow(np.linspace(0,1, len(snapshots)))
		plt.rcParams['legend.fontsize']= 'small'
		for snapshot, color in zip(snapshots,colors):	
			try:
				snapname = snapshot.split('/')[-1]
				print run, snapname
				tcool = cooling_time(snapshot)
				tff = freefall_time(snapshot)
				profile = np.divide(tcool, tff)
				ax.plot(radii, profile, c=color,label=snapname[0:4])
			except (KeyError, ValueError, IndexError) as e:
				print 'oops!'
				continue
		
		ax.set_title(r'$t_{cool}/t_{ff}$ in '+run)
		ax.set_xscale('log')
		ax.set_xlim(10,1000)
		handles, labels = ax.get_legend_handles_labels()
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(handles, labels, bbox_to_anchor=(1.2,1.1))
		ax.set_xlabel('R (kpc)')
		ax.set_ylabel(r'Cooling ratio ')
		fig.savefig('cooling_ratio_'+run+'evoln.png')
