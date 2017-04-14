#self-similarity
#plot all quantities normalized to values as R_200c
#for comparison with McDonald et al
from normalize_properties import *

def Prop_norm_delta(property='entropy', snapshot, delta=200):
	massfile = glob.glob(snapshot+'/halo_profile_ma*')[0]
	with open(massfile,'r') as file:
		header = file.readlines()[10]
	mvir = float(header.split(' ')[9])
	Tvir, Pvir, Kvir = getSelfSimilarValues(mvir, delta, crit=True, aexp=1.0, omega_m=0.27, omega_l=0.73,	omega_b = 0.0469, hubble=0.7 )
	profile = compute_profile('entropy',snapshot)
	return profile/Kvir

def normPlot(property):
	names=['a04','a10']
	snaps = [0.4,1.0]
	zs = [1.5, 0]
	for ind in xrange(len(snaps)):
		plt.clf()
		scalefactor = snaps[ind]
		colours = cm.rainbow(np.linspace(0,1, len(runs)))
		for run, colour in zip(runs, colours):
			snapshot = closest_snap_run(scalefactor, rundir=run)
			
			virfile = glob.glob(snapshot+'/halo_list_500c_a*')[0]
			with open(virfile,'r') as file:
			    header = file.readlines()[-1]
			Rvir = float(header.split(' ')[1])

			normed_entropy = Prop_norm_delta(property, snapshot)
			
			label = runname[runs==run]
			plt.plot(radii/Rvir, normed_entropy, color=colour, label =label)
		plt.yscale('log')
		plt.xscale('log')
		plt.ylim(1e-3, 10)
		plt.xlim(.03, 13)
		plt.ylabel(r'K/K$_{crit}$')
		plt.xlabel(r'r/R$_{500}$')
		plt.legend(loc='best')
		plt.title('z='+str(zs[ind]))
		# name = str(scalefactor*10)
		plt.savefig(homedir+'/'+names[ind])

m_H = 1.67e-27 #g
def E(z, omega_m=0.27, omega_l=0.73):
	return np.sqrt(omega_m*pow(1+z,3)+omega_l)

def selfSimilarity():
	fig = plt.figure()
	ax = fig.add_subplot(111)
	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111)
	colours = cm.rainbow(np.linspace(0,1, len(runs)))
	for run, colour in zip(runs, colours):
		snaps = sorted(glob.glob(os.getcwd()+run+'/profile_analysis/profiles/*'))
		a_s = np.asarray([float(snap.split('/')[-1]) for snap in snaps])
		snaps = np.ma.masked_where(a_s<0.4, snaps).compressed()
		a_s = a_s[a_s>0.4]
		n_z_r_array = np.empty([len(a_s),len(radii)+1])
		for scalefactor in xrange(len(a_s)):
			print run, scalefactor
			n_z_r_array[scalefactor,0] = 1./(1+a_s[scalefactor])
			snapshot = snaps[scalefactor]
			
			density_profile = compute_profile('density', snapshot)
			n_z_r_array[scalefactor,1:] = density_profile/(1.12*m_H)

		#now we have an array of radial number density profiles at all redshifts between 1.5 and 0
		#first column is redshifts, i.e. x-axis of plots
		ne_today = n_z_r_array[-1]
		n_z_r_array[:,1:]/=ne_today[1:]
		alpha_mean = np.empty(len(radii))
		alpha_std = np.empty(len(radii))
		label = runname[runs==run]
		for rbin in xrange(len(radii)):
			n_r = n_z_r_array[:,rbin+1]
			alpha = np.divide(np.log10(n_r), np.log10(E(n_z_r_array[:,0])))
			alpha_mean[rbin] = np.mean(alpha)
			alpha_std[rbin] = np.std(alpha)
		ax.plot(radii, alpha_mean, color=colour,label = label)

		fig2.clf()
		ax2.plot(radii, alpha_mean, c='k')
		ax2.fill_between(radii, alpha_mean-alpha_std, alpha_mean+alpha_std, color='k',alpha=0.5)
		ax2.set_xlabel('R (kpc)')
		ax2.set_ylabel(r'$/alpha$')
		fig2.savefig(homedir+'/'+label+'_selfsimilarity.png')
	ax.set_ylabel(r'$/alpha$')
	ax.set_xlabel('R (kpc)')
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels)
	fig.savefig(homedir+'/selfsimilarity_comparison.png')

if __name__ == '__main__':
	normEntropyPlot()
	selfSimilarity()