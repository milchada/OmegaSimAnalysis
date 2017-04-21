#plot all quantities normalized to values as R_500c
#for comparison with McDonald et al
from massfrac import *
from normalize import *

mH = 1.67e-24 #g
"""
CONFIRM PRESSURE KEY BELOW FROM gas_colnames.txt on omega
"""
keys = ['gas_temperature(mass-weighted) [K','gas_pressure(mass-weighted) [','gas_entropy(mass-weighted) [keV cm ^ 2', 'gas_density(emission-weighted) [g cm^{-3}', 'gas X-ray emission [erg /s']
proplist = ['temperature', 'pressure', 'entropy', 'density']
propnames = ['T', 'P', 'K']

def compute_profile(property, snap):
	gas = read_run(snap)

	with open(homedir+'/gas_colnames.txt','r') as file:
	    colnames=file.readlines()[1]
	gas_colnames = colnames.split('] ')
	
	propkey = [key for key in keys if property in key][0]
	# print gas_colnames.index(propkey), gas.keys()
	return gas[:,gas_colnames.index(propkey)]

def Prop_norm_delta(property='entropy', snapshot, delta=200):
	massfile = glob.glob(snapshot+'/halo_profile_ma*')[0]
	with open(massfile,'r') as file:
		header = file.readlines()[10]
	mvir = float(header.split(' ')[9])
	prop_virs = getSelfSimilarValues(mvir, delta, crit=True, aexp=1.0, omega_m=0.27, omega_l=0.73,	omega_b = 0.0469, hubble=0.7 ) #T, P, K
	propvir = prop_virs[proplist.index(property)] 
	profile = compute_profile(property,snapshot)
	return profile/propvir

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

def delta_r(snapshot, fineness=10, crit=True, omega_m=0.27, omega_l=0.73,    omega_b = 0.0469, hubble=0.7 ) :
    aexp = float(snapshot.split('/')[-1])
    Ez = E(1./aexp - 1)
    rho_crit_a = rho_crit_0*(Ez**2)
    mass = read_run(snapshot, prop='mass')
    gas_index = m.profile_columns['gas_M_cum']
    dm_index = m.profile_columns['dark_M_cum']
    star_index = m.profile_columns['star_M_cum']
    total_mass = mass[:,gas_index]+mass[:,dm_index]+mass[:,star_index]
    fine_radii = np.arange(radii.min, radii.max, len(radii)*fineness)
    fine_mass = np.interp(fine_radii, radii, total_mass)
    volume = 4*np.pi*np.power(fine_radii,3)/3
    return fine_radii, np.divide(fine_mass, volume*rho_crit_a) #unitless
    
def Rdelta(delta_r, delta, fine_radii):
	return fine_radii[abs(delta_r - delta)==min(abs(delta_r - delta))]

def profile_gas(snapshot, fineness = 10, delta=False, omega_m=0.27, omega_l=0.73, omega_b= 0.0469):
	aexp = float(snapshot.split('/')[-1])
	Ez = E(1./aexp - 1)
	rho_crit_gas = rho_crit_0*(Ez**2)*omega_b
	mass = read_run(snapshot, prop='mass')
	gas_index = m.profile_columns['gas_M_cum']
	fine_radii = np.arange(radii.min, radii.max, len(radii)*fineness)
	gas_profile = np.interp(fine_radii, radii, mass[:,gas_index])
	if delta==True:
		delta_gas = gas_profile/rho_crit_gas
	else:
		delta_gas = gas_profile/(Ez**2*mH)
	return delta_gas


def normPlot(property, snapshot):
	scalefactor = float(snapshot.split('/')[-1])
	z = 1./aexp - 1
	plt.clf()
	colours = cm.rainbow(np.linspace(0,1, len(runs)))
	for run, colour in zip(runs, colours):
		snapshot = closest_snap_run(scalefactor, rundir=run)
		fine_radii, overdensities = delta_r(snapshot)
		Rvir = Rdelta(overdensities, 500, fine_radii)
		if property == 'density':
			normed_profile = profile_gas(snapshot) #following Nagai07b. If we want overdensity, set delta=True
		else:
			normed_profile = Prop_norm_delta(property, snapshot)		
		label = runname[runs==run]
		plt.plot(radii/Rvir, normed_profile, color=colour, label =label)
	plt.yscale('log')
	plt.xscale('log')
	plt.xlim(.03, 3)
	propname = propnames[proplist==property]
	if property == 'density':
		plt.ylabel(r'$n_{gas}(r)E^{-2}(z)[cm^{-3}]$')
	else:
		ylabel = propname+'/'+propname
		plt.ylabel(ylabel+r'$_{500}$')
	plt.xlabel(r'r/R$_{500}$')
	plt.legend(loc='best')
	plt.title('z='+str(z))
	plt.savefig(homedir+'/'+property+str(scalefactor)+'_normed.png')

if __name__ == '__main__':
	normPlot()
