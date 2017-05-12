import yt

sciencedir = '/lustre/scratch/client/fas/nagai/projects/L100_0'
homedir = '/home/fas/nagai/uc24'
os.chdir(sciencedir)
runs = ['/CSF','/AGN/CCA/recentering/fiducial','/AGN/CCA/largethermal/fiducial', '/AGN/CCA/largethermal/rfb10', '/AGN/CCA/largethermal/rfb50','/AGN/CCA/largethermal/slope1','/AGN/CCA/largethermal/slope3']
runname = ['No AGN', r'$r_{FB}$=3, $\alpha$=2', r'$r_{FB}$=25, $\alpha$=2', r'$r_{FB}$=10, $\alpha$=2', r'$r_{FB}$=50, $\alpha$=2', r'$r_{FB}$=25, $\alpha$=1', r'$r_{FB}$=25, $\alpha$=3']

def run_snaps(run):
    os.chdir(sciencedir+run)
    snaps = sorted(glob.glob(os.getcwd()+"/DAT/*.art"))
    a_plus = [snap.split('a0')[-1] for snap in snaps[2:-1]] #to leave out a=1
    a_sp = [aplus.split('.')[-2] for aplus in a_plus]
    return ['0.'+asp for asp in a_sp]

def snapsdone(property = "entropy"):
    return len(glob.glob(os.getcwd()+'/yt/'+property+'/*'))

def find_center(ds):
	v, c = ds.find_max("density")
	return c

def formatPlot(plot, center):
    plot.set_center((center[0],center[2]))
    plot.set_width((2,'Mpc'))
    plot.set_axes_unit('Mpc')
    plot.save()

def make_images(run):
    os.chdir(sciencedir+run)
    os.mkdir(yt)
	for a in run_snaps(run): #i think this should do the chdir as well. test.
		name = str(a)
		if len(name) < 6:
			for i in xrange(6-len(name)):
				name += '0'
		ds = yt.load("DAT/L100_a%s.art" %  name)
        center = find_center(ds)
		col_density = yt.ProjectionPlot(ds, "y", 'density')
        formatPlot(col_density, center)
        col_density.save()
        for property in ['temperature','pressure','entropy']:
            plot = yt.SlicePlot(ds, 'y', property)
            formatPlot(plot, center)

if __name__=='__main__':
	for run in runs:
        make_images(run)