import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
#from matplotlib.ticker import PercentFormatter

np.random.seed(42)
n_points = 10000
n_bins   = 50

filenames = ["stat_mass_beforeCR.txt", "stat_mass_afterCR.txt", "gluon_mass_beforeCR.txt",
	     "gluon_mass_afterCR.txt", "stat_stringlength.txt", "gluon_stringlength.txt"]

for name in filenames:
	data = np.loadtxt(name)
	fig, axs = plt.subplots(1,1,tight_layout=True)
	axs.hist(data, bins=n_bins)
	plt.xlabel("total string-length")
	plt.ylabel("entries")
	#plt.show()
	plt.savefig(name.replace(".txt", ".pdf"))
