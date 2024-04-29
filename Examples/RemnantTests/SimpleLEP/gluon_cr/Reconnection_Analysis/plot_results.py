import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
#from matplotlib.ticker import PercentFormatter
import statistics as stats

np.random.seed(42)
n_points = 10000
#n_bins   = 100000
n_bins = 10

#filenames = ["stat_mass_beforeCR.txt", "stat_mass_afterCR.txt", "gluon_mass_beforeCR.txt",
#	     "gluon_mass_afterCR.txt", "stat_stringlength.txt", "gluon_stringlength.txt"]

filenames = [ "stat_stringlength.txt", "gluon_stringlength.txt" ]

for name in filenames:
	data = np.loadtxt(name)
#	data = [x for x in data if x != 0]
	fig, axs = plt.subplots(1,1,tight_layout=True)
#	plt.xlim(0,100)
	axs.hist(data, bins=n_bins)
	plt.xlabel("total string-length")
	plt.ylabel("entries")
	plt.xlim(0,10)
#	plt.xscale('log')
	plt.show()
	plt.savefig(name.replace(".txt", ".pdf"))
	
	print(name + " -------------------------")
	print("     Mean: ", stats.mean(data))
	print("   Median: ",stats.median(data) )
	print("  st. dev: ", stats.stdev(data))
	print("max value: ", max(data))
	print("min value: ", min(data))
	print("\n")
