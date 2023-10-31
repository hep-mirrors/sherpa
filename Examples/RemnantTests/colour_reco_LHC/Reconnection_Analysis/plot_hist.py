import matplotlib.pyplot as plt
import numpy as np

f = np.genfromtxt("stringlength_vec.txt", delimiter="\n")

#bins = np.linspace(100, 10e19, 200)
bins = np.linspace(0, 10e12, 100)

plt.hist(f, histtype='bar', bins=bins)
plt.xlabel("String Length [?]")
plt.show()
