import numpy as np
import matplotlib.pyplot as plt

# Load the data (assumes whitespace-separated columns)
data = np.loadtxt("pn_total_low.dat", usecols=(0, 1))

# If the file has two columns
x = data[:, 0]
y = data[:, 1]

plt.plot(x, y)
plt.xlabel("X")
plt.ylabel("Y")
plt.title(".dat Plot")
plt.show()

