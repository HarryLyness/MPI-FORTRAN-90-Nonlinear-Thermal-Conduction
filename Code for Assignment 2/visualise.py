"""Visualise solution in file 'solution.txt'

It is assumed that solution.txt contains the following data:

  * a line containing the problem size m [single integer]
  * a blank line
  * a line containing (m+1) x (m+1) real numbers encoding the solution U

The solution in assumed to be stored in column-major format, i.e.

  U = (u_{0,0},u_{1,0},...,u_{m,0},u_{0,1},...,u_{m,1,},...,u_{m+1,m+1})

This script produces the plot in the file 'solution.pdf'
"""

import numpy as np
from matplotlib import pyplot as plt

print ("Reading data from file solution.txt")

# Read data
with open("solution.txt","r") as f:
    content = f.readlines()
    # check that file contains exactly three lines
    assert len(content) == 3, "Expected file to contain exactly three lines"
    header, _, raw_data = content
    m = int(header)
    print (f"problem size m = {m:4d}")
    # check that the final line contains exactly (m+1)^2 floating point numbers
    assert len(raw_data.split()) == (m+1)**2, f"Expected third line to contain exactly (m+1)^2 = {(m+1)**2:8d} floating point numbers"
    data = np.asarray([float(x) for x in raw_data.split()]).reshape((m+1,m+1)).T

plt.imshow(data,aspect=1.0,cmap="jet",origin="lower",extent=[0,1,0,1])
plt.colorbar()
ax = plt.gca()
ax.set_xlabel(r"$x_1$")
ax.set_ylabel(r"$x_2$")

print ("Success! solution.pdf was created")
plt.savefig("solution.pdf",bbox_inches="tight")
