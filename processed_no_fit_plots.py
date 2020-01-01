import numpy as np
from scipy.optimize import leastsq
import matplotlib
matplotlib.use('TkAgg') # Actually lets the graph be displayed
import pylab as plt
from process_data import *


fig = plt.figure(1, figsize=(9, 5), dpi=150)
#fig.suptitle('\\textbf{Torque Felt by Driven Gear vs. Difference in Displacements}', fontweight='bold')
fig.subplots_adjust(left=0.08, top=0.97, right=0.98, bottom=0.1)
plt.plot(x1, y1, '.', c='black', label='Trial 1')
plt.plot(x2, y2, '*', c='black', label='Trial 2')
plt.plot(x2, y2, 'x', c='black', label='Trial 3')
plt.ylabel('\\textbf{Torque Felt by Driven Gear (Nm)}')
plt.xlabel('\\textbf{Difference in Displacements (rad)}')


plt.legend(numpoints=1, handletextpad=0.1)
plt.show()
