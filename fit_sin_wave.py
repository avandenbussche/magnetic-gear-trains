import numpy as np
from scipy.optimize import leastsq
import matplotlib
matplotlib.use('TkAgg')
import pylab as plt
from math import sqrt, atan, cos
from process_data import *

guess_mean = np.mean(y1)/2
guess_std = 3*np.std(y1)/(2**0.5)
guess_phase = 0
guess_stretch = 0.3


data_first_guess = guess_std*np.sin(np.sin(guess_stretch**-1 * (x1))) + guess_mean

optimize_func = lambda x: x[0]*np.sin(np.sin(x[1]**-1 *(x1))) - y1
est_std, est_stretch, est_mean = leastsq(optimize_func, [guess_std, guess_stretch, guess_mean])[0]

fig = plt.figure(1, figsize=(9, 5), dpi=150)
fig.suptitle('\\textbf{Torque Felt by Driven Gear vs. Difference in Displacements}', fontweight='bold')
fig.subplots_adjust(left=0.11, top=0.9, right=0.98, bottom=0.1)

plt.plot(x1, y1, '.', label='Processed Data Points', c='black')
plt.plot(x1, est_std*np.sin(est_stretch**-1 *(x1)+est_mean), '--', label='Fitted Sine Wave', c='black')

plt.ylabel('\\textbf{Torque Felt by\\\\Driven Gear (Nm)}')
plt.xlabel('\\textbf{Difference in Displacements (rad)}')
plt.xlim(0, np.pi/2)

plt.legend(numpoints=1)
plt.show()
