from scipy.interpolate import UnivariateSpline # v0.19.1 must be used
import numpy as np
from math import cos, sqrt, atan, floor
import matplotlib
matplotlib.use('TkAgg')
import pylab as plt
import sys
from process_data import *

print('The periods of each spline interpolation are {}, respectively.'.format(period))

fig = plt.figure(1, figsize=(9, 7), dpi=150)
fig.subplots_adjust(hspace=0.76,left=0.11, top=0.94, right=0.74, bottom=0.1)
anchor = (1.4, 1.02)

ax1 = plt.subplot(311)
ax1.plot(x1, y1, '.', c='black', label='Processed\nData Points')
ax1.plot(xs1[:-cutoff_index_1], spl1(xs1[:-cutoff_index_1]), '--', c='black', lw=1, label='Spline\nInterpolation')
plt.ylabel('\\textbf{Torque Felt by\\\\Driven Gear (Nm)}')
plt.xlabel('\\textbf{Difference in Displacements (rad)}')
plt.legend(numpoints=1,bbox_to_anchor=anchor)
ax1.set_title('\\textbf{Spline Interpolation for Trial 1}')

ax2 = plt.subplot(312)
ax2.plot(x2, y2, '*', c='black', label='Processed\nData Points')
ax2.plot(xs2[:-cutoff_index_2], spl2(xs2[:-cutoff_index_2]), '--', c='black', lw=1, label='Spline\nInterpolation')
plt.ylabel('\\textbf{Torque Felt by\\\\Driven Gear (Nm)}')
plt.xlabel('\\textbf{Difference in Displacements (rad)}')
plt.legend(numpoints=1,bbox_to_anchor=anchor)
ax2.set_title('\\textbf{Spline Interpolation for Trial 2}')

ax3 = plt.subplot(313)
plt.plot(x3, y3, 'x', c='black', label='Processed\nData Points')
plt.plot(xs3[:-cutoff_index_3], spl3(xs3[:-cutoff_index_3]), '--', c='black', lw=1, label='Spline\nInterpolation')
plt.ylabel('\\textbf{Torque Felt by\\\\Driven Gear (Nm)}')
plt.xlabel('\\textbf{Difference in Displacements (rad)}')
plt.legend(numpoints=1,bbox_to_anchor=anchor)
ax3.set_title('\\textbf{Spline Interpolation for Trial 3}')

plt.savefig('processed_w_fit.png')
plt.show()
