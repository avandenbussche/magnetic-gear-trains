# Runge-Kutta Method for Predicting Slipping in Accelerating Magnetic Gear Trains
# Script written by Adam Vandenbussche on October 8, 2017
# Based on the method outlined in Elementary Differential Equations with Applications 2nd ed. (Derrick & Grossman, 1981)
#
# Run like this:
# python runge-kutta.py a+-a_unc t h I_0+-I_0_unc I_1+-I_1_unc I_2+-I_2_unc --ignore-tolerance
# where: a+-a_unc -> Driver acceleration and absolute uncertainty
#        t -> Total time interval to compute
#        h -> Runge-Kutta increment
#        I_0+-I_0_unc -> Coefficient of velocity^0 of driven gear (moment of inertia) and absolute uncertainty
#        I_1+-I_1_unc -> Coefficient of velocity^1 of driven gear (moment of inertia) and absolute uncertainty
#        I_2+-I_2_unc -> Coefficient of velocity^2 of driven gear (moment of inertia) and absolute uncertainty
#        --ignore-tolerance -> Ignore maximum displacement per interval warning
#

# Import dependencies
from math import pi, ceil, log
import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
import pylab as plt
import sys
from process_data import *

def main(args):

    #
    # PROCESS PROGRAM PARAMETERS
    #

    h = float(args[3]) # Step size
    total_time = float(args[2]) # Seconds
    current_time = 0

    # Driver acceleration
    driver_a, driver_a_unc = args[1].split('+-')
    driver_a = float(driver_a)
    driver_a_unc = float(driver_a_unc)

    total_provided_unc = driver_a_unc / driver_a

    driver_v = 0

    # Moment of inertia quadratic coefficients
    moi_0, moi_0_unc = args[4].split('+-')
    moi_0 = float(moi_0)
    moi_0_unc = float(moi_0_unc)

    if moi_0 != 0:
        total_provided_unc += moi_0_unc / moi_0

    moi_1, moi_1_unc = args[5].split('+-')
    moi_1 = float(moi_1)
    moi_1_unc = float(moi_1_unc)

    if moi_1 != 0:
        total_provided_unc += moi_1_unc / moi_1

    moi_2, moi_2_unc = args[6].split('+-')
    moi_2 = float(moi_2)
    moi_2_unc = float(moi_2_unc)

    if moi_2 != 0:
        total_provided_unc += moi_2_unc / moi_2


    #
    # MAIN RUNGE-KUTTA ALGORITHM
    #

    iterations = ceil(total_time/h)
    successful_iterations = 0
    already_alerted = False

    time_xs = np.empty([1, iterations])
    displacement_ys = np.empty([1, iterations])
    velocity_ys = np.empty([1, iterations]) ######
    acceleration_ys = np.empty([1, iterations])
    acceleration_ys_unc = np.empty([1, iterations])
    velocity_ys_unc = np.empty([1, iterations])
    displacement_ys_unc = np.empty([1, iterations])

    def driven_acceleration(time, driven_displacement, driven_velocity):
        delta_theta = 0.5 * driver_a * time ** 2 - driven_displacement
        values = [ spl1( delta_theta % period[0]), spl2( delta_theta % period[1]), spl3( delta_theta % period[2]) ]

        average = sum(values) / 3
        range_unc = (max(values) - min(values)) / 2

        moi = moi_0 + moi_1 * driven_velocity + moi_2 * driven_velocity ** 2

        return (average / moi, range_unc / moi, [values[0]/moi, values[1]/moi, values[2]/moi])

    # Initial conditions
    d_current = 0
    v_current = 0

    for i in range(iterations):
        #print('Performing calculation {0:,d} of {1:,d}!'.format(i, iterations))
        driver_v += driver_a * h
        driver_d = 0.5 * driver_a * h ** 2
        if (driver_v * h + driver_d > pi/24):
            if not already_alerted:
                print('Surpassed maximum displacement per iteration threshold after {0} seconds! Stopping simulation.'.format(current_time))
                already_alerted = True
            if '--ignore-tolerance' not in args:
                break

        time_xs[0][i] = current_time
        displacement_ys[0][i] = d_current
        velocity_ys[0][i] = v_current

        dd1 = h * v_current
        dv1, dv1_a_unc, dv1_a_individual = driven_acceleration(current_time, d_current, v_current)
        dv1 *= h

        dd2 = h * (v_current + dv1 / 2)
        dv2, dv2_a_unc, dv2_a_individual = driven_acceleration(current_time + h / 2, d_current + dd1 / 2, v_current + dv1 / 2)
        dv2 *= h

        dd3 = h * (v_current + dv2 / 2)
        dv3, dv3_a_unc, dv3_a_individual = driven_acceleration(current_time + h / 2, d_current + dd2 / 2, v_current + dv2 / 2)
        dv3 *= h

        dd4 = h * (v_current + dv3)
        dv4, dv4_a_unc, dv4_a_individual = driven_acceleration(current_time + h, d_current + dd3, v_current + dv3)
        dv4 *= h

        dd = (dd1 + 2 * dd2 + 2 * dd3 + dd4) / 6
        dv = (dv1 + 2 * dv2 + 2 * dv3 + dv4) / 6
        d_current += dd
        v_current += dv

        acceleration_ys[0][i], acceleration_ys_unc[0][i] = driven_acceleration(current_time + h, d_current, v_current)[:2]

        velocity_ys_unc[0][i] = h * acceleration_ys_unc[0][i]
        displacement_ys_unc[0][i] = h * velocity_ys_unc[0][i] + 0.5 * acceleration_ys_unc[0][i] * (h ** 2)

        current_time += h
        successful_iterations += 1


    difference_ys = np.empty([1, successful_iterations])
    difference_ys_unc = np.empty([1, successful_iterations])

    # Analyze slipping
    done_analyzing_slipping = False
    slipping_time = 0
    earliest_possible_slipping_time = 0
    latest_possible_slipping_time = 0
    for i in range(successful_iterations):
        difference_ys[0][i] = 0.5 * driver_a * (h * i) ** 2 - displacement_ys[0][i]
        difference_ys_unc[0][i] = driver_a_unc / driver_a + displacement_ys_unc[0][i]
        if not done_analyzing_slipping and earliest_possible_slipping_time == 0 and difference_ys[0][i] * (1 + total_provided_unc) + difference_ys_unc[0][i] > average_period - period_unc:
            earliest_possible_slipping_time = h * i
        if not done_analyzing_slipping and slipping_time == 0 and difference_ys[0][i] > average_period:
            slipping_time = h * i
        if not done_analyzing_slipping and difference_ys[0][i] * (1 - total_provided_unc) - difference_ys_unc[0][i] > average_period + period_unc:
            done_analyzing_slipping = True
            latest_possible_slipping_time = h * i
    if done_analyzing_slipping:
        print('Slipping occurs after {0} s, although could occur between {1} and {2} s ({0}+{3}-{4}).'.format(slipping_time, earliest_possible_slipping_time, latest_possible_slipping_time, slipping_time-earliest_possible_slipping_time, latest_possible_slipping_time-slipping_time))
    else:
        print('No slipping occurs!')

    # Determine peak driven gear acceleration
    max_driven_a = np.amax(acceleration_ys)
    max_driven_a_unc = max_driven_a * total_provided_unc + acceleration_ys_unc[0][np.argmax(acceleration_ys)]
    print('Peak Driven Acceleration: {} ± {}'.format( max_driven_a, max_driven_a_unc))



    #
    # GENERATE PLOTS
    #

    # Set up window
    fig = plt.figure(1, figsize=(9, 6), dpi=150)
    fig.subplots_adjust(hspace=0.48, wspace=0.24, left=0.1, top=0.94, right=0.98, bottom=0.18)
    size = 0.2 ** 2 # Size of points in pt^2
    x_time = np.linspace(0, h * successful_iterations, successful_iterations)
    opacity = 0.2

    # Acceleration-Time
    ax1 = plt.subplot(221)
    # Driver uncertainty
    unc = ax1.fill_between(time_xs[0][:successful_iterations - 2], np.linspace(driver_a, driver_a, successful_iterations - 2) + np.linspace(driver_a_unc, driver_a_unc, successful_iterations - 2), np.linspace(driver_a, driver_a, successful_iterations - 2) - np.linspace(driver_a_unc, driver_a_unc, successful_iterations - 2), facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)
    # Driven uncertainty
    ax1.fill_between(time_xs[0][:successful_iterations - 2], acceleration_ys[0][:successful_iterations - 2] * (1 + total_provided_unc) + acceleration_ys_unc[0][:successful_iterations - 2], acceleration_ys[0][:successful_iterations - 2] * (1 - total_provided_unc) - acceleration_ys_unc[0][:successful_iterations - 2], facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)
    # Driver
    ax1.plot(x_time, np.linspace(driver_a, driver_a, successful_iterations), '--', c='black', label='Driver Acceleration')
    # Driven
    ax1.scatter(time_xs[:,:successful_iterations - 2], acceleration_ys[:,:successful_iterations - 2], s=0.01, marker='.', color='black', label='Driven Acceleration', alpha=0.5+log(h, 10)*0.04)

    plt.ylabel('\\textbf{Acceleration (rad s\\textsuperscript{-2})}')
    plt.xlabel('\\textbf{Time Elapsed (s)}')
    plt.xlim(0, (successful_iterations - 2) * h)
    ax1.set_title('\\textbf{Acceleration vs. Time}')


    # Velocity-Time
    ax2 = plt.subplot(222)
    # Driver uncertainty
    ax2.fill_between(x_time, (np.linspace(driver_a, driver_a, successful_iterations) * x_time) * (1 + np.linspace(driver_a_unc/driver_a, driver_a_unc/driver_a, successful_iterations)), (np.linspace(driver_a, driver_a, successful_iterations) * x_time) * (1 - np.linspace(driver_a_unc/driver_a, driver_a_unc/driver_a, successful_iterations)), facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)
    # Driven uncertainty
    ax2.fill_between(time_xs[0][:successful_iterations], velocity_ys[0][:successful_iterations] * (1 + total_provided_unc) + velocity_ys_unc[0][:successful_iterations], velocity_ys[0][:successful_iterations] * (1 - total_provided_unc) - velocity_ys_unc[0][:successful_iterations], facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)

    driver = ax2.plot(x_time, driver_a*x_time, '--', c='black', label='Driver Velocity')
    driven = ax2.scatter(time_xs[:,:successful_iterations], velocity_ys[:,:successful_iterations], s=0.005, marker='.', color='black', label='Driven Velocity')

    plt.ylabel('\\textbf{Velocity (rad s\\textsuperscript{-1})}')
    plt.xlabel('\\textbf{Time Elapsed (s)}')
    plt.xlim(0, successful_iterations * h)
    ax2.set_title('\\textbf{Velocity vs. Time}')


    # Difference in Displacements-Time
    ax3 = plt.subplot(223)

    ax3.fill_between(time_xs[0][:successful_iterations], difference_ys[0][:successful_iterations] * (1 + total_provided_unc) + difference_ys_unc[0][:successful_iterations], difference_ys[0][:successful_iterations] * (1 - total_provided_unc) - difference_ys_unc[0][:successful_iterations], facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)

    ax3.scatter(time_xs[:,:successful_iterations], difference_ys[0], s=0.01, marker='.', color='black', label='Difference in Displacements')
    plt.ylabel('\\textbf{Difference in \\\\ Displacements (rad)}')
    plt.xlabel('\\textbf{Time Elapsed (s)}')
    plt.xlim(0, successful_iterations * h)
    ax3.set_title('\\textbf{Difference in Displacements vs. Time}')


    # Displacement-Time
    ax4 = plt.subplot(224)
    # Driver uncertainty
    ax4.fill_between(x_time, (0.5 * np.linspace(driver_a, driver_a, successful_iterations) * x_time ** 2) * (1 + np.linspace(2*driver_a_unc/driver_a, 2*driver_a_unc/driver_a, successful_iterations)), (0.5 * np.linspace(driver_a, driver_a, successful_iterations) * x_time ** 2) * (1 - np.linspace(2*driver_a_unc/driver_a, 2*driver_a_unc/driver_a, successful_iterations)), facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)
    # Driven uncertainty
    ax4.fill_between(time_xs[0][:successful_iterations], displacement_ys[0][:successful_iterations] * (1 + total_provided_unc) + displacement_ys_unc[0][:successful_iterations], displacement_ys[0][:successful_iterations] * (1 - total_provided_unc) - displacement_ys_unc[0][:successful_iterations], facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)

    ax4.plot(x_time, 0.5*driver_a*x_time**2, '--', c='black', label='Driver Displacement')
    ax4.scatter(time_xs[:,:successful_iterations], displacement_ys[:,:successful_iterations], s=size, marker='.', c='black', label='Driven Displacement')
    plt.ylabel('\\textbf{Displacement (rad)}')
    plt.xlabel('\\textbf{Time Elapsed (s)}')
    plt.xlim(0, successful_iterations * h)
    ax4.set_title('\\textbf{Displacement vs. Time}')


    # Custom legend
    labels = ['Driver Gear', 'Driven Gear']
    markers = ['_', '.']
    colors = ['black', 'black']
    patches = [ plt.plot([],[], marker=markers[i], ms=10, ls='', mec=None, color=colors[i], label='{:s}'.format(labels[i]) )[0] for i in range(len(labels)) ]
    labels.append('Uncertainty')
    patches.append(unc)

    fig.legend(patches, labels, numpoints=1, loc='lower center', ncol=3)

    # Generate separate velocity-time and displacemen-time graphs for driver and driven to emphasize that they are the same

    # Velocity-Time
    fig2 = plt.figure(2, figsize=(9, 2.5), dpi=150)
    fig2.subplots_adjust(hspace=0.48, wspace=0.24, left=0.1, top=0.9, right=0.98, bottom=0.17)
    axA = plt.subplot(121)
    # Driver uncertainty
    axA.fill_between(x_time, (np.linspace(driver_a, driver_a, successful_iterations) * x_time) * (1 + np.linspace(driver_a_unc/driver_a, driver_a_unc/driver_a, successful_iterations)), (np.linspace(driver_a, driver_a, successful_iterations) * x_time) * (1 - np.linspace(driver_a_unc/driver_a, driver_a_unc/driver_a, successful_iterations)), facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)

    driver = axA.plot(x_time, driver_a*x_time, '--', c='black', label='Driver Velocity')

    plt.ylabel('\\textbf{Velocity (rad s\\textsuperscript{-1})}')
    plt.xlabel('\\textbf{Time Elapsed (s)}')
    plt.xlim(0, successful_iterations * h)
    axA.set_title('\\textbf{Driver Gear Velocity vs. Time}')

    axB = plt.subplot(122)

    # Driven uncertainty
    axB.fill_between(time_xs[0][:successful_iterations], velocity_ys[0][:successful_iterations] * (1 + total_provided_unc) + velocity_ys_unc[0][:successful_iterations], velocity_ys[0][:successful_iterations] * (1 - total_provided_unc) - velocity_ys_unc[0][:successful_iterations], facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)

    driven = axB.scatter(time_xs[:,:successful_iterations], velocity_ys[:,:successful_iterations], s=0.005, marker='.', color='black', label='Driven Velocity')

    plt.ylabel('\\textbf{Velocity (rad s\\textsuperscript{-1})}')
    plt.xlabel('\\textbf{Time Elapsed (s)}')
    plt.xlim(0, successful_iterations * h)
    axB.set_title('\\textbf{Driven Gear Velocity vs. Time}')




    fig3 = plt.figure(3, figsize=(9, 2.5), dpi=150)
    fig3.subplots_adjust(hspace=0.48, wspace=0.24, left=0.1, top=0.9, right=0.98, bottom=0.17)

    # Displacement-Time
    axC = plt.subplot(121)
    # Driver uncertainty
    axC.fill_between(x_time, (0.5 * np.linspace(driver_a, driver_a, successful_iterations) * x_time ** 2) * (1 + np.linspace(2*driver_a_unc/driver_a, 2*driver_a_unc/driver_a, successful_iterations)), (0.5 * np.linspace(driver_a, driver_a, successful_iterations) * x_time ** 2) * (1 - np.linspace(2*driver_a_unc/driver_a, 2*driver_a_unc/driver_a, successful_iterations)), facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)

    axC.plot(x_time, 0.5*driver_a*x_time**2, '--', c='black', label='Driver Displacement')
    plt.ylabel('\\textbf{Displacement (rad)}')
    plt.xlabel('\\textbf{Time Elapsed (s)}')
    plt.xlim(0, successful_iterations * h)
    axC.set_title('\\textbf{Driver Gear Displacement vs. Time}')

    # Displacement-Time
    axD = plt.subplot(122)
    # Driven uncertainty
    axD.fill_between(time_xs[0][:successful_iterations], displacement_ys[0][:successful_iterations] * (1 + total_provided_unc) + displacement_ys_unc[0][:successful_iterations], displacement_ys[0][:successful_iterations] * (1 - total_provided_unc) - displacement_ys_unc[0][:successful_iterations], facecolor='black', edgecolor='none', alpha=opacity, interpolate=True)

    axD.scatter(time_xs[:,:successful_iterations], displacement_ys[:,:successful_iterations], s=size, marker='.', c='black', label='Driven Displacement')
    plt.ylabel('\\textbf{Displacement (rad)}')
    plt.xlabel('\\textbf{Time Elapsed (s)}')
    plt.xlim(0, successful_iterations * h)
    axD.set_title('\\textbf{Driven Gear Displacement vs. Time}')

    plt.show()

if __name__=='__main__':
    sys.exit(main(sys.argv))
