"""
Plot the drag coefficient over 4 oscillation cycles.
"""

import pathlib
import numpy
from matplotlib import pyplot
from scipy import signal


# Read the drag force from file.
simu_dir = pathlib.Path(__file__).parents[1]
data_dir = simu_dir / 'output'
filepath = data_dir / 'forces-0.txt'
with open(filepath, 'r') as infile:
    t, fx = numpy.loadtxt(infile, dtype=numpy.float64,
                          unpack=True, usecols=(0, 1))

# Set the parameters of the kinematics.
KC = 5.0  # Keulegan-Carpenter number
D = 1.0  # cylinder diameter
f = 0.2  # frequency of oscillation
w = 2 * numpy.pi * f  # angular frequency
Am = KC * D / (2 * numpy.pi)  # amplitude of oscillation
rho = 1.0  # fluid density
Um = w * Am  # maximum translational velocity of cylinder
V = numpy.pi * D**2 / 4  # volume of cylinder

# Add force due to body acceleration.
ax = w**2 * Am * numpy.sin(w * t)
fx += rho * V * ax
# Get the drag coefficient.
cd = fx / (0.5 * rho * Um**2 * D)

# Compute and print info abount extrema of the drag coefficient.
idx_min = signal.argrelextrema(fx, numpy.less_equal, order=100)[0][1:-1]
t_min = t[idx_min]
print('Non-dimensional time-interval between minima:\n\t{}'
      .format(f * (t_min[1:] - t_min[:-1])))
cd_min = cd[idx_min]
print('Drag coefficient valleys: {}'.format(cd_min))
idx_max = signal.argrelextrema(fx, numpy.greater_equal, order=100)[0][1:]
t_max = t[idx_max]
print('Non-dimensional time-interval between maxima:\n\t{}'
      .format(f * (t_max[1:] - t_max[:-1])))
cd_max = cd[idx_max]
print('Drag coefficient peaks: {}'.format(cd_max))

# Plot the drag coefficient over the 4 cycles.
pyplot.rcParams['font.size'] = 16
pyplot.rcParams['font.family'] = 'serif'
fig, ax = pyplot.subplots(figsize=(8.0, 4.0))
ax.grid()
ax.set_xlabel('$t / T$')
ax.set_ylabel('$C_D$')
ax.plot(f * t, cd)
ax.axis((0.0, 4.0, -6.0, 6.0))
fig.tight_layout()
pyplot.show()

# Save the figure.
fig_dir = simu_dir / 'figures'
fig_dir.mkdir(parents=True, exist_ok=True)
filepath = fig_dir / 'dragCoefficient.png'
fig.savefig(str(filepath), dpi=300)
