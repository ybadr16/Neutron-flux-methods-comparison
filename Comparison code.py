import numpy as np
import matplotlib.pyplot as plt
from scipy.special import k0, i0

# Define the parameters
a = 40
So = 20     # neutrons/second
D = 1
sigma_a = 0.1    #Macroscopic absorption crosssection
L = np.sqrt(D / sigma_a)  # Diffusion length
mean_free_path = 1/ (sigma_a)

# Create a grid of points
r = np.linspace(0, a, num = 100000)
#dr = r[1] - r[0]
# Initialize the solutions as arrays equal to r in lengths
phi_analytical = np.zeros(len(r))
phi_montecarlo = np.zeros(len(r))

# Calculate the analytical solution
phi_analytical = So / (2 * np.pi * D) * (k0(r / L) - i0(r / L) * k0(a / L) / i0(a / L))

# bin boundaries and widths for Monte Carlo simulation
n_bins = 50
bin_boundaries = np.linspace(0, a, n_bins+1)
bin_widths = np.diff(bin_boundaries)

# Compute bin volumes, where pi* (r1^2-r2^2) would equal the volume of a cylindircal shell
bin_volumes = np.pi*(bin_boundaries[1:]**2 - bin_boundaries[:-1]**2)

# Initialize neutron flux array for Monte Carlo simulation
flux = np.zeros(n_bins)

# Initialize positions and directions of particles for Monte Carlo simulation
n_particles = 100000
positions = np.zeros(n_particles)
directions = np.random.uniform(0, 2*np.pi, n_particles)

# Perform Monte Carlo simulation
for i in range(n_particles):
    #Hashing out the section below for now because it's not working
    # Compute distance traveled
    #if i == 0:
    #    distance_traveled = positions[i]
    #else:
    #    distance_traveled = positions[i] - positions[i-1]

    # Move particle a small distance
    pos = positions[i] + mean_free_path * np.cos(directions[i])  # getting the position of the neutron in the x direction and mean free path, maybe try some other approaches

    # Check if particle is outside the cylinder
    if pos > a:
        continue

    # Compute bin index based on position
    bin_index = int(pos / (bin_widths[0]))

    # Compute weight based on distance traveled
    #weight = (sigma_a*np.exp(-sigma_a*pos)) #* dr
    #weight = 1.0               #Removing the weight for trial right now current considerations: not using mean free path, using probabilities
    weight = So/(4*np.pi*D)

    # Add contribution to neutron flux
    flux[bin_index] += weight

    # Scatter particle
#directions[i] += np.random.normal(scale=np.sqrt(D/L))

# Normalize neutron flux by bin width and volume
#flux *= So
flux /= bin_widths*bin_volumes*n_particles

# Compute Monte Carlo solution by averaging over adjacent bins
for i in range(len(r)):
    bin_index = int(r[i]/(bin_widths[0]))
    if bin_index+1 >= n_bins:
        bin_index = n_bins-2
    phi_montecarlo[i] = 0.5*(flux[bin_index] + flux[bin_index+1])*1000  # the 1000 is to make the units similar

# Plot the solutions together
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(r, phi_analytical, label='Analytical')
ax.plot(r, phi_montecarlo, label='Monte Carlo')
ax.set_xlabel('Radial position (cm)')
ax.set_ylabel('Neutron flux (cm$^{-2}$ s$^{-1}$)')
# Plot the solutions side by side
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
axs[0].plot(r, phi_analytical)
axs[0].set_xlabel("r")
axs[0].set_ylabel("phi")
axs[0].set_title("Analytical solution")
axs[1].plot(r, phi_montecarlo)
axs[1].set_xlabel('Radial position (cm)')
axs[1].set_ylabel('Neutron flux (cm$^{-2}$ s$^{-1}$)')
axs[1].set_title('Monte Carlo solution')
plt.show()



