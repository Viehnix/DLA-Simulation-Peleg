#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import random
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Params
#-------------------------------------------------------------------------------
width, height = 500, 500            # Roomsize
num_particles = 700                 # Number of particles
mean_radius = 5                     # Mean value of particle radius
std_dev_radius = 2                  # Standard deviation of radius
mean_speed = 8                      # Mean value of particle speed
std_dev_speed = 2                   # Standard deviation of set speed
center = (width // 2, height // 2)  # Starting point of the first particle
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Subroutines
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def create_particle(speed, radius):
#-------------------------------------------------------------------------------
# @brief   Creates a moving particle with a random direction, initial position,
#          and a random radius.
#
# @param   speed:
#            (float) The magnitude of the particle's velocity.
#          radius:
#            (float) The radius of the particle.
#
# @return  dict:
#            x  => Initial position on the X-axis (integer).
#            y  => Initial position on the Y-axis (integer).
#            dx => Velocity component along the X-axis (float).
#            dy => Velocity component along the Y-axis (float).
#            radius => Radius of the particle (float).
#-------------------------------------------------------------------------------
    angle = random.uniform(0, 2 * np.pi)  # zufällige Richtung
    dx = speed * np.cos(angle)
    dy = speed * np.sin(angle)
    return {'x': random.randint(0, width), 'y': random.randint(0, height), 'dx': dx, 'dy': dy, 'radius': radius}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def distance(p1, p2):
#-------------------------------------------------------------------------------
# @brief   Calculates the distance between two particles considering their radii.
#
# @param   p1:
#            (dict) The first particle with 'x', 'y' coordinates and 'radius'.
#          p2:
#            (dict) The second particle with 'x', 'y' coordinates and 'radius'.
#
# @return  float:
#            Distance between the two particles.
#-------------------------------------------------------------------------------
    return np.sqrt((p1['x'] - p2['x'])**2 + (p1['y'] - p2['y'])**2)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def move_particle_to_distance(particle, center_particle):
#-------------------------------------------------------------------------------
# @brief   Moves a particle to a fixed distance considering the radii from another particle's center.
#
# @param   particle:
#            (dict) The moving particle with 'x', 'y' coordinates and 'radius'.
#          center_particle:
#            (dict) The fixed particle with 'x', 'y' coordinates and 'radius' 
#            to which the distance is measured.
#
# @return  None:
#            Updates position of the moving particle in place to maintain
#            the desired distance from the center_particle.
#-------------------------------------------------------------------------------
    dist = distance(particle, center_particle)
    min_dist = particle['radius'] + center_particle['radius']
    if dist < min_dist:
        
        dx = particle['x'] - center_particle['x']          # Calculate X Direction of vector
        dy = particle['y'] - center_particle['y']          # Calculate Y Direction of vector
        
        length = np.sqrt(dx**2 + dy**2)                    # Normalize vector to set position exactly 2 * radius away
        if length < 1e-10:                                 # Avoid division by zero (if near 0, sets to 10^-10)
            length = 1e-10  
        
        scale = min_dist / length
        
        particle['x'] = center_particle['x'] + dx * scale  # Update particle's X position
        particle['y'] = center_particle['y'] + dy * scale  # Update particle's Y position
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def check_collision_with_centers(particle):
#-------------------------------------------------------------------------------
# @brief   Checks for a collision between a moving particle and all fixed centers.
#
# @param   $particle:
#            The moving particle to check for collision.
#
# @return  True, center_particle  => A collision was detected, and the colliding center is returned.
#          False, None            => No collision detected with any center.
#-------------------------------------------------------------------------------
    for center_particle in fixed_centers:
        dist = distance(particle, center_particle)  
        min_dist = particle['radius'] + center_particle['radius']
        if dist < min_dist:                         # If the distance is less than the combined radii, a collision occurs
            return True, center_particle            # Return True and colliding center
    return False, None                              # No collision detected
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def move_particles():
#-------------------------------------------------------------------------------
# @brief   Moves particles and checks for collisions with fixed centers or walls.
#          If collision with a fixed center is detected, particle is fixed.
#          Otherwise, position is updated on particle's velocity and boundary reflections.
#
# @return  list:
#            Returns updated list of moving particles after considering all collisions and boundary reflections.
#-------------------------------------------------------------------------------
    new_moving_particles = []                       # List for the new positions of the particles
    for p in moving_particles:
        future_x = p['x'] + p['dx']
        future_y = p['y'] + p['dy']

        # Check if particle collides with any fixed center
        collision, collided_center = check_collision_with_centers(p)
        
        if collision:                               # Set particle to the valid distance from the center and add particle to list of fixed centers
            move_particle_to_distance(p, collided_center)
            fixed_centers.append({'x': p['x'], 'y': p['y'], 'radius': p['radius']})
            continue

        # If no collision, update the position of the particle
        p['x'] = future_x
        p['y'] = future_y
        
        # Boundary handling: Reflection when the particle hits the walls
        if p['x'] <= p['radius']:
            p['x'] = p['radius']
            p['dx'] = -p['dx']
        elif p['x'] >= width - p['radius']:
            p['x'] = width - p['radius']
            p['dx'] = -p['dx']
        if p['y'] <= p['radius']:
            p['y'] = p['radius']
            p['dy'] = -p['dy']
        elif p['y'] >= height - p['radius']:
            p['y'] = height - p['radius']
            p['dy'] = -p['dy']

        new_moving_particles.append(p)
    
    return new_moving_particles
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def calculate_fractal_dimension(positions, grid_sizes=[20, 10, 5, 2]):
#-------------------------------------------------------------------------------
# @brief   Calculates the fractal dimension using the Box-Counting method.
#          The method divides the space into smaller boxes of different sizes 
#          and counts how many boxes contain at least one particle. The fractal 
#          dimension is estimated by fitting a power law to the number of occupied 
#          boxes as a function of the box size.
#
# @param   positions:
#            (list of tuples) A list of particle positions as (x, y) coordinate pairs.
#          grid_sizes:
#            (list of int) The different grid divisions to be used for the Box-Counting method.
#            The actual box size is calculated as width // grid_size.
#
# @return  float:
#            The estimated fractal dimension based on the Box-Counting method.
#
# @source   https://fractalfoundation.org/OFC/OFC-10-5.html
#           https://www.sciencedirect.com/science/article/pii/S2590123020300128
#-------------------------------------------------------------------------------
    box_counts = []                     # List to store the number of occupied boxes
    box_sizes = []                      # List to store the actual box sizes

    for grid_size in grid_sizes:
        box_size = width // grid_size   # Calculate the box size for the current grid
        boxes = set()                   # Set to store occupied boxes

        for p in positions:
            box_x = p[0] // box_size
            box_y = p[1] // box_size
            boxes.add((box_x, box_y))   # Add the occupied box to the set

        box_counts.append(len(boxes))   # Store the number of occupied boxes
        box_sizes.append(box_size)      # Store the current box size

    # Calculate logarithms for log-log regression
    log_counts = np.log(box_counts)
    log_sizes = np.log(box_sizes)

    # Perform linear regression (slope of the line is the fractal dimension)
    fit = np.polyfit(log_sizes, log_counts, 1)
    fractal_dimension = -fit[0]         # The slope of the line gives the fractal dimension

    return fractal_dimension
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def plot_particles():
#-------------------------------------------------------------------------------
# @brief   Visualizes moving and fixed particles, and calculates and displays 
#          the fractal dimension of the fixed particles.
#
# @return  None
#-------------------------------------------------------------------------------
    plt.clf()                                        # Clear previous plot
    plt.xlim(0, width)                               # Set x-axis limits
    plt.ylim(0, height)                              # Set y-axis limits
    plt.gca().set_aspect('equal', adjustable='box')  # Equal aspect ratio

    # Plot moving particles as blue circles
    for p in moving_particles:
        particle = plt.Circle((p['x'], p['y']), p['radius'], color='blue')
        plt.gca().add_artist(particle)

    # Plot fixed center particles as green circles
    for center_particle in fixed_centers:
        particle = plt.Circle((center_particle['x'], center_particle['y']), center_particle['radius'], color='green')
        plt.gca().add_artist(particle)

    # Calculate the fractal dimension of the fixed particles
    positions = [(p['x'], p['y']) for p in fixed_centers]
    fractal_dimension = calculate_fractal_dimension(positions)

    # Get particle counts
    total_particles = len(fixed_centers)
    moving_particles_count = len(moving_particles)

    # Display fractal dimension and particle counts
    plt.text(10, -height // 10, f"Fractal Dimension: {fractal_dimension:.3f}", fontsize=12, color='black', va='top')
    plt.text(10, height + height // 10, f"Fixed Particles: {total_particles - 1}", fontsize=12, color='black', va='top')
    plt.text(10, height + height // 20, f"Moving Particles: {moving_particles_count}", fontsize=12, color='black', va='top')
    plt.text(width - width // 4, height + height // 20, f"Steps: {steps}", fontsize=12, color='black', va='top')

    plt.draw()                                       # Update the plot
    plt.pause(0.05)                                  # Pause for animation effect
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------
moving_particles = []                               # List to store moving particles
fixed_centers = []                                  # List to store fixed center particles
steps = 0

# Add the first particle at the center with zero velocity
moving_particles.append({'x': center[0], 'y': center[1], 'dx': 0, 'dy': 0, 'radius': mean_radius})
fixed_centers.append({'x': center[0], 'y': center[1], 'radius': mean_radius})

# Generate normally distributed speeds and radii for the particles
speeds = np.random.normal(mean_speed, std_dev_speed, num_particles)
radii = np.random.normal(mean_radius, std_dev_radius, num_particles)

# Add the moving particles
for i in range(num_particles):
    moving_particles.append(create_particle(speeds[i], radii[i]))
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Simulation
#-------------------------------------------------------------------------------
fig = plt.figure(num="DLA Simulation")
while plt.fignum_exists(fig.number):
    if len(moving_particles) > 0:
        steps += 1
    moving_particles = move_particles() 
    plot_particles()
plt.show()
#-------------------------------------------------------------------------------
