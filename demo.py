import pymunk
import pymunk.matplotlib_util
import matplotlib.pyplot as plt
import numpy as np

# Initialize space
space = pymunk.Space()
space.gravity = (0, -900)  # gravity in the downward direction

# Create mass
mass = 1  # mass of 1 kg
radius = 10  # radius of the mass
moment = pymunk.moment_for_circle(mass, 0, radius)  # moment of inertia for the mass
body = pymunk.Body(mass, moment)  # create the body
body.position = (0, 0)  # position the body at the origin
shape = pymunk.Circle(body, radius)  # create the shape
space.add(body, shape)  # add the body and shape to the space

# Create spring
spring_anchor = (0, 100)  # position of the spring anchor
spring = pymunk.DampedSpring(body, space.static_body, (0, 0), spring_anchor, 100, 1000, 0)
space.add(spring)  # add the spring to the space

# List to store positions
positions = []

# Create a matplotlib figure and axes
fig, ax = plt.subplots()

# Create a DrawOptions object
draw_options = pymunk.matplotlib_util.DrawOptions(ax)

# Simulation loop
for i in range(1000):
    space.step(0.01)  # step the space forward in time
    positions.append(body.position)  # store the position

    # Clear the axes
    ax.clear()

    # Draw the space
    space.debug_draw(draw_options)

    # Plot the position of the mass
    positions_array = np.array(positions)
    ax.scatter(positions_array[:, 0], positions_array[:, 1], color='red')

    # Update the figure
    fig.canvas.draw()

plt.show()
