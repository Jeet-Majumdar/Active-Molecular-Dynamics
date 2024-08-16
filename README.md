# A Brief Introduction to Scalar Activity

Active molecular dynamics is an extension of traditional molecular dynamics that incorporates the effects of active particles or components. In this context, "active" refers to particles that can generate their own motion or forces, often consuming energy to do so. This approach is particularly useful for studying biological systems, soft matter, and other out-of-equilibrium phenomena.
When using scalar activity in active molecular dynamics:

Activity is represented as a scalar quantity associated with each particle or component in the system.
This scalar activity can influence the motion or behavior of particles, typically by modifying their velocity or applied forces.
The activity can be constant or time-dependent, and may vary between different particles in the system.
The equations of motion are modified to include terms that depend on the scalar activity, often in the form of an additional force or velocity component.
This approach allows for the simulation of systems with heterogeneous activity levels and can capture phenomena like collective motion, phase separation, or pattern formation in active matter.
Scalar activity models are often simpler to implement and analyze compared to more complex vectorial activity models, making them a good starting point for studying active systems.

By incorporating scalar activity into molecular dynamics simulations, researchers can investigate a wide range of active matter systems, from bacterial suspensions to artificial microswimmers and active colloids.

## How is scalar coming into the picture over here?
Here we use two temperature model to incorporate scalar activity into the picture.
Particles are divided into two class by defining two different temperatures. This brings two kinds of interaction into the picture:
1. Interaction between particles with same temperatures
2. Interaction between particles with different temperatures

The codes were developed by me to have a hands-on implementation of this simple two temperature activity simulation, and to study any induced phase separation. 
When all the particles are assigned the same temperature, it performs equilibrium dynamics.