# natedrop.jl


**"Simulate 3D printed droplet arrays, quickly and easily"** 

natedrop.jl is a **particle simulation** for printed, liquid coated aqueous droplets that adhere to each other and form stable lipid bilayers. These droplets can be arranged in networks that display large scale emergent behavior. Currently the modelling of the adhesion between droplets, which is approximated by a linear spring of stiffness, k. Additionally, water can flow through the bilayers from an area of low osmolarity to an area of high osmolarity through the circular area of intersection between droplets. Both of these relationships can, in tandem, lead to dynamic folding and unfolding of droplets arrays, which can form shapes not possible through simple printing. The folding and unfolding behavior can also lead to future application in the field of soft robotics. 

## How to use 

each droplet in natedrop is represented as a row in a matrix. Each row will have eight columns, and each column will represent a specific parameter:
Column 1: Radius 
Column 2: Osmolarity
Column 3: X-position
Column 4: Y-position
Column 5: Z-position 
Column 6: X-velocity
Column 7: Y-velocity
Column 8: Z-velocity

Additionally, the natural length between each droplet, **L** needs to be specified for the function. A lower **L** value will have the droplets stick closer together. 
A recommended L-value is 0.8  
The bilayer attraction is represented as a linear spring of stiffness **K**.
A recommended K-value is 1000.0
To make the diffusion between droplets begin when the velocities have hit zero, set boolean_osm to [0.0] 



