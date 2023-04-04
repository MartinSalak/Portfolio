# Architectural lessons from Topology Optimization

This page shows my work for my graduation, where I researched the topic of Topology Optimization and its possible implementations in architecture. This research was performed under the supervision of P. Nourian (_Faculty of Architecture and the Built Environment_) and M. Langelaar (_Faculty of 3Me_). A link to my work in the TU Delft repository can be found [**here**](http://resolver.tudelft.nl/uuid:5dc60528-701c-496c-90a2-a804d7a7aada).

This page will summarize my process and the results of this research. 

### Topology Optimization

Topology Optimization is a mathematical approach in designing geometry, where the volume is minimized, while still reaching a high stiffness. It is often used in mechanical- and aerospace engineering to optimize parts so they require less material. The main idea behind Topology Optimization is to calculate what voxels (or pixels) in an element are important for the stiffness and what voxels can be removed. Because there is no preconceived shape, Topology Optimization can create innovative and high-performance shapes.

### Previous research

The method (originally developed by Bendsøe and Sigmund, 2004) starts with dividing the design space in pixels and preparing the supports and loads, the steps can also be seen in figure (1). Each pixel has 4 nodes that can move in both X and Y directions (their Degrees of Freedom). Supports are defined as nodes that cannot move in one or all of the directions. With the supports and forces defined, the displacements of the nodes can be calculated using Finite Element Analysis (FEA). In optimization problems, the system is always minimized, so instead of maximizing the stiffness, the compliance is minimized. Compliance is a value of how much the whole system can move and is written as in formula (1). To find the optimal solution, Gradient Descent is used, which will set values of each element based on its derivative of the compliance. 

![Compliance](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_1.png?raw=true "Compliance")
[Larger image](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_1.png)

### Topology Optimization in Architecture

Deetman (2019) made changes to this standard method of Topology Optimization to make it work in 3D. The structure of the algorithm and main calculations are identical, but Deetman added a new Stiffness matrix for the FEA and used a strict numbering system. Each element now has 8 nodes and each node can move in 3 directions, which causes the matrices to increase by a power of 3. Furthermore, loads and supports can be set in the same way and this can already be used in some models. Figure 2 shows a simplified version of the QNCC building, where Topology Optimization is used to generate columns for the big roof.

![3D Results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_2.png?raw=true "3D results")
[Larger image](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_2.png)

Architectural cases are different than cases in mechanical engineering, as most forces are not predefined, but are dependent on the density. The geometry itself is heavier than the forces on the system. Placing a voxel will then not always improve the stiffness, as it can also generate new forces which decrease the stiffness. Density-dependent forces have to be added to the algorithm, which are new forces on nodes where the element around the nodes exists. This can be mathematically described using a sigmoid function, which sets values to either 0 or 1, depending on x. 

Previously the force was pre-set and constant, so the derivative of the force was 0. Now that the force is dependent on the density, the force should be included in the derivative. Rewriting the derivative previously found, but without the derivative of the force being 0, gives a new derivative, as shown in figure 3. When looking at this formula, one can see that the result of the derivative is no longer always negative. The graph of an element is no long monotonic and therefore Gradient Descent can no longer be used. Therefore, another optimizer is implemented, called the Method of Moving Asymptotes. 

![Selfweight](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_Self.png?raw=true "Selfweight")
[Larger image](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_Self.png)

Solving the system now still relies on a pre-set force, which is usually not the case in architectural problems. Instead of forces in the system, the constraint in architectural cases is that each void has to have a roof over itself. Or in other words; for each void, the sum of the elements in the column (above the void) should be larger than 1. Figure 4 shows the roof-constraint working, for columns where there is no roof, the algorithm sets all the values above the void to a certain value. Note that “grey” values are punished, in order to get black and white results. Another sigmoid function is implemented to that values that are grey are considered as a 0 and values close to 1 are considered as a 1. After many iterations, the shape is black and white and shows to be a dynamically relaxed structure. However, it is very thin. The optimizer will minimize the volume and increasing the thickness will generate more forces.

![Roof Constraint](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_Optimization.gif?raw=true "Project examples")
[Larger image](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_Optimization.gif)

When this would be built, it can be quickly seen that any forces on the roof will make it collapse. Another type of density-dependent force has to be added, namely a snow-load on the roof. This is a force that is placed on the highest element, which should make sure the roof will not fail when forces are placed on it. In other words; the element will gain a force if the sum of the elements above it, including itself, is equal to 1.  This can be mathematically described using a smooth-Heaviside function, which sets the value of y to 1 if the value of x is in the range of 0.5 and 1.5. 
The total forces can now be written as the sum of a preset force, the self-weight and the snow-load.  Summarized, the derivative can then be taken from this and written as follows:

![Selfweight results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_Total.png?raw=true "Project examples")
[Larger image](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_Total.png.png)

### Results

In order to validate these findings, several tests were performed, comparing the algorithm to the Topology Optimization software in Ansys. The results were comparable for simple topology optimization problems, but no roof-constraint could be added in the software. It requires further research to find the feasibility of the generated geometry. However, the results of certain configurations can be analyzed and compared with existing architecture. 

During this research, toy-problems were used to solve each step. Figure (6) shows the result of 2 configurations of the final toy-problem when all the constraints were added. It shows the roof-constraint being added to columns where voids exist and this causes a roof to be created. The main problem in this configuration is the number of voxels. Doubling the resolution of a 3D problem will increase the Degrees of Freedoms by a factor of 8. The most time the algorithm takes is spent solving the system and calculating the displacements. The result of the toy-problem shows the generation of dome-like structures above the large void and the beginning of arc-like shapes above the doors. Domes and arches are very common in masonry structures, as they allow for building materials with high compression, but low tension quality. These results show that the algorithm generates geometry that is representing some elements in architecture. However, it can be noticed that the dome has many inaccuracies, due to the resolution and specifics in the optimization process. Another conclusion that can be drawn is that cubic voids (currently the only possibility in the algorithm) are a poor choice to use, as they don’t allow for optimal geometry.


![Selfweight results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_3.png?raw=true "Project examples")
[Larger image](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_3.png.png)

Lastly, the question is “Can we generate buildings using Topology Optimization?”. To reflect on this research, an example of a Bauhaus building is taken as a configuration and its geometry can be seen in figure (7). The configuration of the Haus am Horn is used as input, together with the roof-constraint and self-weight. The results are still quite poorly because of calculation time, but still, some shapes can be seen. Walls in between rooms are always built, as they are needed to carry the roofs. However, above all the doors, arches are generated to save material above them. The rooms themselves all are generated with domes above them, with the large room having a high dome. The section of the dome is promising and hallways are starting to be optimized, however being subject to the low resolution.

![Selfweight results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_4.png?raw=true "Project examples")
[Larger image](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Graduation_4.png.png)

To answer the question; yes, I think eventually Topology Optimization could be a method to generate buildings. The process itself is directly linked to FEA and even more constraints can be added in the optimizer itself. However, more research is needed to generate more useful geometry, that eventually could also be tested in structural calculations. Essential for this useful geometry is the resolution of the voxels, which allow for more accurate results. One other constraint the could be added is the voids themselves; instead of starting with cubic voids, starting with a 2D layout could be better, where the optimizer is allowed to generate optimal voids as well. Concluding, Topology Optimization could generate and shape masonry architecture, but a higher resolution and more optimal configurations are needed.

[**Back to main**](https://rickvandijk1.github.io/PortFolio)
