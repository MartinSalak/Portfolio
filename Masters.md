# Masters

_During my Master's Building Technology I focussed largely on Structural Design and Computational Design. On this page, I will go over some projects and what I learned in these projects._

## Earthy

For this project I documented everything at another MarkDown website, click here to visit [**this**](https://rickvandijk1.github.io/Earthy/) site. 

In Earthy students will dive into generating geometry for masonry buildings in refugee camps. In each aspect of the design of the project, the computational power of GrassHopper is used. In the early stages of the project, this meant finding optimal locations for the project and filtering out locations on basis of preset preferences. My team and I developed a design for a Hammam in the camp, so we analyzed other Hammam and using bubble diagrams tried to find an optimal layout. After taking on the challenge of using hexagonal grids, we could map the optimal layout to the grid and therefore had our location and program requirements.

![Rulesets](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Ruleset-Gif.gif?raw=true "Ruleset-Gif")

To translate this to a buildable Hammam, the layout was extruded to make usable rooms. In this step, not only the shape was important, but also how this would be constructed on site. After experimenting with rectangular bricks, it was found that triangular bricks would support the hexagonal grid better. Makeability was important as the Hammam should be constructed by local people. A method was designed for how the walls and its openings could be made, more info in [**this**](https://rickvandijk1.github.io/Earthy/walls/) and [**this**](https://rickvandijk1.github.io/Earthy/openings/) link. The domes were more challenging, as they should be dynamically relaxed and makeable. To make the domes, we followed the traditional method of making domes with muqarnas, small elements in similar shapes that make the dome. The elements can be seen in the picture below, which was a model we made. Lastly, we made a second roof to improve its structural performance, as well as to improve the rain flow from the roof. 

There were a few other interesting features of our Hammam that were developed by dynamically relaxing the shape. The stoa around the garden was made by dynamically relaxing the shape with points at the garden and at the wall. The thickness of the stoa was dependent on the stresses in the material. Secondly more information on the systematic generation of our shape using pseudocode, more information can be found [**here**](https://rickvandijk1.github.io/Earthy/forming/).

![Results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Hammam1.jpg?raw=true "Hammam")

![Results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/section_hammam1.jpg?raw=true "Hammam section")

## Robotic Manufacturing

In the spring of 2019, I followed Design Informatics and 1:1 Interactive Architecture Prototypes Workshop, two courses focusing on robotic manufacturing of geometry. In these courses, we explored the possibilities of using robots in architecture. In the Prototypes Workshop, I focussed on making a pavilion out of wood, generating the shape, and generating toolpaths for the robot. The shape can be seen in the picture below, which was based on light, rainwater, and forces. We were able to test out our toolpaths in the Hogeschool van Amsterdam and learned a lot about how robots work and wood as a material for this. More info can be found [**here**](http://100ybp.roboticbuilding.eu/index.php/Msc2G2:Main)

![Results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Wood.png?raw=true "Wood")

During the same time, I followed the course of Design Informatics, where the focus was again on robotic manufacturing. In this course instead of removing material, additive manufacturing was used. In combination with the University of Prague, we developed a technique where hemp rope is infused with bone glue and then placed on panels. Using a robot arm the shape was generated and could then be dried. After a week of intense work, we managed to create a dome. In this project, I was responsible with another student for any movement the robot made, from toolpaths to managing it on runtime. More information on this project can be found here. More info can be found [**here**](https://tudelftaet.wordpress.com/2019/04/12/dreamcatcher-pavilion-unveiled/?fbclid=IwAR1C8cO_yEERQ9zClqRViW4g9BaVEHM7qKbFZeLZ_iMqZviT2Gz4BRT33oE)

![Results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Dreamcatcher.jpg?raw=true "Dream")

![Results](https://github.com/RickvanDijk1/PortFolio/blob/gh-pages/assets/img/Dreamcatcherme.jpg?raw=true "Dream")

## 3D Computer graphics

Following my personal interests in game design, I followed a course at the faculty of EWI, on the subject of 3D Computer Graphics. The course consisted of in-depth explanations of how 3D graphics are generated and how shaders are written. Together with some other students, we created a small game from scratch, using C++ and the OpenGL library. In this game, the player rides his car over hills and tries to hit enemies. We wrote everything ourselves, from shaders, and terrain generation to game logic. In this project I learned how 3D graphics are generated and how to work in informatics projects, using programs like Git. In my free time, I took this knowledge and further kept developing small games for fun.

## Other Courses

During my Master's I followed other courses on structural design, building physics, and architectural design. In Extreme I designed a new library for the country of Sint Maarten, taking into account extreme weather conditions. In Bucky Lab, we created an alternative way of heating a hospital using old CPU coolers. This turned out to be not very effective, but the advantage was that rooms could be heated locally, resulting in lower heat demand. 

[**Back to main**](https://rickvandijk1.github.io/PortFolio)
