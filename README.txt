##############################################################
###   C++ program of vertex-spring method for simulation   ###
###   of 3D-morphogenesisthrough nonuniform sheet growth   ###
##############################################################
Written by Hiroshi C. Ito (2024)
Email: hiroshibeetle@gmail.com

This directory contains three C++ program files used for producing numerical results in a paper titled "Growth regulation bringing modularity to morphogenesis of complex three-dimensional exoskeletons" (see the end of this file for the abstract):
vertex1.cc (for figure 6)
vertex2.cc (for figure 7)
vertex3.cc (for figure 8)

See the description below for execution of these programs.

<Tested environment>
OS: Utunbu 20.04
Application: g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0, mfem-3.3.2, yorickvis-0.2, rlwrap-0.43

OS: MacOS (12.6.9)
Application: g++ (Apple clang version 13.1.6), mfem-3.3.2, XQuartz-2.8.5, yorickvis-0.2, rlwrap-0.46.1

<Requirement>
Execution of the programs requires C++ compiler.
Visualization of the simulation output requires a visualization tool "yorickvis" written in Yorick language, and rlwrap.

Under MacOS, XQuartz is also required for visualization.


<Local install of yorickvis (not needed if already installed somewhere)>
yorickvis can be quickly installed under armorbuilderVERTEX/
by following commands (using git):

cd armorbuilderVERTEX/
git clone https://github.com/yorickuser/yorickvis.git
cd yorickvis/
./install.sh
cd ../

Install of yorickvis produces a symbolic link "Yorick" at home directory, by which you can easily run yorickvis by typing "~/Yorick/yorickvis".

yorickvis is uninstalled by removing "~/Yorick" and "armorbuilderVERTEX/yorickvis".

<Compilation and execution (in the case of vertex2.cc)>

g++ vertex2.cc
./a.out -o output

The state of simulation is recorded in three data files:
output.dat (data of sheet shape)
output_met.dat (data of growth metric)
output_outcount.dat (count of output)

State of the ongoing simulation is visualized simultaneously by running yorikvis at another console:

cd armorbuilderVERTEX/
~/Yorick/yorickvisf output

Pressing the Mouse-Left button at the right-top "Start/End" in the window titled "Yorick 0" starts following the simulation.
Pressing the Mouse-Right button at the right-top "Start/End" or "Pause/End" ends visualization, entering yorick-prompt.
Typing "quit" ends yorick, returning to the normal console prompt.




################################################
###           Paper information              ###
################################################

Title: Growth regulation bringing modularity to morphogenesis of complex three-dimensional exoskeletons

Author: Hiroshi C. Ito and Yu Uchiumi

DOI: 10.1098/rspb.2024.1943

Abstract:
Diverse three-dimensional morphologies of arthropods’ outgrowths, including beetle horns, are formed through the non-uniform growth of epidermis. Prior to moulting, epidermal tissue peels off from the old cuticle and grows non-uniformly to shape protruding structures, which are often branching, curving, or twisting, from the planar epidermis. This non-uniform growth is possibly regulated by the distribution of morphogens on the epidermal cell sheet. Previous studies have identified molecules and signalling pathways related to such morphogenesis; however, how local regulation of cell sheet growth can transform planar epidermis globally into complex three-dimensional structures, such as beetle horns, remains unclear. To reveal the relationship between epidermal growth regulation and generated structures, this study theoretically examined how various shapes can be generated from planar epidermis under a deductive growth model that corresponds morphogen distributions to non-uniform growth on tissue. The results show that the heterochronic expression of multiple morphogens can flexibly fuse multiple simple shapes to generate various structures emulating complex outgrowths of beetles. These findings indicate that morphogenesis through such a mechanism may have developmental stability and modularity, providing insights into the evolution of the diverse morphology of arthropods.
