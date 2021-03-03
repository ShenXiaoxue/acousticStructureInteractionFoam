# acousticStructureInteractionFoam
A time-domain solver for acoustic-structure interaction problem based on OpenFOAM.

However, this solver is not mature and it was my work during master's. And the two cases have not been validated. 

This program is implemented on ubuntu 14.04 based on the foam version foam-extend-3.1, there may be some bugs in other OpenFOAM versions, so foam-extend-3.1 is recommended.

In order to run the case, it is necessary to  build the solver first. cd acousticStructureInteractionFoam, and run "wmake".

There are two cases, plateCavity and stiffenedPlateCavity. The former one is simulation about the 200Hz 2Pa acoustic plane wave passing through a rectangular plate in time domain, while the latter one is spherical acoustic source generating acoustic waves and passing through stiffened plates. In the second case, I use python to obtain the nodes,blocks,etc,the mesh can be seen in "stiffenedPlateCavity/origin mesh", just copy it to "acoustic/constant".

