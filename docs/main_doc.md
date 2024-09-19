\mainpage Overview

The HPC AMR framework is constructed on top of the p4est framework. 
This framework is designed for arbitrary physics 
numerical discretizations that make use of 2D/3D structured Cartesian or unstructured 
quad/hexahedral elements. The goal of this work is provide automated parallel computing 
infrastructure and adaptive mesh refinement capabilities.


## Introduction
This is the introduction.

## Code Modules
### [AMR](group__amr__group.html) 
+ main program which choreographs all functions related to the adaptive mesh refinement data structure

### [PHYSICS](group__physics__group.html)
+ interface physics functions called from AMR code module

### [SOLVER](group__solver__group.html)
+ external example numerical discretization kernels that interface calls from PHYSICS code module
