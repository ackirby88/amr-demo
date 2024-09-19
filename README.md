Finite Volume AMR Solver
=========
This is an demo program of solving the compressible Euler equations using a finite volume method on non-conforming adaptive meshes. The adaptive mesh refinement infrastructure is provided by [p4est](https://github.com/cburstedde/p4est).

---
# 1. Obtaining the DEMO  
    git clone https://github.com/ackirby88/amr-demo.git

---
# 2. Compilation
    To compile the demo and the 3rd party libraries (p4est, metis):  
    ./makescript.sh <options>

**Usage**: ./makescript.sh `{OPTIONS}`...`{COMPILER OPTIONS}`...`{DEMO OPTIONS}`...`{3PL OPTIONS}`  

| OPTIONS:             | Shortcut    | Description                                        |
|:---------------------|-------------|:---------------------------------------------------|
| --3pl                | -3pl        | build the 3rd party libraries: metis, p4est, gtest |
| --demo               | -demo       | build the demo solver                              |
|                      |             |                                                    |
| --help               | -h          | displays this help message                         |
| --clean              | -c          | removes local build directories                    |
| --distclean          | -dc         | removes builds and install directories             |
| --release            | -opt        | compile the project in optimized mode              |
| --debug              | -deb        | compile the project in debug mode                  |
|                      |             |                                                    |
| **COMPILER OPTIONS**:|             |                                                    |
| CC=`<arg>`           | cc=`<arg>`  | set the C compiler                                 |
| CXX=`<arg>`          | cxx=`<arg>` | set the C++ compiler                               |
| FC=`<arg>`           | fc=`<arg>`  | set the Fortran compiler                           |
|                      |             |                                                    |
| **3PL OPTIONS**:     |             |                                                    |
| --ALL3P              | -all3p      | compile all 3rd party libraries                    |
| --P4EST              | -p4est      | compile p4est                                      |
| --METIS              | -metis      | compile metis                                      |

# Recommended Build Options:
**Default Option (-go)**: Sets CC=mpicc CXX=mpicxx FC=mpif90  
    ./makescript.sh -go

**Intel MPI (-impi)**: Sets CC=mpiicc CXX=mpiicpc FC=mpiifort  
    ./makescript.sh -impi

---
# 3. Dependencies
    Software:
        CMake (Version 2.8 or newer): We also recommend ccmake or cmake-gui.
        MPI: optional MPI-3 one-sided functions available
        libtool: apt install libtool
        automake: apt install automake
        zlib: apt install zlib1g zlibc
        
    Documentation:
        doxygen: apt install doxygen
        graphviz: apt install graphviz
        latex: apt install texlive-full texmaker
---
# 4. License

The MIT License (MIT)

Copyright (c) 2024 Andrew C. Kirby

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
