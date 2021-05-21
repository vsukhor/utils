#  utils
A library of loosely related constructs of general use.  
Intended mainly for the use in other repositories in this collection.

## Overview

Functionality of this library can be best summarized by listing the namespaces provided 
(see doxygen documentation for more details):

#### Arrays 
Lightweight  array types implementing convenient arithmetics as well as some functionaity commonly used in geometric applications.

#### Biochemical
Reader and writer for files in Protein Data Bank format.

#### Common
General-usability typedefs and functions.

### [Config](utils/config/readme.md) 
A simple configuration flie reader.

#### Graph
Abstruct graphs and some common operations on them.

#### Random
Pseudo-random number factories.

## Installlation

The library contains both the classes that require separate compilation, and the header-only includes.
No external dependencies are necessary in order to build the library. 
However, random number factories (despite being header-only), 
would need the availability of [boost](https://www.boost.org/) (only the headers) or [NVIDIA cuRAND](https://developer.nvidia.com/curand), 
depending on the factory type chosen, in the system where an executable using these classes is built.

The library can be compiled with a C++17-capable compiler (e.g. on macOS, either gcc 7.3.0 or clang 10.0.0 would be sufficient).
There are two ways for building the library:  

1. On systems having cmake (ver. 3.15 or higher) installed, the build can be performed e.g. as follows:  
    * cd utils
    * mkdir build
    * cmake -S . -B build
    * cmake --build build
    
2. cmake-independent Makefile is available for for a more targeted manulal build. 
    Please see beginning of the Makefile for instructions.


