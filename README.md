#  utils
A library of loosely related constructs of general use. 
Intended mainly for the use in other ropositories in this collection.

## Overview

Functionality of this library can be best summerized by listing the namespaces available. 
See doxygen documentation for more details.

#### Arrays 
Lightweight  array types with convenient arithmetics as well as some functionaity commonly used in geometric applications.

#### Biochemical
Reader and writer for files in Protein Data Bank format.

#### Common
General-usability typedefs and functions.

#### Graph
Abstruct graphs and some common operations on them.

#### Random
Pseudo-random number factories.

## Installlation

The library contains both classes require building and header-only includes.
No external dependencies are required in order to build the library. 
However, random number factories, which are themselves header-only, 
would need [boost](https://www.boost.org/) (only the headers) or [NVIDIA cuRAND](https://developer.nvidia.com/curand), 
depending on the factory type chosen.

There are two ways for building the library:  
* Using cmake, which also creates documentation if doxygen (ver. > 3.14) is installed.
* With a simple Makefile for a more direct manual build. Please see beginning of the Makefile for instructions regarding the variables required.


