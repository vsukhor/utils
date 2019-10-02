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

The library contains both classes require building and header-only includes.
No external dependencies are required in order to build the library. 
However, random number factories, which are themselves header-only, 
would need [boost](https://www.boost.org/) (only the headers) or [NVIDIA cuRAND](https://developer.nvidia.com/curand), 
depending on the factory type chosen.

There are two ways for building the library:  
* Using cmake, which also creates documentation if doxygen (ver. > 3.14) is installed.
* With a simple Makefile for a more direct manual build. Please see beginning of the Makefile for instructions regarding the variables required.
A C++17-capable compiler is required (e.g. on macOS, either gcc 7.3.0 or clang 10.0.0 would be sufficient).

