# config
A simple configuration flie reader.

This library provides a simple interpreter for reading parameter names and numerical values from configuration files.

See also:

[- How to structure the configuration file?](conf_file_structure.md)

[- How to use the library in your code?](code_example.md)

### Installation

Although the library files are all headers, it makes use of [utils](https://github.com/vsukhor/utils) repository as a submodule. 
So, when cloning make sure to use '--recursive' option:

    git clone --recursive https://github.com/vsukhor/config

and, in the case of an update, do this explicitly on the submodule:

    git submodule init
    git submodule update


