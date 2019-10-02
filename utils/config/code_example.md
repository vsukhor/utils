# Example of the source code 

Consider the following code. It defines parameters
- 'init_mass'  
- 'assign_separately'
- 'spring_k'   
of different types, and uses 'Reader' functor provided by the library to import their values from a config file.
 
```
#include <array>
#include <numeric>
#include "config/reader.h"

class ConfigReader {

	using Reader = Config::Reader;
	
	// Definitions of acceptable range boundaries. 
	// (Alternatively, these and other definitions may be included from the 'utils' library.)
    
	template <typename T> constexpr auto huge {std::numeric_limits<T>::max()};
    template <typename T, auto N> constexpr std::array<T,N> zeros = {};
    template <typename T, auto N> constexpr std::array<T,N> huges = {filled_array<T,N>(huge)};
    template <typename T> constexpr std::array<T,2> zerohuge = {static_cast<T>(0.L), huge};
	constexpr std::array<bool,2> bools = {false, true};

public:

	const Reader read;			// The library-supplied reader.

	// The parameters.
	
    const std::array<double,3> init_mass;	// Initial particle masses, by type.
    const bool 	assign_separately;	    // Assign each type separatley.
    const real	spring_k;				    // Spring coefficient.

	// Constructor.
	
    ConfigReader( const Reader& read )
	: read {read}
	, init_mass {read("init_mass", {zeros<double,3>, huges<double,3>})}
	, assign_separately {read("assign_separately", bools)}
	, spring_k {read("spring_k", zerohuge<real>)}
{}

};
```
Here, **read()** takes two arguments: the parameter name and the allowed boundaries on its value.


