# Configuration file structure

Configuration files are expected to consist of independent text-formatted lines.
Blank lines and lines starting from the '#' character indicate start of a comment are ignored.
All other lines are treated as parameter records containing pairs 'parameter_name' 'parameter_value'.
The records are to follow the syntaxis:

parameter_name = parameter_vaules [# comment]
 
The field in brackets is optional.
'parameter_name' is extracted as the character sequence composed of all characters preceding the first space character. 
The equality sign should be followed by either one or several 'parameter_vaules' in decimal notation. 
In the former case, the parameter is read in as a scalar. 
Alternatively, array-type parameters (implemented as std::array) of predefined length are read in from a sequence of numerical parameter values.
The actual length and type of the parameter is set by that expected from the reader in ther source code rather than the configuration file itself.

### Examples

Legitimate records                             | Remarks
--------------------------------------------------------------- | -----------------------------------------------------
myValue = 30.2     						| scalar, has neither the description nor the  comment
edge_length (mm) = 3.e+4   # starting from the previous node	| scalar, has both the description and the comment
reaction_rates = 0.87 33.4 64.326 50.22 #(1/sec)		| array of length 4 (of real type)
dimensionality = 2 2  # in x- and y-directions respectivley	| array of length 2 (of integer type)

Illegitimate records  | Remarks
--------------------------------------------------------------- | -----------------------------------------------------
thisNumber =  77 items          			| missing number sign '#' separating the comment
a_parameter= 29.3					| missing  space delimiting the parameter name
threeValues as 3 numbers 53 82 09			|  missing equality sign preceding the value(s)



