# H2CO_QRG
Formaldehyde vibrational energy level calculation using Quasi Regular Grids


## SRC 
Fortran implementations for generating a QRG and and evaluating the 
associated vibrational eigenspectra. 


### Internal_Coordinates
Calculations are performed using 6 Degrees of Freedom defined using internal 
coordinates.

### Cartesian_Coordinates
Calculations are performed using 6 Degrees of Freedom defined using Cartesian
coordinates.

### Morse
Example calculations performed on the morse potential. Test case with 
analytic results for development and comparison. For a more detailed 
understanding of the morse model see our 
[previous work](https://github.com/swflynn/Quasi_Regular_Grids).

# External Codes
External Fortran modules were used during this project.
A detailed list is provided below along with links to the appropriate source code.

#### Quasi-Random Sequences:
The [Sobol Sequence Generator](https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html)
was used for the Quasi-Random Calculations (sobol.f90).
This code is under a [GNU LGPL License](https://www.gnu.org/licenses/lgpl-3.0.en.html)
and has been included here for convenience.

Special thanks to John Burkardt for taking the time to discuss various
quasi-random number generators and their implementation.

# Scientific Context
Accurate RoVibrational Energy calculations have a long history in computational
chemistry/physics.

The project was motivated in part by the work of
[Manzhos and Carrington](https://aip.scitation.org/doi/full/10.1063/1.4971295).
Special thanks to Tucker Carrington for discussing his results with us.

For an explination of Quasi Regular Grids see our 
[publication](https://aip.scitation.org/doi/full/10.1063/1.5134677).


## Authors
Shane W. Flynn, Vladimir A. Mandelshtam. 2019. UCI
