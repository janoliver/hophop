==================================
 A monte carlo hopping simulation
==================================

About
=====

This software simulates charge carrier movement through disordered
 systems. It is capable of using the variable range mechanism as well as
 some lattice structure. The energetic disorder can be specified as an
 input parameter.
The software determines some key features of the charge transport such as
 mobility, diffusivity, certain important energies and so on. 

compilation
===========
Run the following commands to compile the program::

    cmake .
    make
    make install

The following libraries and software is required to compile successfully:

gsl
    The gnu scientific library provides many mathematical algorithms,
    functions, constants and so on. It is heavily used within the program,
    mainly for pseudo random number generation and probability
    distributions.

gengetopt
    The gnu getopt generator is used to provide nicely usable command line
    options and to generate the usage string of the software. (see ``Usage``)

gomp
    This is a openmp library for parallel computing. Parts of the software
    are parallelized which is why this is needed. This might become
    optional in the near future. 

Usage
=====

See ``hop --help``

Units
=====

Length
    The parameters specifying the length of the sample, are given in units
    of N^(-1/3). So choosing -l20 would result in 20*20*20 sites. N^(-1/3)
    is thus also the length scale in the simulation. Internally, the sample
    is split into cells of with 2 * --rc * --llength (cutoff radius *
    loc. length).  
    
    
Energies/Temperatures
    All energies are measured in units of the disorder parameter \sigma.
    Since:: 
    
        \sigma/kT \approx 3.0
        
    is a realistic value, the temperature T is usually around 0.3.

Coding style
============
    
We use a modified gnu coding style. The code can and should be formatted
using the indent tool like this: :: 
 
       indent -gnu -fc1 -i4 -bli0 -nut -cdb -sc -bap -l80

Block comments should look like this and precede functions etc.::

      /*
       * I am a block comment!
       */

For one-lined comments, use ``//`` 
