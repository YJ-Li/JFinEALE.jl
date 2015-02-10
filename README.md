# JFinEALE

Finite element toolkit in Julia


![Alt Sample mesh](http://hogwarts.ucsd.edu/~pkrysl/site.images/ScreenHunter_31%20Feb.%2009%2020.54.jpg "JFinEALE.jl")

This toolkit is a redesign of the Matlab toolkit

[FinEALE](https://github.com/PetrKryslUCSD/FinEALE): Finite Element Analysis Learning Environment
  

JFinEALE is at this point very incomplete, at least compared 
to the original Matlab toolkit.  If something is missing, 
one can probably find it in the Matlab source code.

The author appreciates feedback and interaction/collaboration. 
Please don't hesitate to e-mail pkrysl@ucsd.edu.
 
## Get JFinEALE
 
Pkg.clone("https://github.com/PetrKryslUCSD/JFinEALE.jl")

## Testing

Pkg.test("JFinEALE")
[![Build Status](https://travis-ci.org/PetrKryslUCSD/JFinEALE.jl.png)](https://travis-ci.org/PetrKryslUCSD/JFinEALE.jl)

## Usage

The examples are available in their own repository (https://github.com/PetrKryslUCSD/JFinEALEexamples). It is not a Julia package.
In case you were wondering: Julia packages live in the .julia folder, whereas
the examples can be anywhere in your directory tree.

To run an example, "include" it (after using JFinEALE).
