# mr-hrtf
A simple 3d sound head related transfer function (HRTF) library implementation.

A HRTF filter can be used to simulate the direction a sound is coming from.
That is if used correct you can close your eyes and hear where a sound is
coming from.

This uses data from the CIPIC database [1], and fast fourier code from the
kiss_fft library [2] to implement a hrtf filter.

It is presented as a simple to use C interface.

## Simple code example

TODO: (sorry)

## Dependencies

* Python 3 + scipy
* Cmake 3.3 or later
* C++11 compatable compiler

## Getting ready

### Getting the CIPIC data

Download the HRTF data from the CIPIC website [2] and save it somewhere.

### Creating the data files

Use `generate_hrtf_database.py` to process the CIPIC matlab data to a data
format the library can use.

### Building

```
    cd mr-hrtf
    mkdir build
    cd build
    cmake ..
    make
```

## This is a simple library

* It doesn't use a real-only FFT
* It doesn't use SIMD
* It doesn't interpolate between HRTFs
* It uses float based sample mixing
* It uses simple linear blending when crossing HRTF boundaries

## License
This code is licesned under the AGPL V3 license, execpt for files for kiss_fft,
which are licensed under the BSD license.

[1]: http://interface.cipic.ucdavis.edu/sound/hrtf.html
[2]: http://kissfft.sourceforge.net/
