# Flagser-laplacian

Copyright © 2023 Benjamin Jones.

### Description

Flagser-laplacian computes the persistent spectra of directed flag complexes using a persistent topological Laplacian.

This includes substantial portions of the flagser package by Daniel Lütgehetmann, which is available at https://github.com/luetge/flagser.

### Building

Flagser-Laplacian requires a C++17 compiler, cmake, and a MATLAB installation. Here is how to obtain, build, and run flagser-Laplacian:


### Running

To call flagser-laplacian you run

```sh
./flagser-laplacian --out-prefix testprefix test/a.flag
```

For more detailed instructions, see `docs/flagser-laplacian-documentation.pdf`. 

Included in this package is also the original program `ripser`, modified only
in that the features of `flagser` are supported.


### License

flagserlaplacian is licensed under the MIT license (COPYING.txt), with an extra clause (CONTRIBUTING.txt) clarifying the license for modifications released without an explicit written license agreement.
