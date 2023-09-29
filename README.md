# Flagser-laplacian

Copyright © 2023 Benjamin Jones.

### Description

Flagser-laplacian computes the persistent spectra of directed flag complexes using a persistent topological Laplacian.

This includes substantial portions of the flagser package by Daniel Lütgehetmann, which is available at https://github.com/luetge/flagser.

### Building

Flagser-Laplacian requires a C++17 compiler, cmake, and a MATLAB installation. 

1. You should clone or download this git repository: ```git clone https://github.com/bdjones13/flagser-laplacian```. 
2. You May need to edit CMakeLists.txt line 21 and possibly lines 27 and 33 to point to your MATLAB installation.
3. In the repository root, create and enter a directory ```build```.
4. From within ```build``` run ```cmake ..``` followed by ```make```.

This will produce an executable ```flagser-laplacian```. If you want to use flagser-laplacian from a different directory, you may want to add the directory containing the ```flagser-laplacian``` executable file to your path.

### Running

To call flagser-laplacian you run

```sh
./flagser-laplacian test/a.flag
```

For more detailed instructions and available options, see `docs/flagser-laplacian-documentation.pdf`. 

Included in this package is also the original program `ripser`, modified only
in that the features of `flagser` are supported.


### License

flagser-laplacian is licensed under the MIT license (COPYING.txt), with an extra clause (CONTRIBUTING.txt) clarifying the license for modifications released without an explicit written license agreement.

### Citation

When the paper for this software is available on the Arxiv and/or published, we will provide an appropriate bibtex entry for those who would like to cite this software.