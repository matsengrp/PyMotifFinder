# PyMotifFinder

![Docker Build](https://img.shields.io/docker/build/jfukuyama/pymotiffinder.svg)

PyMotifFinder is a python implementation of MotifFinder, an algorithm for identifying mutations that are potentially due to gene conversion.

The package depends on `numpy`, `pandas`, and `biopython`.
To install the package, clone the repository and run
```
pip install .
```
from the base directory of the package.

The package can also be used inside a docker image.
To build the docker image, clone the repository, run
```
docker build . -t pmf
```
Once the image is built, calling
```
docker run pmf
```
will run all the tests, and
```
docker run -it pmf /bin/bash
```
will run an interactive session with PyMotifFinder installed.
