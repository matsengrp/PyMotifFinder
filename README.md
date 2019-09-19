# PyMotifFinder

[![Docker Build](https://img.shields.io/docker/build/jfukuyama/pymotiffinder.svg)](https://hub.docker.com/r/jfukuyama/pymotiffinder/builds)

PyMotifFinder is a python implementation of MotifFinder, an algorithm for identifying mutations that are potentially due to gene conversion.

The package depends on `numpy`, `pandas`, and `biopython`.
To install the package, clone the repository and run
```
pip install .
```
from the base directory of the package.

The package can also be used inside a docker image.
You can pull from dockerhub, with
```
docker pull matsengrp/pymotiffinder
```

Once the image is pulled, calling
```
docker run matsengrp/pymotiffinder
```
will run all the tests, and
```
docker run -it matsengrp/pymotiffinder /bin/bash
```
will run an interactive session with PyMotifFinder installed.
