#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='PyMotifFinder',
      version='1.0',
      description='Functions for identifying gene conversion',
      author='Julia Fukuyama',
      author_email='julia.fukuyama@gmail.com',
      url='http://www.github.com/matsengrp/PyMotifFinder',
      packages=find_packages(exclude='test')
)
