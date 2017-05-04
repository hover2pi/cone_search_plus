#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
try:
    from setuptools import setup, find_packages
    setup
except ImportError:
    from distutils.core import setup
    setup

from codecs import open
from os import path

# Compile the query_2mass Fortran code
os.system('gfortran all_sky_target_tool/query_2mass.F -o all_sky_target_tool/query_2mass')

setup(
    name='All-Sky Target Tool',
    version='0.1.0',
    description='Search the sky for targets using a variety of tunable constraints',
    url='https://grit.stsci.edu/ins-tel/ote_target_tool',
    author='Neil Zimmerman, Joe Filippazzo',
    author_email='jcfilippazzo@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
    ],
    keywords='astrophysics',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['numpy','astropy','matplotlib'],

)