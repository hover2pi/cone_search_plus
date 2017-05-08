#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
try:
    from setuptools import setup, find_packages
    setup
except ImportError:
    from distutils.core import setup
    setup

setup(
    name='Cone Search Plus',
    version='0.1.0',
    description='Search the sky for targets using a variety of tunable constraints',
    url='https://github.com/hover2pi/cone_search_plus',
    author='Joe Filippazzo, Neil Zimmerman',
    author_email='jfilippazzo@stsci.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
    ],
    keywords='astrophysics',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['numpy','astropy','matplotlib','ephem','astroquery'],

)