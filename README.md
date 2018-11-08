# cone_search_plus

[![Build Status](https://travis-ci.org/hover2pi/cone_search_plus.svg?branch=master)](https://travis-ci.org/hover2pi/cone_search_plus)
[![Coverage Status](https://coveralls.io/repos/github/hover2pi/cone_search_plus/badge.svg?branch=master&service=github)](https://coveralls.io/github/hover2pi/cone_search_plus?branch=master&service=github)
[![Documentation Status](https://readthedocs.org/projects/cone_search_plus/badge/?version=latest)](https://cone_search_plus.readthedocs.io/en/latest/?badge=latest)

A whiz-bang Python package for all-sky cone searches with tunable constraints.

Requirements:
- numpy
- astropy
- matplotlib
- ephem
- astroquery

Or... check out the Web application at [https://csp.stsci.edu](https://csp.stsci.edu)

## Installation

Install via PYPI with

```
pip install cone_search_plus
```

or via Github with

```
git clone https://github.com/hover2pi/cone_search_plus.git
python cone_search_plus/setup.py install
```

## Documentation

Full documentation for the latest build can be found on [ReadTheDocs](https://cone_search_plus.readthedocs.io/en/latest/).

The package also contains detailed Jupyter notebooks highlighting the core functionality of its primary classes, including

- [cone_search_plus.csp.SourceList](https://github.com/hover2pi/cone_search_plus/blob/master/notebooks/csp_demo.ipynb)

## Demo

Here is a demo of the software with some code snippets:

```
from cone_search_plus import csp
mc = csp.SourceList('Great package!')
print(mc.words)
```

And at least one nice plot:

![png](figures/plot.png)

## Licensed

This project is Copyright (c) Joe Filippazzo and licensed under the terms of the BSD 3-Clause license. See LICENSE for more information.
