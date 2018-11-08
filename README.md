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

Or... check out the Web application (locally for now) with

```
python App/csp_app/app_csp.py
```

And then navigate to `http://0.0.0.0:5000/` in your browser.

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

Here is a demo of the software:

```
# Imports
from cone_search_plus import csp
import astropy.units as q

# Coordinates of Trappist-1
ra = 346.6223683553692
dec = -05.0413976917903

# Make the SourceList object
sl = csp.SourceList([ra,dec], 2*q.arcmin)
```

```
  _r     USNO-A2.0    RAJ2000    DEJ2000   ACTflag Mflag Bmag Rmag  Epoch  
 arcm                   deg        deg                   mag  mag     yr   
------ ------------- ---------- ---------- ------- ----- ---- ---- --------
1.7118 0825-19856975 346.595528  -5.031456               18.3 17.2 1953.680
0.8099 0825-19857185 346.610139  -5.035592               17.4 15.5 1953.680
0.1138 0825-19857340 346.621587  -5.039667               18.8 18.2 1953.680
1.9556 0825-19857424 346.627275  -5.073620               18.7 17.7 1953.680
1.1837 0825-19857643 346.641906  -5.044656               17.9 17.0 1953.680
1.6607 0825-19857776 346.650070  -5.039175               15.2 14.7 1953.680
6 sources found within 2.0 arcmin
```

And plot it:

```
sl.proximity_plot()
```


![png](figures/csp_demo.png)

Nice!

## Licensed

This project is Copyright (c) Joe Filippazzo and licensed under the terms of the BSD 3-Clause license. See LICENSE for more information.
