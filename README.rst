FFANS
=====
.. image:: https://sites.google.com/site/btanygin/_/rsrc/1430314739953/research/physics/simulation/ferrofluids/doc-images/logo-header-48px.png
.. image:: http://img.shields.io/travis/psci2195/ffans/master.svg?style=flat
        :target: http://travis-ci.org/psci2195/ffans
.. image:: http://img.shields.io/badge/license-GPL-yellow.svg?style=flat
        :target: https://github.com/psci2195/ffans/blob/master/LICENSE.txt
.. image:: http://img.shields.io/badge/arXiv-1503.05854-orange.svg?style=flat
        :target: http://arxiv.org/abs/1503.05854

Ferrofluid Aggregates Nano Simulator (FFANS) is a ferrofluid simulation platform.

About
-----
The FFANS is a molecular dynamics and dipole-dipole micromagnetics based simulation package. It is powered by an interactive orthographic and perspective 3D presentation layer with the automatic screenshoting and state saving allowing researchers to investigate physical process continuously in an obvious and visual way. High numerical stability gives an ability of a long-term run. Atomic object of simulation is a nanoparticle. Implicit solvent is modeled by Langevin thermostat. The Langevin dynamics has been implemented by Analytical Dissipative Integrator approach. Units are given in SI.

Product is written in C/C++ (physical process simulation) and OpenGL (3D graphics). Currently supported OS is only Windows. The cross-platform Qt layer implementation was started within other `project <https://github.com/psci2195/qt-ffans>`_. Also, the key feature of the given project (Analytical Dissipative Integrator) is being merged to other `project <https://github.com/psci2195/espresso-ffans>`_. However, the visualization capabilities are still very flexible and rich in the present project and its Qt port. If you have questions regarding usage and/or collaboration, feel free to `contact us <b.m.tanygin@gmail.com>`_. Original publication for `citation <http://cpb.iphy.ac.cn/EN/abstract/abstract65596.shtml>`_: ::

  @article{Tanygin2015,
  author = {Tanygin, B. M. and Shulyma, S. I. and Kovalenko, V. F. and Petrychuk, M. V.},
  eprint = {1503.05854},
  journal = {Chinese Physics B},
  number = {10},
  pages = {104702},
  title = {{Ferrofluid nucleus phase transitions in an external uniform magnetic field}},
  volume = {24},
  year = {2015}
  }

.. image:: https://sites.google.com/site/btanygin/_/rsrc/1430314550490/research/physics/simulation/ferrofluids/ffans-downloads/FFANS-logo192px.png

License
-------
Copyright (C) 2006,2011,2013-2017 Dr. Bogdan Tanygin, Mr. Dmytro Matskevych, present repository contributors, and authors specified at the source files beginning.

The FFANS is a free software made available under the GPL License. For details see the LICENSE file.
