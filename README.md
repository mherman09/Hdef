# Hdef

*Hdef* is designed to calculate displacements, strains, and stresses for various geologically relevant sources in an elastic half-space (hence, **H**alfspace **def**ormation). The currently available sources are:

- Point dislocations (shear and tensile)
- Rectangular dislocations (shear and tensile)
- Triangular dislocations (shear and tensile)
- Point volume change (Mogi) sources

*Hdef* also includes several programs that I wrote to complement and contextualize the half-space outputs (so in that sense, maybe *Hdef* stands for "**H**erman" **def**ormation). Actually, I got tired of looking up and recalculating things like plate velocity vectors, earthquake seismic moments, etc., so I wrote programs to have these results at my fingertips. These are included in *Hdef*, and I hope you find them useful as well.

I have used the codes in *Hdef* to make calculations for several peer-reviewed publications, and many other people have used *Hdef* for their own applications. However, I cannot guarantee that there are no bugs. If you think you have found a bug, please report it to me at matthew.w.herman@gmail.com


## Installation

Please see [INSTALL.md](INSTALL.md) for directions on how to download and install Hdef.


## Tutorials & Documentation

Tutorials can be found at my [website](http://www.matthewwherman.com/software.html).

Man pages for most of the programs and some of the scripts are in the man/ directory.


## License

*Hdef* is distributed under the MIT license. Please see [LICENSE](LICENSE) file.



## Citing Hdef

At a minimum, please cite the *Hdef* code with the following DOI: https://doi.org/10.5281/zenodo.3894137, as well as the following papers:

Herman, M.W., Furlong, K.P., Hayes, G.P., Benz, H.M. (2016). [Foreshock triggering of the 1 April 2014 Mw 8.2 Iquique, Chile, earthquake](https://doi.org/10.1016/j.epsl.2016.04.020). *Earth and Planetary Science Letters* 447, 119-129.


Herman, M.W., Govers, R. (2020). [Locating fully locked asperities along the South America subduction megathrust: a new physical inter‐seismic inversion approach in a Bayesian framework](https://doi.org/10.1029/2020GC009063). *Geochemistry, Geophysics, Geosystems* 21, e2020GC009063.

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Programs in Hdef:

COMPILED
	anneal_post            post-process an annealing run from fltinv
	clip                   determine whether x-y points are inside a polygon
	colortool              generate perceptually linear color maps
	dateutil               convert from/to dates and number of days
	distaz2lola            compute coordinates given distance/azimuth from a geographic point
	ff2gmt                 convert finite fault parameters to GMT-usable format
	fitutil                fit a set of x-y points using least-squares
	fltinv                 invert geodetic data for fault slip
	grid                   generate evenly spaced Cartesian and geographic points
	interpolate            interpolate between x-y points
	lola2distaz            compute distance/azimuth between coordinate pairs
	mtutil                 convert between seismic source characteristics
	numint                 integrate input function numerically
	o92util                compute displacement, stress, strain from rectangular fault source
	perspective            project 3-D object in perspective view
	platemotion            compute velocity for plate pair or Euler pole
	rangen                 generate random numbers
	readGCMT               read and parse GCMT catalog
	readkik                read strong motion in KiK-net format
	sphfinrot              make finite rotation on sphere
	stereo_project         compute forward and inverse stereographic projections
	tri_tool               manipulate triangular meshes
	triutil                compute displacement, stress, strain from triangular fault source
	vec2los                project displacement vector onto line of sight
	wraplos                wrap LOS displacement given observation frequency

BASH SCRIPTS
	coul_dip.sh            plot Coulomb stress changes on a dipping plane (map view)
	coul_hor.sh            plot Coulomb stress changes on a horizontal plane (map view)
	coul_xsec.sh           plot Coulomb stress changes on a vertical cross-section
	gmtcpt.sh              make GMT color palette from (nice) existing color maps
	insar.sh               plot line-of-sight and wrapped displacements (map view)
	simplify_ffm.sh        remove slip less than threshold value from FFM (pretty outdated...)
	surf_disp.sh           plot surface displacements (map view)
	ternary.sh             plot earthquakes on ternary discriminant diagram
	trg_schem.sh           plot target fault kinematics

INCOMPLETE OR OUTDATED
	eqempirical            implement earthquake empirical relation
	eventfrequency         count frequency of a list of values
	multifit               fit data series with multiple functions
	polyfit                fit data series with a polynomial function
	utm2geo                convert UTM to geographic coordinates

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Libraries in Hdef:

	algebra_module.f90
	annealing_module.f90
	calendar_module.f90
	earth_module.f90
	elast_module.f90
	eq_module.f90
	ffm_module.f90
	geom_module.f90
	interpolation_module.f90
	io_module.f90
	map_module.f90
	misfit_module.f90
	okada92_module.f90
	random_module.f90
	solver_module.f90
	test_module.f90
	tri_disloc_module.f90
	trig_module.f90 -->
