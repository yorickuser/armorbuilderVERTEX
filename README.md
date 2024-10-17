-------
GUI tool for visualization and simlation controrol written in [Yorick language](https://github.com/LLNL/yorick)
-------

* Directory "Yorick" contains additional functions for 2D- and 3D-plot and for event-driven control by mouse (function "idl"). Some functions are the same with those in "myYorick" in [my yorick fork](https://github.com/yorickuser/yorick)
* Currently, documents for "help, function_name" may be incomplete and source codes are not tidy, but most functions listed below have examples, shown by "help, function_name" at yorick prompt.
* For examples of the evnt-driven control function, see [my yorick website](https://translate.google.co.jp/translate?sl=ja&tl=en&js=y&prev=_t&hl=ja&ie=UTF-8&u=http%3A%2F%2Fbeetle.starfree.jp%2Fyorick%2Fyorick_basic_idl.html&edit-text=&act=url) (translated from japanese to english by google translation)

Requirement
-------
rlwrap

Usage
-------

* At yorickvis/ execute ./install.sh, which automatically get yorick-2.2.04 and compile it, and then make a simbolic link "Yorick" at your home directory.
* Type "~/Yorick/yorickvis" for starting yorick with yorickvis functions.
* Type "~/Yorick/yorickvis  ~/Yorick/idl_demo1.i" for running "idl_demo1.i" directory
* Type "~/Yorick/yorickvisf datafile" for visualization of ongoing simulation.
* Type "help,function_name" at yorick prompt show its example.

Functions in idl.i
-------

* idl: idling function for event driven programing controlled by mouse click.

Functions in idl_demo0.i, idl_demo1.i, idl_demo2.i, idl_demo3.i
-------

* idl_demo0: demo for idl
* idl_demo1: demo for idl
* idl_demo2: demo for idl
* idl_demo3: demo for idl

Functions in graph2d.i
-------

* win2: changes window properties such as size, shape and ticks.

* plgl:"plg" function with graph-legends on window.
* plgf:"plg" function filling under the line with specified color.
* plh: plots histogram from frequency data (with graph-legends as line segments).
* plhf: plots histogram from frequency data (with graph-legends as filled boxes).
* pler: plots error bars.

* fp: plots line of specified function y=f(x).
* fp2: plots a surface of 2D function f(x,y) with plwf.
* fpfc: plots filled contour of 2D function z=f(x,y) with plfc.
* cbfc: plots color bar with tickes for plfc function.
* make_hist: make frequency data from raw data, and plot histogram.
* hist2d: aggregates 2D dataset.
* make_xy_plf: make x-y coordinates for plf.
* xyt: plots titles for x- and y-axes with a slightly finer adjustment than xytitles.

* pal: sets color palette.
* cb: draws a color bar to the right of the plot by "plf" function.
* fmal: resets plottings and graph-legend.
* outpng: saves plottings in png format.
* hcpng: saves plottings in png and/or eps format.

* scalexy: scales plotting in horizontal and vertical directions separately.
* scalex: scales plotting in horizontal direction.
* scaley: scales plotting in vertical direction.
* scale: scales plotting in horizontal and vertical directions by the same ratio.
* I:short expression of "include" function.


Functions in graph3d.i
-------

* ro: rotates 3D-plottings with mouse.

* win3: initializes the window to 3D-mode
* changech: changes properties of cage3.
* read_poly: reads polygon data for pl3surf and pl3tree to a text file.
* write_poly: writes polygon data for pl3surf and pl3tree to a text file.

* lim3: a short-name version of limit3.
* lim32xz: adjusts the rotation and scale of plottings so that the scales of x- and z-axes in "win2" coincide with those of x- and z-axes on the 3D-cage.
* lim32xy: adjusts the rotation and scale of plottings so that the scales of x- and y-axes in "win2" coincide with those of x- and y-axes on the 3D-cage.
* lim32yz: adjusts the rotation and scale of plottings so that the scales of y- and z-axes in "win2" coincide with those of y- and z-axes on the 3D-cage.

* spal: prepares two different palettes (mainly for pl3tree).
* make_mypal : gets palette data.
* make_xyz: transform data for plwf into data for pl3surf and pl3tree.
* get_light:calculates brightness of each polygon consisting surface in the same manner with pl3surf.


Functions in tube3.i
-------

* tube3: 3D-tube plot.

Functions in lcontour3.i
-------

* lcontour3: cuts polygon into parts between contours for values defined on the surface of the polygon.
* lcon3: emurates smooth shading.
* plfc3: plots 3D filled contours for function z=f(y,x)

