OPTICAL FIBRE TOOLBOX 
Created in Dieter Meschede research group
by K. Karapetyan et al. (2008--2011)

Contact: kotya.karapetyan@gmail.com

OFT is licensed under the BSD and CC-BY licenses.

Optical Fibre Toolbox (OFT) provides functions for fast automatic calculation of guided modes in simple optical fibres. Developed with tapered microfibres (aka nanofibres) in mind. Exact solutions for weak and strong guidance cases are provided. Material dispersion is taken into account.

Main functionality: 
- Find the guided modes.
- Calculate the effective refractive index of each mode for the given diameter and wavelength or versus variable diameter or wavelength (modal dispersion).
- Calculate the electric and magnetic fields of the modes (only two-layer modes).
- Find phase-matching points for harmonic generation in fibres.

A number of additional useful functions is provided.

For installation instruction, see install.txt.
For examples of use, see the demo/ directory.
For the full list of functions, see Contents.m.
For the list of authors, see authors.txt.

REVISIONS
==== Version 2.1 (2011-12-08) ====
Bug fix:
- Erdogan modes didn't work for neff vs. wavelength --- fixed.
- Several small stability improvements.

New functionality:
- skipmodes parameter can be specified in mode task. The listed modes will be found by buildModes but not traced
- New materials in refrIndex.

==== Version 2.0 (2011-10-31) ====
This version includes several major changes:
1. First of all, starting with this release OFT supports three-layer structures (double-clad or core-cladding-surrounding). See demo/tutorial3ls.m. 
NOTE: The field of three-layer modes cannot yet be calculated.
2. Material refractive index (function refrIndex.m) can now be calculated taking temperature into account.
3. An attempt to start documenting the toolbox is made. Type 'doc contents' to see a full list of functions or 'doc functionname' to see documentation for a specific function.
4. Also the first tutorial created --- demo/tutorial3ls.m.
5. Starting with this release, the latest stable version of the code will be available at http://code.google.com/p/optical-fibre-toolbox/. This is mainly done to speed up the process of publishing debug versions and to keep track of the history of releases.

Other changes:
- General code clean up.
- V-parameter calculation added.
- White sector in displayField2.m suppressed in a more correct way.
- Many other small improvements and corrections.

==== 2010-06-12 ====
Several bugs in field display corrected.

==== 2010-06-02 ====
The first public release of OFT (version 1.0). Only two layer structures are supported.
