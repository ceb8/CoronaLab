

***************
Getting Started
***************


``CoronaLab`` is a package for generating synthetic data products from a model stellar corona, it does not  require  the  user to use a particular method for modelling the corona, as long as a set of open and closed magnetic field lines  is  produced.

        
====================
Installing CoronaLab
====================


Using pip
---------

The easiest way to install Coronalab is using pip::

    pip install corona_lab


From source
-----------

To install the bleeding edge version from github without downloading,
run the following command::

  pip git+https://github.com/St-andrews-cool-stars/CoronaLab.git

The latest development version of Coronalab can be cloned from github
and then install with the following::

    git clone https://github.com/St-andrews-cool-stars/CoronaLab.git
    cd CoronaLab
    pip install .

For active developement intall in develop mode

    pip install -e .


=================================
Initializing a ModelCorona object
=================================


The main ``CoronaLab`` class is  `~corona_lab.corona.ModelCorona`, which holds the model corona from which synthetic observations will be created. This object  can  be initialized in two ways: by supplying a valid Table with all required columns and meta data, or through the `~corona_lab.build_corona.FieldlineProcessor` class where various attributes  of the model are set and then magnetic field  measurements  along open and closed field lines are supplied so  that  all other required properties can be calculated.



Using `~corona_lab.corona.ModelCorona.from_field_lines`
-------------------------------------------------------

The minimum required columns and metadata for a valid `~corona_lab.corona.ModelCorona` object are listed below. The metadata can be supplied directly as arguments to `~corona_lab.corona.ModelCorona.from_field_lines`, or as metadata on the input table. All values should be `~astropy.units.Quantity` objects with appropriate units.
    
**Required Columns**

- ``Radius``, ``Theta``, ``Phi`` Position in spherical coordinates (rotation axis is assumed to be at θ = 90°).
- ``dV`` Volume of field line segment
- ``dA_c`` Cross-sectional area  of field line segment
- ``Bmag`` Magnetic field strength of field line segment
- ``ndens`` Number density of field line segment
- ``Temperature`` Temperature of field line segment
- ``prom`` Is the segment part of a prominence (boolean)
- ``Mprom`` Prominence mass of segment (0 if segment is not part of a prominence)
- ``line_num`` Field line ID number positive for closed corona, negative for wind

**Required Metadata**

- ``Radius`` Stellar radius
- ``Mass`` Stellar mass
- ``Source Surface Radius`` Model outer radius
- ``Period`` Stellar rotation period
- ``Distance`` Distance from observer to star (optional)



Using `~corona_lab.build_corona.FieldlineProcessor`
---------------------------------------------------

The required inputs for  building a `~corona_lab.corona.ModelCorona` object from a set of  magnetic field lines are given  below (note that all values should be `~astropy.units.Quantity` objects with appropriate units).

**Open and closed field lines, as individual tables with columns:**
    
- ``radius``, ``theta``, ``phi`` Position in spherical coordinates (rotation axis is at θ = 90°).
- ``ds`` Segment length for each field line cell.
- ``Brad``, ``Btheta``, ``Bphi`` Magnetic field components.
- ``Bmag``  Magnetic field strength ( √(B\ :sub:`r`\ ² + B\ :sub:`θ`\ ² + B\ :sub:`φ`\ ²) )
  

**Stellar Characteristics**
    
- ``mass`` Stellar mass
- ``radius`` Stellar radius
- ``period`` Stellar rotation period
- ``mean_ptc_mass`` Mean particle mass in the stellar corona (optional, defaults to pur hydrogen)
- ``distance`` Distance from observer to star (optional)

    
**Model Parameters**
    
- ``T_cor`` Corona temperature (assumed isothermal excepting prominences)
- ``T_prom`` Prominence temperature (optional, default is 85000 K)
- ``kappa_power`` Pressure scaling factor, relates base pressure to magnetic field strength as p\ :sub:`0` = 10\ :sup:-`kappa_power` B\ :sub:`0`\ ²
- ``rss`` Source surface, considered the outer radius of the model.
- ``dtheta``, ``dphi`` Model grid resolution in θ and  φ (determines cross-sectional area of field line flux tube)

:ref:`build-corona` works through an example usage of this functionality, and information on the specific algorithms used can be found :doc:`here <algorithms>`.


===========
Basic Usage
===========

Here we work through a basic example that shocases all of the available functionality.

    
.. include:: basic_example.rst


===========
Limitations
===========

``CoronaLab`` is still under active development, and in its current version is intended primarily 
for the study of M Dwarfs and solar-type stars in the radio regime. Use beyond these domains should
be undertaken with caution and with the knowledge that software modification may be required.

- **Stellar wind is not modelled.** Users interested in assessing the contribution of wind or anticipating its influence on their results must implement an appropriate treatment.
- **Scattering effects in free-free image synthesis are neglected.** This assumption is based on the large mean free paths in the coronae we modelled, which rendered scattering negligible. However, this condition will not hold universally.
- **Emission mechanisms are limited to free-free and ECM processes.** For many stellar types and observational wavelengths, additional mechanisms such as gyroresonance, gyrosynchrotron emission, and plasma radiation may contribute significantly and are not currently included.
-  **Wavelength coverage is restricted,** The software assumes that all synthetic observational products lie within the radio regime. Users seeking to model emission at other wavelengths will need to extend the relevant portions of the code accordingly.


