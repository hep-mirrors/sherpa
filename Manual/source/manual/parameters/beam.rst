 .. _Beam Parameters:

***************
Beam parameters
***************

.. index:: BEAMS
.. index:: BEAM_ENERGIES

Mandatory settings to set up the colliding particle beams are

* The initial beam particles specified through :option:`BEAMS`, given
  by their PDG particle number. For (anti)protons and (positrons)
  electrons, for example, these are given by :math:`(-)2212` or
  :math:`(-)11`, respectively.  The code for photons is 22. If you
  provide a single particle number, both beams will consist of that
  particle type. If the beams consist of different particles, a list
  of two values have to be provided.

* The energies of both incoming beams are defined through
  :option:`BEAM_ENERGIES`, given in units of GeV. Again, single values
  apply to both beams, whereas a list of two values have to be given
  when the two beams do not have the same energy.


Examples would be

.. code-block:: yaml

   # LHC
   BEAMS: 2212
   BEAM_ENERGIES: 7000

   # HERA
   BEAMS: [-11, 2212]
   BEAM_ENERGIES: [27.5, 820]

More options related to beamstrahlung and intrinsic transverse momentum can
be found in the following subsections.

.. contents::
   :local:

.. _Beam Spectra:

Beam Spectra
============

.. index:: BEAM_SPECTRA
.. index:: SPECTRUM_FILES
.. index:: BEAM_SMIN
.. index:: BEAM_SMAX

If desired, you can also specify spectra for beamstrahlung through
``BEAM_SPECTRA``. The possible values are

:option:`Monochromatic`
  The beam energy is unaltered and the beam
  particles remain unchanged.  That is the default and corresponds to
  ordinary hadron-hadron or lepton-lepton collisions.

:option:`Laser_Backscattering`
  This can be used to describe the
  backscattering of a laser beam off initial leptons. The energy
  distribution of the emerging photon beams is modelled by the CompAZ
  parameterisation, see :cite:`Zarnecki2002qr`.  Note that this
  parameterisation is valid only for the proposed TESLA photon
  collider, as various assumptions about the laser parameters and the
  initial lepton beam energy have been made. See details below.

:option:`Simple_Compton`
  This corresponds to a simple light backscattering
  off the initial lepton beam and produces initial-state
  photons with a corresponding energy spectrum.  See details below.

:option:`EPA`
  This enables the equivalent photon approximation for colliding
  protons, see :cite:`Archibald2008aa`. The resulting beam particles
  are photons that follow a dipole form factor parameterisation,
  cf. :cite:`Budnev1974de`.  The authors would like to
  thank T. Pierzchala for his help in implementing and testing the
  corresponding code. See details below.

:option:`Pomeron`
  This enables the Proton--Pomeron flux for diffractive jet production, see
  details below.

:option:`Reggeon`
  This enables the Proton--Reggeon flux, see details below.

The ``BEAM_SMIN`` and ``BEAM_SMAX`` parameters may be used to specify
the minimum/maximum fraction of cms energy squared after
Beamstrahlung. The reference value is the total centre of mass energy
squared of the collision, *not* the centre of mass energy after
eventual Beamstrahlung.

The parameter can be specified using the internal interpreter, see
:ref:`Interpreter`, e.g. as ``BEAM_SMIN: sqr(20/E_CMS)``.

Laser Backscattering
--------------------

.. index:: E_LASER
.. index:: P_LASER
.. index:: LASER_MODE
.. index:: LASER_ANGLES
.. index:: LASER_NONLINEARITY

The energy distribution of the photon beams is modelled by the CompAZ
parameterisation, see :cite:`Zarnecki2002qr`, with various assumptions
valid only for the proposed TESLA photon collider. The laser energies
can be set by ``E_LASER``. ``P_LASER`` sets their polarisations,
defaulting to ``0.``.  Both settings can either be set to a single
value, applying to both beams, or to a list of two values, one for
each beam.  The ``LASER_MODE`` takes the values ``-1``, ``0``, and
``1``, defaulting to ``0``.  ``LASER_ANGLES`` and
``LASER_NONLINEARITY`` can be set to ``true`` or to ``false``
(default).

Simple Compton
--------------

This corresponds to a simple light backscattering off the initial
lepton beam and produces initial-state photons with a corresponding
energy spectrum.  It is a special case of the above Laser
Backscattering with ``LASER_MODE: -1``.

EPA
---

.. index:: EPA:Q2Max
.. index:: EPA:ThetaMax
.. index:: EPA:Use_old_WW
.. index:: EPA:PTMin
.. index:: EPA:Form_Factor
.. index:: EPA:AlphaQED

The equivalent photon approximation, cf. :cite:`Archibald2008aa`,
:cite:`Budnev1974de`, has a few free parameters, listed below.
Each of these parameters has to be set in the subsetting ``EPA``, like so

.. code-block:: yaml

   EPA:
     Q2Max: 4.5

The usual rules for yaml structure apply, c.f. :ref:`Input structure`.

:option:`Q2Max`
  Parameter of the EPA spectra of the two beams, defaults to ``3.`` in
  units of GeV squared. For the electron, the maximum virtuality is taken
  to be the minimum of this value and the kinematical limit, given by

  .. math::

    Q^2_{max,kin} = \frac{(m_e x)^2}{1-x} + E_e^2 (1-x) \theta^2_{max}

  with :math:`m_e` the electron mass, :math:`E_e` the electron energy,
  :math:`x` the energy fraction that the photon carries and
  :math:`\theta_{max}` the maximum electron deflection angle, see below.

:option:`ThetaMax`
  Parameter of the EPA spectrum of an electron beam, cf. :cite:`Frixione:1993yw`.
  Describes the maximum angle of the electron deflection, which
  translates to the maximum virtuality in the photon spectrum. It defaults to ``0.3``.

:option:`Use_old_WW`
  In Sherpa version 3, a more accurate Weizsäcker-Williams weight for electron beams is used, as described in
  :cite:`Schuler:1996qr` and :cite:`Frixione:1993yw`. By default, Sherpa uses this improved version of the formula,
  if you would like to use the previous version, set this switch to ``true``.

:option:`PTMin`
  Infrared regulator to the EPA beam spectra. Given in GeV, the value
  must be between ``0.`` and ``1.`` for EPA approximation to hold.
  Defaults to ``0.``, i.e. the spectrum has to be regulated by cuts on
  the observable, cf :ref:`Selectors`.

:option:`Form_Factor`
  Form factor model to be used on the beams. The options are ``0``
  (pointlike), ``1`` (homogeneously charged sphere, ``2`` (gaussian
  shaped nucleus), and ``3`` (homogeneously charged sphere, smoothed
  at low and high x). Applicable only to heavy ion beams.  Defaults to
  ``0``.

:option:`AlphaQED`
  Value of alphaQED to be used in the EPA. Defaults to ``0.0072992701``.

``Q2Max``, ``PTMin``, ``Form_Factor``, ``XMin`` can either be set to
single values that are then applied to both beams, or to a list of two
values, for the respective beams.

Pomeron
-------

The Pomeron flux is implemented as used in :cite`H1:2006zyl` :cite:`Goharipour:2018yov` :cite:`H1:2006uea` and, integrating out the momentum transfer, is given by

.. math::

    f_{\mathbb{P}}(x) = \int^0_{-t_\mathrm{max}} A_\mathbb{P} \frac{e^{B_\mathbb{P} t}}{{x}_\mathbb{P}^{2 \alpha_\mathbb{P}\left(t\right) -1}}
    = A_\mathbb{P} x^{1 - 2 \alpha\left(0\right)}
    \frac{1-\mathrm{e}^{-B_\mathbb{P} t_\mathrm{max}} x^{2 \alpha^\prime t_\mathrm{max}}}
         {B_\mathbb{P} - 2 \alpha^\prime \mathrm{log}(x)}

where :math:`t` is the squared transferred four-momentum and :math:`\alpha` is assumed to be
linear, :math:`\alpha_\mathbb{P}\left(t\right) = \alpha\left(0\right) + \alpha^\prime t`. The default values are set
to the ones obtained in Fit A in :cite:`H1:2006zyl` and can each be changed like so:

.. code-block:: yaml

    Pomeron:
      tMax: 1.
      xMax: 1.
      xMin: 0.
      B: 5.5
      Alpha_intercept: 1.111
      Alpha_slope: 0.06

where ``Alpha_intercept`` and ``Alpha_slope`` are :math:`\alpha\left(0\right)` and :math:`\alpha^\prime`, respectively.
Please note that ``tMax`` is the absolute value, i.e. a positive number.
``xMax`` denotes the fraction of the proton momentum taken by the Pomeron.

Other fluxes can be implemented upon request.

Reggeon
-------

The Reggeon flux, defined in complete analogy to the Pomeron flux above.
Default values taken from :cite:`H1:2006zyl`, set to:

.. code-block:: yaml

    Reggeon:
      tMax: 1.
      xMax: 1.
      xMin: 0.
      B: 1.6
      Alpha_intercept: 0.5
      Alpha_slope: 0.3
      n: 1.4e-3

The parameter ``n`` is the relative normalization of the Reggeon flux with
respect to the Pomeron flux.

.. _Beam Polarization:

Beam Polarization
=============================

Sherpa can also provide cross-sections for polarized beams.
These calculations can only be provided using the  ``AMEGIC`` ME generator.
The value for the beam polarization can be given as a percentage e.g. 80 or in decimal form e.g. 0.8 .
The flavour of :option:`BEAM_1/BEAM_2` follows the definition given to  :option:`BEAMS`.

.. code-block:: yaml

   POLARIZATION:
     BEAM_1: 0.8
     BEAM_2: -0.3

