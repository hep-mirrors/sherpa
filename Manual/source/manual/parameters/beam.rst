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
  parametrization, see :cite:`Zarnecki2002qr`.  Note that this
  parametrization is valid only for the proposed TESLA photon
  collider, as various assumptions about the laser parameters and the
  initial lepton beam energy have been made. See details below.

:option:`Simple_Compton`
  This corresponds to a simple light backscattering
  off the initial lepton beam and produces initial-state
  photons with a corresponding energy spectrum.  See details below.

:option:`EPA`
  This enables the equivalent photon approximation for colliding
  protons, see :cite:`Archibald2008aa`. The resulting beam particles
  are photons that follow a dipole form factor parametrization,
  cf. :cite:`Budnev1974de`.  The authors would like to
  thank T. Pierzchala for his help in implementing and testing the
  corresponding code. See details below.

:option:`Spectrum_Reader`
  A user defined spectrum is used to describe the energy spectrum
  of the assumed new beam particles. The name of the corresponding
  spectrum file needs to be given through the keywords
  ``SPECTRUM_FILES``.

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
parametrisation, see :cite:`Zarnecki2002qr`, with various assumptions
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
.. index:: EPA:XMin
.. index:: EPA:Use_old_WW
.. index:: EPA:PTMin
.. index:: EPA:Form_Factor
.. index:: EPA:AlphaQED

The equivalent photon approximation, cf. :cite:`Archibald2008aa`,
:cite:`Budnev1974de`, has a few free parameters, listed below.
Each of these parameters has to be set in the subsetting ``EPA``, like so

.. code-block:: yaml

   EPA:
     XMin: 0.01

The usual rules for yaml structure apply, c.f. :ref:`Input structure`.

:option:`Q2Max`
  Parameter of the EPA spectra of the two beams, defaults to ``3.`` in
  units of GeV squared. For the electron, the maximum virtuality is taken
  to be the minimum of this value and the kinematical limit, given by

  .. math::

    Q^2_{max,kin} = \frac{(m_e x)^2}{1-x} + E_e^2 (1-x) \theta^2_{max}

  with :math:`m_e` is the electron mass, :math:`E_e` the electron energy,
  :math:`x` the energy fraction that the photon carries and
  :math:`\theta_{max}` the maximum electron deflection angle, see below.

:option:`ThetaMax`
  Parameter of the EPA spectrum of an electron beam, c.f. :cite:`Frixione:1993yw`.
  Describes the maximum angle of the electron deflection, which
  translates to the maximum virtuality in the photon spectrum. It defaults to ``0.3``.

:option:`XMin`
  Restricts the phase space by imposing a minimum energy fraction that the photon must have with respect to the beam energy.
  Its default value is ``0``.

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

.. _Intrinsic Transverse Momentum:

Intrinsic Transverse Momentum
=============================

.. index:: INTRINSIC_KPERP
.. index:: K_PERP_SCHEME
.. index:: K_PERP_MEAN
.. index:: K_PERP_SIGMA
.. index:: K_PERP_EXP
.. index:: K_PERP_EREF
.. index:: BEAM_REMNANTS

The intrinsic transverse momentum of the colliding particles can be
set by subsettings of the ``INTRINSIC_KPERP`` setting:

.. code-block:: yaml

   INTRINSIC_KPERP:
     Parameter_1: <value_1>
     Parameter_2: <value_2>
   ...

The possible parameters and their defaults are

:option:`ENABLED (default: true)`
  Setting this to ``false`` disables the intrinsic transverse momentum
  altogether.

:option:`SCHEME (default for protons: 0)`
  This parameter specifies the scheme to calculate the intrinsic
  transverse momentum of the beams in case of hadronic beams, such as
  protons.

:option:`MEAN (default: 1.1)`
  This parameter specifies the mean intrinsic transverse momentum in GeV
  for the beams in case of hadronic beams, such as protons.

  If two values are provided, the intrinsic momenta means of the two
  beams are set to these two values, respectively.

:option:`SIGMA (default: 0.85)`
  This parameter specifies the width of the Gaussian distribution in GeV
  of the intrinsic transverse momenta for the beams in case of
  hadronic beams, such as protons.

  If two values are provided, the intrinsic momenta widths of the two
  beams are set to these two values, respectively.

:option:`EXP (default: 0.55)`
  This parameter specifies the energy extrapolation exponent of the
  width of the Gaussian distribution of intrinsic transverse momentum
  for the beams in case of hadronic beams, such as protons.

  If two values are provided, the exponents for each of the two beams
  are set to these two values, respectively.

:option:`EREF (default: 7000)`
  This parameter specifies the reference scale in GeV in the energy
  extrapolation of the width of the Gaussian distribution of intrinsic
  transverse momentum for the beams in case of hadronic beams, such as
  protons.

  If two values are provided, the reference scales for each of the two
  beams are set to these two values, respectively.

If the option :option:`BEAM_REMNANTS: false` is specified, pure
parton-level events are simulated, i.e. no beam remnants are
generated. Accordingly, partons entering the hard scattering process
do not acquire primordial transverse momentum.
