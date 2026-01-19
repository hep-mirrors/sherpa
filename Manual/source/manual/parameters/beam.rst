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
  of two values has to be provided.

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

.. _Laser Backscattering:

Laser Backscattering
====================

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

.. _Simple Compton:

Simple Compton
==============

This corresponds to a simple light backscattering off the initial
lepton beam and produces initial-state photons with a corresponding
energy spectrum.  It is a special case of the above Laser
Backscattering with ``LASER_MODE: -1``.

.. _EPA:

EPA
===

.. index:: EPA

The Equivalent Photon Approximation (EPA) converts incoming charged particles
(leptons, protons, ions) into photons plus remnants. This is particularly
suited for ultra-peripheral collisions at the LHC (pp, pPb, PbPb).
The implementation is based on Budnev et al. :cite:`Budnev1974de`
and includes extensive form factor models for hadrons and nuclei.

EPA is enabled through :option:`BEAMSPECTRA: EPA`. Parameters are set in an
:option:`EPA` block, with single values applying to both beams or lists
``[beam1, beam2]`` for asymmetric setups.

.. code-block:: yaml

    BEAMS: 2212
    BEAM_ENERGIES: 6500
    BEAM_SPECTRA: EPA
    PDF_SET: None

    EPA:
      Form_Factor: Dipole

.. code-block:: yaml

    # LHC PbPb example
    BEAMS: 1000822080
    BEAM_ENERGIES: 2510
    BEAM_SPECTRA: EPA
    PDF_SET: None

    EPA:
      FormFactor: Woods-Saxon
      WoodsSaxon_R: 6.49
      WoodsSaxon_d: 0.54

.. contents::
   :local:

.. _EPAPhysics-details:

Physics details
---------------

The electromagnetic field of a fast moving
charged particle can be treated as a flux of quasi-real photons.
For simple expression, e.g. for leptons, the corresponding flux can be computed analytically
and we factorise the cross-section the following way:

.. math::

    \sigma_{A B \to X}
    = \int \mathrm{d}x_1\, f_{\gamma/A}(x_1)
    \int \mathrm{d}x_2\, f_{\gamma/B}(x_2)
    \,\mathrm{d}\hat{\sigma}_{\gamma\gamma\to X}(x_1, x_2, \mu^2) \\
    \mbox{ with }\ f_{\gamma/A,B}(x) = \frac{\alpha_\mathrm{em}}{\pi} N(x)

However, for nucleons and ions, the spatial extent of the source needs to be taken
into account. In these cases, we need to Fourier transform the electromagnetic
form factor with the following formula :cite:`Vidovic:1992ik`:

.. math::

   f_{\gamma}(x,\mathbf{r})
   =& \frac{Z^2 \alpha_{\text{em}}}{\pi^2}\,\frac{1}{x}\,
     \left|
       \int_0^\infty \mathrm{d}k_\perp\,
       \frac{k_\perp^2}{Q^2}\,
       F\big(k_\perp^2 + (y m_p)^2\big)\,
       J_1(|\mathbf{r}|\,k_\perp)
     \right|^2 \\
   \implies f_{\gamma}(x, b)
   =& \frac{2 Z^2 \alpha_{\text{em}}}{\pi}\,\frac{1}{x}\,
     \left|
       \int_0^\infty \mathrm{d}k_\perp\,
       \frac{k_\perp^2}{Q^2}\,
       F\big(Q^2\big)\,
       J_1(|\mathbf{r}|\,k_\perp)
     \right|^2,

where

* :math:`y` is the photon longitudinal momentum fraction,
* :math:`\mathbf{r}` is the distance vector from the emitting beam particle,
* :math:`m_p` is the proton mass or the nucleon mass,
* :math:`J_1` is the Bessel function of the first kind of order one,
* :math:`b = |\mathbf{r}|` is the transverse distance from the emitting beam particle,
* :math:`Q^2` is depending on both :math:`x` and :math:`k_\perp` with the relation :math:`Q^2 = \frac{k_\perp^2 + m_p^2 x^2}{1 - x}`.

In order to accelerate the computation, this convolution is precomputed during
initialisation for a 2-dimensional grid in :math:`x` and :math:`b`.
Settings to control this grid can be found below.

Thus, we arrive at the following factorisation formula:

.. math::

    \sigma_{A B \to X}
    = \int \mathrm{d}x_1\,\int \mathrm{d}b_1\, f_{\gamma/A}(x_1, b_1)
    \int \mathrm{d}x_2\,\int \mathrm{d}b_1\, f_{\gamma/B}(x_2, b_2)
    \,\mathrm{d}\hat{\sigma}_{\gamma\gamma\to X}(x_1, x_2, \mu^2, b_1, b_2) \\
    \mbox{ with }\ f_{\gamma/A,B}(x, b) = \frac{\alpha_\mathrm{em}}{\pi} N(x, b)

The fluxes :math:`N(x)` and :math:`N(x, b)` will be defined with the settings below.
The impact parameter can then be computed as :math:`b = |\mathbf{r_1} - \mathbf{r_2}|`.

.. _EPA Parameters:

EPA Parameters
--------------

.. index:: EPA:Q2Max
.. index:: EPA:Q2Min
.. index:: EPA:xMin
.. index:: EPA:xMax
.. index:: EPA:Form_Factor
.. index:: EPA:ThetaMax
.. index:: EPA:bMin
.. index:: EPA:bThreshold
.. index:: EPA:bMax
.. index:: EPA:xBins
.. index:: EPA:bBins
.. index:: EPA:AlphaQED
.. index:: EPA:Q02
.. index:: EPA:MagneticMu
.. index:: EPA:WoodsSaxon_R
.. index:: EPA:WoodsSaxon_d
.. index:: EPA:WoodsSaxonApprox_a
.. index:: EPA:OutputSpectra

:option:`Q2Max`
  Maximum photon virtuality :math:`Q^2_\mathrm{max}` in GeV^2.
  Defaults to ``1.0``. The actual maximum virtuality used in the cal

:option:`Q2Min`
  Minimum photon virtuality :math:`Q^2_\mathrm{min}` in GeV^2. Defaults to ``-1.0``,
  i.e. use kinematic limit :math:`Q^2_\mathrm{min} = \frac{m^2 x^2}{1 -x}`.
  If set to a value greater than 0, it is dynamically set to the maximum of
  the user-specified value and the kinematic limit.

:option:`xMin`
  Minimum photon energy fraction :math:`x_\mathrm{min}`. Defaults to ``1e-5``.

:option:`xMax`
  Maximum photon energy fraction :math:`x_\mathrm{max}`. Defaults to ``1.0``.

:option:`Form_Factor`
  Form factor model. Auto-selected by beam type, or set explicitly:

  ``Lepton``
    Point-like flux for leptons :cite:`Frixione:1993yw,Schuler:1996qr`:

    .. math::

        N(x) = \frac{1}{2x} \left[ (1+(1-x)^2) \ln\frac{Q^2_\mathrm{max}}{Q^2_\mathrm{min}}
        - 2m^2 x^2 \left( \frac{1}{Q^2_\mathrm{min}} - \frac{1}{Q^2_\mathrm{max}} \right) \right]

  ``Approx_Lepton``
    Point-like flux for leptons without the correction term for :math:`x \to 1`
    :cite:`Budnev1974de`:

    .. math::

        N(x) = \frac{1+(1-x)^2}{2x} \ln\frac{Q^2_\mathrm{max}}{Q^2_\mathrm{min}}

  ``Proton``
    Proton Sachs form factor :cite:`Budnev1974de`, given by

    .. math::

       N(x)
       = \frac{1-x}{x}
         \left[
           \phi\left(x, \frac{Q^2_\mathrm{max}}{Q_0^2}\right)
           - \phi\left(x, \frac{Q^2_\mathrm{min}}{Q_0^2}\right)
         \right]

    with

    .. math::

       \phi(x,z)
       &= \left(1 + a\,y\right)
          \left[
            -\ln\left(1 + \frac{1}{z}\right)
            + \frac{1}{1+z}
            + \frac{1}{2(1+z)^2}
            + \frac{1}{3(1+z)^3}
          \right] \\
       &\quad
          + \frac{1-b}{4}\,\frac{y}{z}\,\frac{1}{(1+z)^3} \\
       &\quad
          + c\left(1 + \frac{y}{4}\right)
          \left[
            \ln\left(\frac{1+z-b}{1+z}\right)
            + \frac{b}{1+z}
            + \frac{b^2}{2(1+z)^2}
            + \frac{b^3}{3(1+z)^3}
          \right],

    and the auxiliary quantities

    .. math::

       y = \frac{x^2}{1-x}, \
       a = \frac{1+\mu^2}{4} + \frac{4 m_p^2}{Q_0^2}, \
       b = 1 - \frac{4 m_p^2}{Q_0^2}, \
       c = (\mu^2 - 1)\,b^{-4}, \
       z = \frac{Q^2}{Q_0^2}.

    Here :math:`m_p` is the proton mass, :math:`\mu` the proton magnetic moment, set by :option:`MagneticMu`, and :math:`Q_0^2` the dipole scale, set by :option:`Q02`.

  ``Approx_Proton``
    Proton Sachs approximation for :math:`Q^2\to0` :cite:`Budnev1974de`:

    .. math::

        N(x) = \frac{1}{x} (1 - x + \frac{\mu^2 x^2}{2}) \ln\frac{Q^2_\mathrm{max}}{Q^2_\mathrm{min}}
        - \frac{1-x}{x} \left( 1 - \frac{Q^2_\mathrm{min}}{Q^2_\mathrm{max}} \right)

    Again, :math:`\mu` is the proton magnetic moment, set by :option:`MagneticMu`.


  ``Point-like_Integrated``
    Point-like nucleus approximation, with impact-parameter dependence integrated out :cite:`Bertulani:1987tz`:

    .. math::

        N(x) = \frac{2}{x} \left[ \chi K_0(\chi) K_1(\chi) - \frac{\chi^2}{2}
        \left( K_1^2(\chi) - K_0^2(\chi) \right) \right] \
        \mbox{ with }\ \chi = x m_N R \mathrm{max}(1, b_{min})

  The above form factors do **not** depend on the impact parameter.
  Following the formulae above in :ref:`EPAPhysics-details`, we implement form factors :math:`F(Q^2)`
  and perform a numerical Fourier transfrom into impact-parameter space.
  In this procedure, the following form factors are available:

  ``Dipole``
    Dipole form factor :cite:`Budnev1974de` :math:`F(Q^2) = (1+Q^2/Q_0^2)^{-2}`, with :math:`Q_0^2` the dipole scale, set by :option:`Q02`.


  ``Approx_Dipole``
    Approximation :math:`Q^2 \to 0` and hence :math:`F(Q^2) = 1` :cite:`Budnev1974de`.

  ``Gaussian``
    Gaussian form factor :cite:`Budnev1974de` :math:`F(Q^2) = \exp(-\frac{Q^2}{2Q_0^2})`, with :math:`Q_0^2` the dipole scale, set by :option:`Q02`.

  ``HCS``
    Homogeneously charged sphere :cite:`Vidovic:1992ik`:

    .. math::

        F(q = \sqrt{Q^2}) = \frac{3}{ (qR)^3 } \left[ \sin(qR) - qR \cos(qR) \right]

  ``Woods-Saxon``
    Woods-Saxon form factor, derived from the nuclear density:

    .. math::

        \rho(r) = \frac{\rho_0 }{ 1 + \exp\left( \frac{r-R}{d} \right) }

        F(Q^2) = 4\pi \int_0^\infty \mathrm{d}r r^2 \; \rho(r) \frac{\sin(qr)}{qr}

    Here, :math:`R` and :math:`d` are the nuclear radius and the nuclear thickness,
    set by :option:`WoodsSaxon_R` and :option:`WoodsSaxon_d`, respectively.
    The normalisation :math:`\rho_0` is chosen such that :math:`F(0) = 1` is fulfilled.

    As this form factor involves a non-trivial integration itself, we precompute
    and tabulate it, for a faster computation of the Fourier transfrom into
    impact-parameter space.

  ``Approx_Woods-Saxon``
    Approximation of the Woods-Saxon form factor as hard sphere, with radius :math:`R`,
    convoluted with a Yukawa potential with range :math:`a`, giving :cite:`Klein:1999qj`

    .. math::

        F(q = \sqrt{Q^2}) = \frac{4\pi\rho_0}{q^3} \left(\sin(qR)-qR\cos(qR)\right) \left(\frac{1}{1+a^2q^2}\right)


    The radius :math:`R` and range :math:`a` are set by :option:`WoodsSaxon_R` and :option:`WoodsSaxonApprox_a`, respectively.
    The normalisation :math:`\rho_0` is chosen such that :math:`F(0) = 1` is fulfilled.

  ``Point-like``
    Ion electric dipole approximation :cite:`Bertulani:1987tz`

    .. math::

        N(x, b) = 2 b x m_N^2 \left( K_1^2(\chi) + \frac{2 m_N}{E_N} K_0^2(\chi) \right)\
        \mbox{ with }\ \chi = x m_N R \mathrm{max}(1, b_{min})

    Please note that with this form factor, there is no need to numerically
    Fourier-transform, as it has been computed analytically.
    However, because it diverges for :math:`b \to 0`, the minimal impact parameter
    is automatically set equal to the radius (or larger,
    depending on user-input for :math:`b_\mathrm{min}`).

  Defaults are ``Lepton`` for leptons, ``Dipole`` for protons and ``Woods-Saxon`` for nuclei.

:option:`ThetaMax`
  Maximum lepton scattering angle in radians, sets
  :math:`p^2_{T,\mathrm{max}}=E^2\theta^2_\mathrm{max}` used to
  determine :math:`Q^2_\text{max}`. Defaults to ``0.3``.

:option:`bMin`
  Minimum impact parameter **relative to the emitter's radius**. Defaults to ``0.1``.

:option:`bMax`
  Maximum impact parameter **relative to the emitter's radius**. Defaults to ``1e3``.

:option:`bThreshold`
  For :math:`b>b_\mathrm{threshold} R`, use point-like approximation. Defaults to ``10.0``.

:option:`xBins`
  Number of :math:`x`-bins for :math:`N(x,b)` grids (logarithmic scaling). Defaults to ``100``.

:option:`bBins`
  Number of :math:`b`-bins for :math:`N(x,b)` grids (logarithmic scaling). Defaults to ``100``.

:option:`AlphaQED`
  QED fine-structure constant :math:`\alpha_\text{em}` used for the EPA calculation. Defaults to ``1/137.036``.

:option:`Q02`
  Dipole scale :math:`\Lambda^2` in GeV^2 for ``Dipole`` and ``Gaussian``. Defaults to ``0.71``.

:option:`MagneticMu`
  Proton magnetic moment for the ``Proton`` form factor. Defaults to ``2.79``.

:option:`WoodsSaxon_R`
  Woods-Saxon nuclear radius :math:`R` in fm. Defaults to ``6.49``.

:option:`WoodsSaxon_d`
  Woods-Saxon skin depth :math:`d` in fm. Defaults to ``0.54``.

:option:`WoodsSaxonApprox_a`
  Yukawa range :math:`a` in fm in the ``Approx_Woods-Saxon`` form factor. Defaults to ``0.7``.

:option:`OutputSpectra`
  If ``true``, output CSV files with :math:`N(x,b)` and :math:`F(Q^2)` for the photon fluxes.
  Defaults to ``false``.

.. _Pomeron:

Pomeron
=======

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

.. _Reggeon:

Reggeon
=======

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
=================

Sherpa can also provide cross-sections for polarized beams.
These calculations can only be provided using the  ``AMEGIC`` ME generator.
The value for the beam polarization can be given as a percentage e.g. 80 or in decimal form e.g. 0.8 .
The flavour of :option:`BEAM_1/BEAM_2` follows the definition given to  :option:`BEAMS`.

.. code-block:: yaml

   POLARIZATION:
     BEAM_1: 0.8
     BEAM_2: -0.3

.. _BEAM_OVERLAP_REJECTION:

BEAM_OVERLAP_REJECTION
======================

.. index:: BEAM_OVERLAP_REJECTION

This parameter controls the rejection of events based on the impact parameter
:math:`b` between the colliding beams, modeling the survival probability
against hadronic interactions. This is particularly important for
ultra-peripheral collisions, where it implements the survival factor.

The following modes are available:

:samp:`{0}`
  No impact-parameter rejection (default). **Not recommended for UPCs.**

:samp:`{1}`
  Hard-sphere rejection based on nuclear radii. Events with impact parameter
  :math:`b < R_1 + R_2` are rejected, where :math:`R_1` and :math:`R_2` are
  the radii of the two beam particles.

:samp:`{2}`
  Model-dependent hadronic survival probability. The specific model is
  automatically selected based on the beam configuration:

  **Proton-proton collisions:**
    Uses the :math:`\sqrt{s}`-dependent parameterization from
    :cite:`Bertulani:2021umy` (Section 3):

    .. math::

       P_{\text{surv}}(b) = \left[ 1 - \exp\left(-\frac{b^2}{2 b_0(\sqrt{s})}\right) \right]^2

    with

    .. math::

       b_0(\sqrt{s}) = 9.81 + 0.211 \log(\sqrt{s}) + 0.0185 \log^2(\sqrt{s})
       \quad [\text{GeV}^{-2}]

    where :math:`\sqrt{s}` is the nucleon-nucleon center-of-mass energy. Events
    are accepted with probability :math:`P_{\text{surv}}(b)`.

  **Proton-nucleus and nucleus-nucleus collisions:**
    Currently not implemented. Please contact the authors if you require these
    modes for your analysis.

The default value is ``0`` (no rejection).

.. note::

   For UPC processes with impact-parameter-dependent fluxes (see above), setting
   ``BEAM_OVERLAP_REJECTION: 1`` (or `2` for proton-proton)
   is **strongly recommended** to obtain
   physically meaningful cross sections. Without this, the calculated cross
   section will include unphysical contributions from small impact parameters
   where hadronic interactions dominate.
