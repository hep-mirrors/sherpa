.. _Beam Remnant Parameters:

*********************
Beam Remnants
*********************

Details for the handling of the beam remnants in Sherpa will be described
in our forthcoming publication.

Broadly speaking, the beam remnants include the parametrization of the
form factors for hadrons or the hadronic components of photons and the
treatment of the beam break-up, most importantly the intrinsic
transverse momentum distribution of the partons and how the recoils
are distributed.

The following parameters are used to steer the beam remnant handling:

.. contents::
   :local:


.. _BEAM_REMNANTS:

BEAM_REMNANTS
==========

.. index:: BEAM_REMNANTS

Specifies whether beam remnants are taken into account, with possible
values 'On' and 'Off'.

.. _Remnants:

REMNANTS
======

.. index:: REMNANTS:KT_FORM
.. index:: REMNANTS:KT_RECOIL
.. index:: REMNANTS:SHOWER_INITIATOR_MEAN
.. index:: REMNANTS:SHOWER_INITIATOR_SIGMA
.. index:: REMNANTS:SHOWER_INITIATOR_Q2
.. index:: REMNANTS:SHOWER_INITIATOR_KTMAX
.. index:: REMNANTS:SHOWER_INITIATOR_KTEXPO
.. index:: REMNANTS:BEAM_SPECTATOR_MEAN
.. index:: REMNANTS:BEAM_SPECTATOR_SIGMA
.. index:: REMNANTS:BEAM_SPECTATOR_Q2
.. index:: REMNANTS:BEAM_SPECTATOR_KTMAX
.. index:: REMNANTS:BEAM_SPECTATOR_KTEXPO
.. index:: REMNANTS:REFERENCE_ENERGY
.. index:: REMNANTS:ENERGY_SCALING_EXPO
.. index:: REMNANTS:MATTER_FRACTION1
.. index:: REMNANTS:MATTER_RADIUS1
.. index:: REMNANTS:MATTER_RADIUS2
.. index:: REMNANTS:MATTER_FORM


Sherpa organises the remnant handling by particle, with the PDG code as
tag-line.  

.. code-block:: yaml

   REMNANTS:
     2212:
       KT_Form: Gauss_limited

The usual rules for yaml structure apply, c.f. :ref:`Input structure`.
Longitudinal momenta for sea partons in hadrons are distributed according
to a probablity distribution in their light-cone momentum :math:`x` given by
:math:`P(x)=x^{-1.5}`. If there are two valence partons left in the beam remnant
after the shower initiators have been treated, the first of the two (usually the
quark) will have a longitudinal momentum with :math:`P(x)=\exp(-1/x)`, while
the last remaining valence parton (usually the di-quark for nucleons) carries
the remaining longitudinal momentum.

For the intrinsic transverse momentum, Sherpa differentiates between the
transverse momentum for shower initiators (``SHOWER_INITIATOR_MEAN`` etc.)
and for beam spectators (``BEAM_SPECTATOR_MEAN`` etc.), and it offers
different strategies to compensate the transverse momentum between the
two sets of partons per beam, see below (``KT_RECOIL``).

:option:`KT_FORM (default: Gauss_Limited)`
  This parameter specifies the scheme to calculate the intrinsic transverse
  momentum of partons within beams.  Available options are:
  
  * ``Gauss``: a simple Gaussian with mean and width;
  * ``Dipole``: a dipole form parametrised by :math:`Q^2`;
  * ``Gauss_Limited``, ``dipole_Limited``: as above but further modifed by a polynomial function of the form :math:`1-(k_{T}/k_{T,\rm{max}})^\eta`, where :math:`k_{T,\rm{max}}` and :math:`\eta` are given by the ``KTMAX`` and ``KTEXPO`` tags;
  * ``None``: no intrinsic transverse momentum is assigned.

:option:`KT_RECOIL (default: Beam_vs_Shower)`
  Transverse momenta for all partons inside the beam are generated
  independently from each other according to the form and parametrization
  specified for them in ``KT_FORM`` and ``SHOWER_INITIATOR_MEAN`` etc., or
  ``BEAM_SPECTATOR_MEAN`` etc..  This will lead to a net residual transverse
  momentum of partons that needs to be compensated within the beams, to
  guarantee that the remnants do not create a total beam transverse
  momentum.  Sherpa has implemented two strategies to achieve this:
  
  * ``Democratic``: the overall residual transverse momentum is distributed over all partons in the beam according to their energies.  
  * ``Beam_vs_Shower``: the residual transverse momentum of all spectators is distributed over the shower initiators according to their energies and vice versa.
	
:option:`SHOWER_INITIATOR_MEAN (default for nucleons: 1.0)`
  This parameter specifies the mean in GeV for the intrinsic
  transverse momentum in case of a limited or unlimited
  Gaussian distribution.

:option:`BEAM_SPECTATOR_MEAN   (default for nucleons: 0.0)`
  Same as for ``SHOWER_INITIATOR_MEAN``.

:option:`SHOWER_INITIATOR_SIGMA (default for nucleons: 1.1)`
  This parameter specifies the sigma in GeV for the intrinsic
  transverse momentum in case of a limited or unlimited
  Gaussian distribution.

:option:`BEAM_SPECTATOR_SIGMA   (default for nucleons: 0.25)`
  Same as for ``SHOWER_INITIATOR_SIGMA``.

:option:`SHOWER_INITIATOR_Q2 (default for nucleons: 1.1)`
  This parameter specifies the :math:`Q^2` in :math:`{\rm GeV}^2`
  of the limited or unlimited dipole distribution for the
  intrinsic transverse momentum.

:option:`BEAM_SPECTATOR_Q2   (default for nucleons: 0.25)`
  Same as for ``SHOWER_INITIATOR_Q2``.

:option:`SHOWER_INITIATOR_KTMAX (default for nucleons: 2.7)`
  This parameter specifies the :math:`k_{T,\rm{max}}` in
  :math:`{\rm GeV}` of the limited dipole or Gaussian distributions
  for the intrinsic transverse momentum. 

:option:`BEAM_SPECTATOR_KTMAX   (default for nucleons: 1.0)`
  Same as for ``SHOWER_INITIATOR_KTMAX``.

:option:`SHOWER_INITIATOR_KTEXPO (default for nucleons: 5.12)`
  This parameter specifies the :math:`\eta` in the equation above
  that limits the intrinsic transverse momentum distribution. 

:option:`BEAM_SPECTATOR_KTEXPO   (default for nucleons: 5.0)`
  Same as for ``SHOWER_INITIATOR_KTEXPO``.

:option:`REFERENCE_ENERGY (default: 7000)`
  This parameter specifies the reference scale in GeV in the energy
  extrapolation of the mean and width of the Gaussian distribution
  and of the :math:`Q^2` of the dipole distribution of intrinsic
  transverse momentum, and of the maximally allowed :math:`k_T`
  in the case of limited distributions.

:option:`ENERGY_SCALING_EXPO (default: 0.08)`
  This parameter specifies the energy extrapolation exponent.

:option:`MATTER_FORM (default: Single_Gaussian)`
  ``Double_Gaussian`` can be used to model the overlap between
  the colliding particles.  ``None`` switches this off.

:option:`MATTER_RADIUS1 (default for nucleons: 0.86, for mesons/photons: 0.75)`
  The radius of the (inner) Gaussian in fm. If used with the
  double-Gaussian matter form, this value must be smaller than ``MATTER_RADIUS2``.

:option:`MATTER_FRACTION1`
  Only to be used for double-Gaussian matter form, where it will control the
  distribution of matter over the two Gaussians. It assumes that a fraction
  :math:`f^2` is distributed by the inner Gaussian :math:`r_1`, another fraction
  :math:`(1-f)^2` is distributed by the outer Gaussian :math:`r_2`,
  and the remaining fraction :math:`2f(1-f)` is distributed by the combined radius
  :math:`r_\text{tot} = \sqrt{\frac{r_1^2+r_2^2}{2}}`. Defaults to ``0.5``.

:option:`MATTER_RADIUS2`
    Defaults to ``1.0``. It is only used for the case of a double-Gaussian
    overlap, see below.


If the option :option:`BEAM_REMNANTS: false` is specified at top level, pure
parton-level events are simulated, i.e. no beam remnants are
generated. Accordingly, partons entering the hard scattering process
do not acquire primordial transverse momentum.

