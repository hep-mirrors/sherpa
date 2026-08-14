.. _Hadronization:

*************
Hadronization
*************

The hadronisation setup covers the fragmentation of partons into
primordial hadrons as well as the decays of unstable hadrons into
stable final states.

.. contents::
   :local:

.. _Fragmentation:

Fragmentation
=============


Fragmentation models
--------------------


.. index:: FRAGMENTATION

The ``FRAGMENTATION`` parameter sets the fragmentation module to be
employed during event generation.

* The default is :option:`Ahadic`, enabling Sherpa's native
  hadronisation model AHADIC++ :cite:`Chahal2022rid`, based on
  the cluster fragmentation model introduced in :cite:`Field1982dg`,
  :cite:`Webber1983if`, :cite:`Gottschalk1986bv`, and :cite:`Marchesini1987cf`.

* The hadronisation can be disabled with the value :option:`None`.

* To evaluate uncertainties stemming from the hadronisation, Sherpa
  also provides an interface to the Lund string fragmentation in
  Pythia 8.3 :cite:`Bierlich:2022pfr` by using the setting
  :option:`Pythia8`.  In this case, the standard Pythia settings
  can be used to steer the behaviour of the Lund string,
  see :cite:`Bierlich:2022pfr`. They are specified in their usual
  form in Pythia in a dedicated settings block. Additionally
  a choice can be made to let Pythia directly handle hadron
  decays via the :option:`DECAYS` setting (separate from the
  :option:`Model` switch mentioned below) and whether Pythias or
  Sherpas default masses and widths should be used through the
  :option:`SHERPA_MASSES` setting. By default the choice of generator
  for the masses and widths setting aligns with the decay setting.

.. code-block:: yaml

   FRAGMENTATION: Pythia8
   PYTHIA8:
     PARAMETERS:
       - StringZ:aLund: 0.68
       - StringZ:bLund: 0.98
         ...
     DECAYS: true
     SHERPA_MASSES: false

Hadron constituents
-------------------

.. index:: M_UP_DOWN
.. index:: M_STRANGE
.. index:: M_CHARM
.. index:: M_BOTTOM
.. index:: M_DIQUARK_OFFSET
.. index:: M_BIND_0
.. index:: M_BIND_1

The constituent masses of the quarks and diquarks are given by

* ``M_UP_DOWN`` (0.3 GeV),

* ``M_STRANGE`` (0.4 GeV),

* ``M_CHARM`` (1.8 GeV), and

* ``M_BOTTOM`` (5.1 GeV).

The diquark masses are composed of the quark
masses and some additional parameters,

with

* ``M_DIQUARK_OFFSET`` (0.3 GeV),

* ``M_BIND_0`` (0.12 GeV), and

* ``M_BIND_1`` (0.5 GeV).

Like all settings related to cluster fragmentation these
are grouped under ``AHADIC``.

.. code-block:: yaml

   AHADIC:
     - M_UP_DOWN: 0.3
       ...
     - M_DIQUARK_OFFSET: 0.3


Hadron multiplets
-----------------

.. index:: MULTI_WEIGHT_R0L0_PSEUDOSCALARS
.. index:: MULTI_WEIGHT_R0L0_VECTORS
.. index:: MULTI_WEIGHT_R0L0_TENSORS2
.. index:: MULTI_WEIGHT_R0L1_SCALARS
.. index:: MULTI_WEIGHT_R0L1_AXIALVECTORS
.. index:: MULTI_WEIGHT_R0L2_VECTORS
.. index:: MULTI_WEIGHT_R0L0_N_1/2
.. index:: MULTI_WEIGHT_R1L0_N_1/2
.. index:: MULTI_WEIGHT_R2L0_N_1/2
.. index:: MULTI_WEIGHT_R1_1L0_N_1/2
.. index:: MULTI_WEIGHT_R0L0_DELTA_3/2
.. index:: SINGLET_MODIFIER
.. index:: SINGLETBARYON_MODIFIER
.. index:: CHARMBARYON_ENHANCEMENT
.. index:: BEAUTYBARYON_ENHANCEMENT
.. index:: CHARMSTRANGE_ENHANCEMENT
.. index:: BEAUTYSTRANGE_ENHANCEMENT
.. index:: BEAUTYCHARM_ENHANCEMENT
.. index:: Mixing_0+
.. index:: Mixing_1-
.. index:: Mixing_2+
.. index:: Mixing_3-
.. index:: Mixing_4+
.. index:: ETA_MODIFIER
.. index:: ETA_PRIME_MODIFIER

For the selection of hadrons emerging in such cluster transitions and decays,
an overlap between the cluster flavour content and the flavour part of the
hadronic wave function is formed.  This may be further modified by production
probabilities, organised by multiplet and given by the parameters

* ``MULTI_WEIGHT_R0L0_PSEUDOSCALARS`` (default 1.0),

* ``MULTI_WEIGHT_R0L0_VECTORS`` (default 2.5 for the CSS shower, 2.2 for Dire),

* ``MULTI_WEIGHT_R0L0_TENSORS2`` (default 1.5),

* ``MULTI_WEIGHT_R0L1_SCALARS`` (default 0.0),

* ``MULTI_WEIGHT_R0L1_AXIALVECTORS`` (default 0.0),

* ``MULTI_WEIGHT_R0L2_VECTORS`` (default 0.5),

* ``MULTI_WEIGHT_R0L0_N_1/2`` (default 1.0),

* ``MULTI_WEIGHT_R1L0_N_1/2`` (default 0.1),

* ``MULTI_WEIGHT_R2L0_N_1/2`` (default 0.0),

* ``MULTI_WEIGHT_R1_1L0_N_1/2`` (default 0.0),

* ``MULTI_WEIGHT_R0L0_DELTA_3/2`` (default 0.15),

Note that ``MULTI_WEIGHT_R0L2_VECTORS`` only switches on the L=2 vector
(spin 1) meson multiplet; the corresponding L=2 tensor (spin 2) mesons,
e.g. :math:`K_2(1820)`, have no weight parameter of their own yet and are
currently never produced by AHADIC++.

In addition, there is a suppression factor applied to meson singlets,

* ``SINGLET_MODIFIER`` (default 2.0),

and enhancement/suppression factors applied to specific baryon and
heavy-flavour combinations,

* ``SINGLETBARYON_MODIFIER`` (default 1.80),

* ``CHARMBARYON_ENHANCEMENT`` (default 8.00),

* ``BEAUTYBARYON_ENHANCEMENT`` (default 0.80),

* ``CHARMSTRANGE_ENHANCEMENT`` (default 2.00),

* ``BEAUTYSTRANGE_ENHANCEMENT`` (default 1.40), and

* ``BEAUTYCHARM_ENHANCEMENT`` (default 1.00).

For the singlet suppression, Sherpa also allows to redefine the mixing angles
through parameters such as

* ``Mixing_0+`` (default -14.1/180*M_PI),

* ``Mixing_1-`` (default 36.4/180*M_PI),

* ``Mixing_2+`` (default 27.0/180*M_PI),

* ``Mixing_3-`` (default 0.5411), and

* ``Mixing_4+`` (default 0.6283).

The latter two, ``Mixing_3-`` and ``Mixing_4+``, currently have no effect:
AHADIC++ does not yet construct any meson multiplet with the corresponding
J=3/J=4 quantum numbers they would apply to.

And finally, some modifiers are applied to individual hadrons, with defaults
depending on whether the CSS or the Dire parton shower is used:

* ``ETA_MODIFIER`` (default 2.2 for the CSS shower, 2.82 for Dire), and

* ``ETA_PRIME_MODIFIER`` (default 4.5 for the CSS shower, 2.03 for Dire).

Cluster transition to hadrons - flavour part
--------------------------------------------

.. index:: STRANGE_FRACTION
.. index:: BARYON_FRACTION
.. index:: P_QS_by_P_QQ_norm
.. index:: P_SS_by_P_QQ_norm
.. index:: P_QQ1_by_P_QQ0
.. index:: DIRECT_TRANSITIONS

The phase space effects due to these masses govern to a large extent
the flavour content of the non-perturbative gluon splittings at the
end of the parton shower and in the decay of clusters.  They are
further modified by relative probabilities with respect to the
production of up/down flavours through the parameters

* ``STRANGE_FRACTION`` (default 0.46), the fraction of strange vs.
  up/down quarks popped from the vacuum,

* ``BARYON_FRACTION`` (default 0.17), the fraction of diquark-antidiquark
  vs. quark-antiquark pairs popped from the vacuum,

* ``P_QS_by_P_QQ_norm`` (default 0.56), entering the probability
  :math:`P_{QS}/P_{QQ}` to pop a diquark containing one strange quark
  relative to an all up/down diquark, given by this parameter multiplied
  by ``STRANGE_FRACTION``,

* ``P_SS_by_P_QQ_norm`` (default 0.056), entering the probability
  :math:`P_{SS}/P_{QQ}` to pop a doubly-strange diquark relative to an all
  up/down diquark, given by this parameter multiplied by the square of
  ``STRANGE_FRACTION``, and

* ``P_QQ1_by_P_QQ0`` (default 0.60), the relative probability of popping a
  spin-1 vs. a spin-0 diquark.

Whether a cluster may transition directly onto its lightest same-flavour
hadron at all (rather than always going through a two-cluster or
two-hadron splitting) is controlled by

* ``DIRECT_TRANSITIONS`` (default 1, i.e. enabled).

The transition of clusters to hadrons is governed by the following
considerations:

* Clusters can be interpreted as excited hadrons, with a continuous
  mass spectrum.

* When a cluster becomes sufficiently light such that its mass is
  below the largest mass of any hadron with the same flavour content,
  it must be re-interpreted as such a hadron.  In this case it will be
  shifted on the corresponding hadron mass, and the recoil will be
  distributed to the "neighbouring" clusters or by emitting a soft
  photon.  This comparison of masses clearly depends on the multiplets
  switched on in AHADIC++.

* In addition, clusters may becomes sufficiently light such that they
  should decay directly into two hadrons instead of two clusters.
  This decision compares the cluster mass to a threshold mass interpolated
  between the lightest and heaviest two-hadron final state accessible to
  its flavour content, with the interpolation fraction given by

  * ``DECAY_THRESHOLD`` (default 0.02, dimensionless: 0 picks the lightest
    accessible two-hadron mass, 1 the heaviest).

* If both options, transition and decay, are available, there is a
  competition between them, decided analogously via an interpolation
  between the lightest and heaviest single-hadron transition mass,

  * ``TRANSITION_THRESHOLD`` (default 0.51, dimensionless, same
    interpolation convention as ``DECAY_THRESHOLD``).

* A small number of further mass thresholds steer rarely-relevant edge
  cases in the radiative treatment of very light two-gluon clusters,
  deciding between a two-photon, a pion-plus-photon, or a two-pion final
  state:

  * ``PI_PHOTON_THRESHOLD`` (default 0.150 GeV),

  * ``DI_PION_THRESHOLD`` (default 0.300 GeV), and

  * ``OPEN_THRESHOLD`` (default 0.100 GeV) -- currently has no effect in
    the code.


Cluster transition and decay weights
------------------------------------

.. index:: MASS_EXPONENT
.. index:: PROMPT_DECAY_EXPONENT

The probability for a cluster C to be transformed into a hadron H is given by
a combination of weights, obtained from the overlap with the flavour part of
the hadronic wave function, the relative weight of the corresponding multiplet
and a kinematic weight taking into account the mass difference of cluster
and hadron and the width of the latter.

For the direct decay of a cluster into two hadrons the overlaps with the
wave functions of all hadrons, their respective multiplet suppression weights,
the flavour weight for the creation of the new flavour q and a kinematical
factor are relevant.  Here, yet another tuning parameter enters,

* ``MASS_EXPONENT`` (default 0.0),

which partially compensates phase space effects favouring light hadrons. A
related exponent enters the kinematic weight for a cluster's direct, prompt
transition onto a single hadron (cf. ``DIRECT_TRANSITIONS`` above),

* ``PROMPT_DECAY_EXPONENT`` (default -1.0).

Cluster decays - kinematics
---------------------------

.. index:: KT_0
.. index:: PT_MAX
.. index:: KT_ORDER
.. index:: MIN_MASS2

Cluster decays are generated by firstly emitting a non-perturbative
"gluon" from one of the quarks, using a transverse momentum
distribution as in the non-perturbative gluon decays, see below, and
by then splitting this gluon into a quark--antiquark of
anti-diquark--diquark pair, again with the same kinematics.  In the
first of these splittings, the emission of the gluon, though, the
energy distribution of the gluon is given by the quark splitting
function, if this quark has been produced in the perturbative phase of
the event.  If, in contrast, the quark stems from a cluster decay, the
energy of the gluon is selected according to a flat distribution.

In clusters decaying to hadrons, the magnitude of the transverse momentum
is drawn from a half-Gaussian (folded normal) distribution of width

* ``KT_0`` (default 1.21 GeV),

resampled until it falls below an upper limit set by the energy available
to the splitting and, by default, by

* ``PT_MAX`` (default 0.68 GeV).

Since ``PT_MAX`` is smaller than ``KT_0`` by default, it is usually this
cutoff, rather than the width of the underlying Gaussian, that determines
the bulk of the accepted transverse-momentum range. If transverse-momentum
ordering of successive splittings is switched on via

* ``KT_ORDER`` (default 0, i.e. off),

each new transverse momentum is additionally required not to exceed the
transverse momentum of the splitting that produced its parent parton.
Finally,

* ``MIN_MASS2`` (default 0.10 GeV²)

enters as a small mass-squared offset in determining the lightest cluster
that can still be processed as such, rather than transitioning directly
onto a hadron.

Splitting kinematics
--------------------

.. index:: GLUON_DECAY_MODE
.. index:: ALPHA_G
.. index:: CLUSTER_SPLITTING_MODE
.. index:: REMNANT_CLUSTER_MODE
.. index:: ALPHA_L
.. index:: BETA_L
.. index:: GAMMA_L
.. index:: ALPHA_D
.. index:: BETA_D
.. index:: GAMMA_D
.. index:: ALPHA_B
.. index:: BETA_B
.. index:: GAMMA_B
.. index:: ALPHA_H
.. index:: BETA_H
.. index:: GAMMA_H

In each splitting, the kinematics is given by the transverse momentum,
the energy splitting parameter and the azimuthal angle.  The latter,
the azimuthal angle is always selected according to a flat
distribution, while the energy splitting parameter will either be
chosen according to the quark-to-gluon splitting function (if the
quark is a leading quark, i.e. produced in the perturbative phase), to
the gluon-to-quark splitting function, or according to a flat
distribution.  The transverse momentum is given by the same
distribution as in the cluster decays to hadrons.

For the non-perturbative splitting of a gluon at the end of the parton
shower into a quark-antiquark or diquark-antidiquark pair, the energy
splitting parameter is accepted with a weight following one of two forms,

* ``GLUON_DECAY_MODE`` (default 0): :math:`z^\alpha+(1-z)^\alpha`, which
  does not vanish at the kinematic endpoints :math:`z\to0,1` unless
  :math:`\alpha` is large, or

* mode 1: :math:`[z(1-z)]^\alpha`, which always vanishes at both endpoints
  for :math:`\alpha>0`,

using a single exponent ``ALPHA_G`` (default 0.97).

For the splitting of a cluster into two new clusters, a richer,
Beta-distribution-like weight :math:`z^\alpha(1-z)^\beta` is used, with the
exponents chosen depending on the type of constituent at each end of the
splitting -- a light quark, a "leading" parton produced in the perturbative
phase, a diquark, or a beam-remnant particle:

* ``ALPHA_L``/``BETA_L``/``GAMMA_L`` (defaults 3.9/0.18/0.48) for light
  quarks (the default category),

* ``ALPHA_H``/``BETA_H``/``GAMMA_H`` (defaults -0.6/1.8/0.024) for leading
  partons,

* ``ALPHA_D``/``BETA_D``/``GAMMA_D`` (defaults 3.4/0.72/0.77) for diquarks,
  and

* ``ALPHA_B``/``BETA_B``/``GAMMA_B`` (defaults 14.2/1.6/8.1) for
  beam-remnant particles.

In the default mode, an additional exponential term proportional to
:math:`\gamma\,(k_T^2+m^2)/\texttt{KT\_0}^2` further suppresses
configurations where the daughter cluster's transverse mass is large
compared to ``KT_0``, in the same spirit as the transverse-mass suppression
term in the Lund symmetric fragmentation function. Which of several
different ways of sampling the cluster-splitting energy fractions is used
is chosen by

* ``CLUSTER_SPLITTING_MODE`` (default 2: a flat sample with no additional
  weighting beyond the standard acceptance cuts), with a separate value
  used whenever one of the two split constituents is a beam-remnant
  particle,

* ``REMNANT_CLUSTER_MODE`` (default 2).

.. _Hadron decays:

Hadron decays
=============

.. index:: PARTICLE_DATA_Mass
.. index:: PARTICLE_DATA_Width
.. index:: PARTICLE_DATA_Stable
.. index:: HADRON_DECAYS
.. index:: HADRON_DECAYS_Model
.. index:: Max_Proper_Lifetime
.. index:: Mass_Smearing
.. index:: QED_Corrections
.. index:: Spin_Correlations

The treatment of hadron and tau decays is steered by the parameters in a block
named ``HADRON_DECAYS``, e.g.

.. code-block:: yaml

   HADRON_DECAYS:
     Model: HADRONS++
     Max_Proper_Lifetime: 10.0
     QED_Corrections: 1

* Hadron properties like mass, width, and active can
  be set in full analogy to the settings for fundamental particles
  using :option:`PARTICLE_DATA`, cf. :ref:`Models`.

* ``Max_Proper_Lifetime: [mm]`` (default: 10.0) Parameter for maximum proper lifetime
  (in mm) up to which hadrons are considered unstable.
  This will make long-living particles stable, even if they are set
  unstable by default or by the user. If you do not want to set this globally,
  set this to a value of -1 and steer the stability
  through :option:`PARTICLE_DATA:<id>:Stable`, cf. :ref:`Models`.

* ``QED_Corrections: [0,1]`` (default: 1) Whether to dress hadron decays
  with QED corrections.


* ``Model: [HADRONS++, Off]`` (default: :option:`HADRONS++`)
  It defaults to :option:`Hadrons` to employ Sherpa's built-in hadron
  decay module HADRONS++ described below.
  Another option is to use the hadron decays from Pythia8 directly in the
  corresponding hadronisation interface, cf. :ref:`Fragmentation` above.
  To disable hadron decays completely, it can be disabled with the option :option:`Off`.

:option:`HADRONS++` is the built-in module within the Sherpa framework which is
responsible for treating hadron and tau decays.  It contains decay
tables with branching ratios for approximately 2500 decay channels, of
which many have their kinematics modelled according to a matrix
element with corresponding form factors.  Especially decays of the tau
lepton and heavy mesons have form factor models similar to dedicated
codes like Tauola :cite:`Jadach1993hs` and EvtGen :cite:`Lange2001uf`.

Its settings are also steered within the ``HADRON_DECAYS`` block as follows:

* ``Mass_Smearing: [0,1,2]`` (default: 1) Determines whether
  particles entering the hadron decay event phase should be put
  off-shell according to their mass distribution. It is taken care
  that no decay mode is suppressed by a potentially too low
  mass. HADRONS++ determines this dynamically from the chosen
  decay channel. Choosing option 2 instead of 1 will only set
  unstable (decayed) particles off-shell, but leave stable particles
  on-shell.

* ``Spin_Correlations: [0,1]`` (default: 0)
  A spin correlation algorithm is implemented and can be switched on with
  this setting. This might slow down event generation slightly.

* ``Channels:``
  Many aspects of the decay tables and individual decay channels can be adjusted
  within this sub-block. The default settings of the Sherpa hadron decay data
  can be found in ``<prefix>/share/SHERPA-MC/Decaydata.yaml`` and can be
  overwritten individually in the run card, e.g. as follows:

  .. code-block:: yaml

     HADRON_DECAYS:
       Channels:
         111:
           22,22:
             BR: [0.98823, 0.00034]
             Origin: PDG2023
         15:
           16,-12,11:
             BR: [0.1782, 0.0004]
             Status: [1, 2, 1]

  The levels are structured first by decaying particle and then by decay
  products. For each decay channel the following settings are available:

  * ``BR: [<br>, <deltabr>]`` branching ratio and its uncertainty. Note that the effective BR in the simulation might differ from this value if the specified BRs of a given decayer do not add up to 100\%.

  * ``Origin: <...>`` origin of BR for documentation purposes

  * ``Status:`` specifies whether the decay channel is enabled or disabled. It can take on values of ``Status: 0`` (disabled), ``Status: 1`` (enabled) and ``Status: 2`` (forced, i.e. all other decay channels of the given decayer will be disabled). If you specify multiple values, the n'th status will apply to the n'th decay of this particle flavour within the event.

  * ``ME:`` lists the matrix elements used for the decay kinematics
    and the permutation that maps the external momenta of the decay into the
    internal convention in the ME implementation.
    Additionally, parameters for the ME calculation can be specified. Example:

    .. code-block:: yaml

       HADRON_DECAYS:
         Channels:
           521:
             321,11,-11:
               BR: [5.5e-07, 7e-08]
               Origin: PDG
               ME:
                 - B_K_Semileptonic[0,1,2,3]:
                     Factor: [1.0, 0.0]
                     LD: 0
                     C1: -0.248
                     C2: 1.107
                     C3: 0.011
                     C4: -0.026
                     C5: 0.007
                     C6: -0.031
                     C7eff: -0.313
                     C9: 4.344
                     C10: -4.669

    If no ME information is specified, Sherpa will fall back to a generic
    matrix element based on the spins of the external particles.

    One special type of ME used very often is :option:`Current_ME` which
    corresponds to the contraction of two (V-A) currents that then have to
    be specified separately and can contain form factors etc. This structure
    allows to combine known currents flexibly without needing to implement
    a dedicated ME for each of these decays.
    Examples are semileptonic B/D-decays which can contain a leptonic current
    and a hadronic one or tau decays which can contain either two leptonic
    currents or also one hadronic one. Syntax example:

    .. code-block:: yaml

       HADRON_DECAYS:
         Channels:
           521:
             -423,12,-11:
               BR: [0.0558, 0.0022]
               Origin: PDG2022
               ME:
                 - Current_ME:
                     J1:
                       Type: VA_F_F
                       Indices: [2,3]
                     J2:
                       Type: VA_P_V
                       Indices: [0,1]
                       FORM_FACTOR: 3

             -411,211,12,-11:
               BR: [0.0002, 0.0002]
               Origin: FS
               ME:
                 - Current_ME:
                     J1:
                       Type: VA_F_F
                       Indices: [3,4]
                     J2:
                       Type: VA_B_DPi
                       Indices: [0,1,2]
                       Vxx: 0.04

  * ``PhaseSpace`` lists the phase-space mappings and optionally their (relative) weights.
    Example:

    .. code-block:: yaml

       PhaseSpace:
         - TwoResonances_a(1)(1260)+_2_rho(770)+_13:
             Weight: 0.5
         - TwoResonances_a(1)(1260)+_3_rho(770)+_12:
             Weight: 0.5

  * ``CPAsymmetryS:`` For CP violation in the interference between mixing and decay, cf. below.
  * ``CPAsymmetryC:`` For CP violation in the interference between mixing and decay, cf. below.
  * ``IntResults:``   This line stores the results from the phase space
    integration of the decay channel (width, MC uncertainty, maximum for unweighting).
    If they are missing, HADRONS++ integrates this channel during the initialization.

    Consequently, if some parameters are changed (also masses of
    incoming and outgoing particles) the maximum might change such that
    a new integration is needed in order to obtain correct kinematical
    distributions. In this case the ``IntResults`` line should be removed and
    replaced by the new one printed out to screen after integration.

* ``Constants`` Some globally used constants

* ``Aliases`` Create alias particles, e.g. to enforce specific decay chains. Example:

  .. code-block:: yaml

     HADRON_DECAYS:
       Aliases:
         999521: 521

       Channels:
         300553:
           999521,-999521:
             BR: 0.5
             [...]
           511,-511:
             BR: 0.5
             [...]

         999521:
           -423,12,-11:
             BR: [0.0558, 0.0022]
             Status: 2
             [...]


* ``Mixing:`` This block contains globally needed parameters for neutral meson mixing.
  Setting ``Mixing_<...> = 1`` enables explicit mixing in the event record according
  to the time evolution of the flavour states.
  The ``Interference_X = 1`` switch would enable rate asymmetries due to CP
  violation in the interference between mixing and decay (cf. ``CPAsymmetry``
  settings below). By default, the mixing parameters are set to the following values:

  .. code-block:: yaml

     HADRON_DECAYS:
       Mixing:
         Mixing_D: 1
         Interference_D: 0
         x_D: 0.0032
         y_D: 0.0069
         qoverp2_D: 1.0

         Mixing_B: 1
         Interference_B: 0
         x_B: 0.770
         y_B: 0.0
         qoverp2_B: 1.0

         Mixing_B(s): 1
         Interference_B(s): 0
         x_B(s): 26.72
         y_B(s): 0.130
         qoverp2_B(s): 1.0

  If one wants to include time dependent CP asymmetries through interference
  between mixing and decay one can set the coefficients of the cos and sin terms
  respectively for each decay channel as described above (``CPAsymmetryS/C``).
  HADRONS++ will then respect these asymmetries between particle and
  anti-particle in the choice of decay channels.

* ``Partonics`` Some partonic decay tables (for c and b) that will be used to
  complement the decay table of hadrons if they don't contain 100% BR and have
  spectators specified in their own setup like:

  .. code-block:: yaml

     521:
       Spectators: [ 2: { Weight: 1.0 } ]

* ``CreateBooklet: true`` to create a Latex booklet of all decay channels read in.

