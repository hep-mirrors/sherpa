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
.. index:: MSTJ
.. index:: MSTP
.. index:: MSTU
.. index:: PARP
.. index:: PARJ
.. index:: PARU

The ``FRAGMENTATION`` parameter sets the fragmentation module to be
employed during event generation.

* The default is :option:`Ahadic`, enabling Sherpa's native
  hadronisation model AHADIC++ :cite:`Chahal2022rid`, based on
  the cluster fragmentation model introduced in :cite:`Field1982dg`,
  :cite:`Webber1983if`, :cite:`Gottschalk1986bv`, and :cite:`Marchesini1987cf`.

* The hadronisation can be disabled with the value :option:`None`.

* To evaluate uncertainties stemming from the hadronisation, Sherpa
  also provides an interface to the Lund string fragmentation in
  Pythia 8.3 :cite:`Sjostrand2015` by using the setting
  :option:`Pythia8`.  In this case, the standard Pythia settings
  can be used to steer the behaviour of the Lund string,
  see :cite:`Sjostrand2015`. They are specified in their usual
  form in Pythia in a dedicated settings block. Additionally
  a choice can be made to let Pythia directly handle hadron
  decays via the :option:`DECAYS` setting (separate from the
  :option:`Model` switch mentioned below) and whether Pythias or
  Sherpas default masses and widths should be used through the
  :option:`SHERPA_MASSES` setting. By default the choice of generator
  for the masses and widths setting aligns with the decay setting.

.. code-block:: yaml

   SHERPA_LDADD: SherpaPythia
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
.. index:: SINGLET_SUPPRESSION
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

* ``MULTI_WEIGHT_R0L0_VECTORS`` (default 1.0),

* ``MULTI_WEIGHT_R0L0_TENSORS2`` (default 0.75),

* ``MULTI_WEIGHT_R0L1_SCALARS`` (default 0.0),

* ``MULTI_WEIGHT_R0L1_AXIALVECTORS`` (default 0.0),

* ``MULTI_WEIGHT_R0L2_VECTORS`` (default 0.0),

* ``MULTI_WEIGHT_R0L0_N_1/2`` (default 1.0),

* ``MULTI_WEIGHT_R1L0_N_1/2`` (default 0.0),

* ``MULTI_WEIGHT_R2L0_N_1/2`` (default 0.0),

* ``MULTI_WEIGHT_R1_1L0_N_1/2`` (default 0.0),

* ``MULTI_WEIGHT_R0L0_DELTA_3/2`` (default 0.25),

In addition, there is a suppression factors applied to meson singlets,

* ``SINGLET_SUPPRESSION`` (default 1.0).

For the latter, Sherpa also allows to redefine the mixing angles
through parameters such as

* ``Mixing_0+`` (default -14.1/180*M_PI),

* ``Mixing_1-`` (default 36.4/180*M_PI),

* ``Mixing_2+`` (default 27.0/180*M_PI),

* ``Mixing_3-`` (default 0.5411),

* ``Mixing_4+`` (default 0.6283),

And finally, some modifiers are applied to individual hadrons:

* ``ETA_MODIFIER`` (default 0.12),

* ``ETA_PRIME_MODIFIER`` (default 1.0),

Cluster transition to hadrons - flavour part
--------------------------------------------

.. index:: STRANGE_FRACTION
.. index:: BARYON_FRACTION
.. index:: CHARM_BARYON_MODIFIER
.. index:: BEAUTY_BARYON_MODIFIER
.. index:: P_{QS}/P_{QQ}
.. index:: P_{SS}/P_{QQ}
.. index:: P_{QQ_1}/P_{QQ_0}

The phase space effects due to these masses govern to a large extent
the flavour content of the non-perturbative gluon splittings at the
end of the parton shower and in the decay of clusters.  They are
further modified by relative probabilities with respect to the
production of up/down flavours through the parameters

* ``STRANGE_FRACTION`` (default 0.42),

* ``BARYON_FRACTION`` (default 1.0),

* ``CHARM_BARYON_MODIFIER`` (default 1.0),

* ``BEAUTY_BARYON_MODIFIER`` (default 1.0),

* ``P_{QS/P_{QQ}}`` (default 0.2),

* ``P_{SS/P_{QQ}}`` (default 0.04), and

* ``P_{QQ_1/P_{QQ_0}}`` (default 0.20).


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
  This decision is based on the heaviest hadrons accessible in a
  decay, modulated by another offset parameter,

  * ``DECAY_THRESHOLD`` (default 500 MeV).

* If both options, transition and decay, are available, there is a
  competition between


Cluster transition and decay weights
------------------------------------

.. index:: MassExponent_C->HH

The probability for a cluster C to be transformed into a hadron H is given by
a combination of weights, obtained from the overlap with the flavour part of
the hadronic wave function, the relative weight of the corresponding multiplet
and a kinematic weight taking into account the mass difference of cluster
and hadron and the width of the latter.

For the direct decay of a cluster into two hadrons the overlaps with the
wave functions of all hadrons, their respective multiplet suppression weights,
the flavour weight for the creation of the new flavour q and a kinematical
factor are relevant.  Here, yet another tuning parameter enters,

* ``MASS_EXPONENT`` (default 4.0)

which partially compensates phase space effects favouring light hadrons,

Cluster decays - kinematics
---------------------------

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

In clusters decaying to hadrons, the transverse momentum is chosen according
to a distribution given by an infrared-continued strong coupling and a
term inversely proportional to the infrared-modified transverse momentum,

constrained to be below a maximal transverse momentum.

Splitting kinematics
--------------------

In each splitting, the kinematics is given by the transverse momentum,
the energy splitting parameter and the azimuthal angle.  The latter,
the azimuthal angle is always selected according to a flat
distribution, while the energy splitting parameter will either be
chosen according to the quark-to-gluon splitting function (if the
quark is a leading quark, i.e. produced in the perturbative phase), to
the gluon-to-quark splitting function, or according to a flat
distribution.  The transverse momentum is given by the same
distribution as in the cluster decays to hadrons.

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

  * ``BR: [<br>, <deltabr>]`` branching ratio and its uncertainty

  * ``Origin: <...>`` origin of BR for documentation purposes

  * ``Status:`` TODO

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

