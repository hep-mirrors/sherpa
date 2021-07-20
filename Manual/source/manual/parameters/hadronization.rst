.. _Hadronization:

*************
Hadronization
*************

The hadronization setup covers the fragmentation of partons into
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
  hadronization model AHADIC++, based on the cluster fragmentation
  model introduced in :cite:`Field1982dg`, :cite:`Webber1983if`,
  :cite:`Gottschalk1986bv`, and :cite:`Marchesini1987cf` and
  implementing some modifications discussed in :cite:`Winter2003tt`.

* The hadronization can be disabled with the value :option:`None`.

* To evaluate uncertainties stemming from the hadronization, Sherpa
  also provides an interface to the Lund string fragmentation in
  Pythia 8.3 :cite:`Sjostrand2015` by using the setting
  :option:`Pythia8`.  In this case, the standard Pythia settings
  can be used to steer the behaviour  of the Lund string,
  see :cite:`Sjostrand2015`. They are specified in their usual
  form in Pythia in a dedicated settings block. Additionally
  a choice can be made to let Pythia directly handle hadron
  decays via the :option:`DECAYS` setting (separate from the
  DECAYMODEL switch mentioned below) and whether Pythias or
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
     DECAYS: false
     SHERPA_MASSES: true

* Alternatively, Sherpa  also provides an interface to
  Pythia 6.4 :cite:`Sjostrand2006za` by using the setting
  :option:`Lund`.  In this case, the standard Pythia switches
  :option:`MSTJ`, :option:`MSTP`, :option:`MSTU`, :option:`PARP`,
  :option:`PARJ` and :option:`PARU` can be used to steer the behaviour
  of the Lund string, see :cite:`Sjostrand2006za`. They can be
  specified as a 2xN matrix:

.. code-block:: yaml

   FRAGMENTATION: Lund
   MSTJ:
   - [<number1>, <value1>]
   - [<number2>, <value2>]
     ...
   MSTP:
   - [<number1>, <value1>]
     ...

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

* Clusters can be interpreted as excited hadrons, with a continous
  mass spectrum.

* When a cluster becomes sufficiently light such that its mass is
  below the largest mass of any hadron with the same flavour content,
  it must be re-iterpreted as such a hadron.  In this case it will be
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
factor are relevant.  Here, yet another tuning paramter enters,

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
term inversemly proportional to the infrared-modified transverse momentum,

constrained to be below a maximal transverse momentum.

Splitting kinematics
--------------------

In each splitting, the kinematics is given by the transverse momentum,
the energy splitting parameter and the azimuthal angle.  The latter,
the azimuthal angle is always seleectred according to a flat
distribution, while the energy splitting parameter will either be
chosen according to the quark-to-gluon splitting function (if the
quark is a leading quark, i.e. produced in the pertrubative phase), to
the gluon-to-quark splitting function, or according to a flat
distribution.  The transverse momentum is given by the same
distribution as in the cluster decays to hadrons.

.. _Hadron decays:

Hadron decays
=============

.. index:: DECAYMODEL
.. index:: WIDTH[<id>]
.. index:: MASS[<id>]
.. index:: STABLE[<id>]
.. index:: DECAYPATH
.. index:: SOFT_MASS_SMEARING
.. index:: MAX_PROPER_LIFETIME

The treatment of hadron and tau decays is specified by
``DECAYMODEL``. Its allowed values are either the default choice
:option:`Hadrons` (the default if :option:`FRAGMENTATION: Ahadic`),
:option:`Lund` (the default if :option:`FRAGMENTATION: Lund`), or it
can be disabled with the option :option:`Off`.

HADRONS++ is the module within the Sherpa framework which is
responsible for treating hadron and tau decays.  It contains decay
tables with branching ratios for approximately 2500 decay channels, of
which many have their kinematics modelled according to a matrix
element with corresponding form factors.  Especially decays of the tau
lepton and heavy mesons have form factor models similar to dedicated
codes like Tauola :cite:`Jadach1993hs` and EvtGen :cite:`Lange2001uf`.

Some general switches which relate to hadron decays are


.. _DECAYPATH:

* ``DECAYPATH`` The path to the parameter files for the hadron and tau
  decays (default: ``Decaydata/``). It is important to note that the
  path has to be given relative to the current working directory.  If
  it doesn't exist, the default Decaydata directory
  (``<prefix>/share/SHERPA-MC/Decaydata``) will be used.

* Hadron properties like mass, width, stable/unstable and active can
  be set in full analogy to the settings for fundamental particles
  using :option:`PARTICLE_DATA`, cf. :ref:`Models`.

.. _SOFT_MASS_SMEARING:

* ``SOFT_MASS_SMEARING = [0,1,2]`` (default: 1) Determines whether
  particles entering the hadron decay event phase should be put
  off-shell according to their mass distribution. It is taken care
  that no decay mode is suppressed by a potentially too low
  mass. While HADRONS++ determines this dynamically from the chosen
  decay channel, for ``Pythia`` as hadron decay handler its ``w-cut``
  parameter is employed. Choosing option 2 instead of 1 will only set
  unstable (decayed) particles off-shell, but leave stable particles
  on-shell.

.. _MAX_PROPER_LIFETIME:

* ``MAX_PROPER_LIFETIME = [mm]`` Parameter for maximum proper lifetime
  (in mm) up to which particles are considered unstable. If specified,
  this will make long-living particles stable, even if they are set
  unstable by default or by the user.

Many aspects of the above mentioned "Decaydata" can be adjusted.
There exist three levels of data files, which are explained in the following
sections.
As with all other setup files, the user can either employ the default
"Decaydata" in ``<prefix>/share/SHERPA-MC/Decaydata``, or
overwrite it (also selectively) by creating the appropriate files in the
directory specified by ``DECAYPATH``.

HadronDecays.dat
----------------

:file:`HadronDecays.dat` consists of a table of particles that are to
be decayed by HADRONS++. Note: Even if decay tables exist for the
other particles, only those particles decay that are set unstable,
either by default, or in the model/fragmentation settings. It has the
following structure, where each line adds one decaying particle:

+-------------------+---------------------+------------------+
| <kf-code>         |   <subdirectory>    |   <filename>.dat |
+-------------------+---------------------+------------------+
| decaying particle | path to decay table | decay table file |
+-------------------+---------------------+------------------+
| default names:    |     <particle>/     |       Decays.dat |
+-------------------+---------------------+------------------+

It is possible to specify different decay tables for the particle
(positive kf-code) and anti-particle (negative kf-code). If only one
is specified, it will be used for both particle and anti-particle.

If more than one decay table is specified for the same kf-code, these
tables will be used in the specified sequence during one event. The
first matching particle appearing in the event is decayed according to
the first table, and so on until the last table is reached, which will
be used for the remaining particles of this kf-code.

Additionally, this file may contain the keyword ``CREATE_BOOKLET`` on
a separate line, which will cause HADRONS++ to write a LaTeX document
containing all decay tables.

Decay table files
-------------------

The decay table contains information about outgoing particles for each channel,
its branching ratio and eventually the name of the file that stores parameters
for a specific channel. If the latter is not specified HADRONS++ will produce it and
modify the decay table file accordingly.

Additionally to the branching ratio, one may specify the error associated with
it, and its source. Every hadron is supposed to have its own decay table in
its own subdirectory. The structure of a decay table is

+--------------------+----------------------+--------------------+
| {kf1,kf2,kf3,...}  | BR(delta BR)[Origin] | <filename>.dat     |
+--------------------+----------------------+--------------------+
| outgoing particles | branching ratio      | decay channel file |
+--------------------+----------------------+--------------------+

It should be stressed here that the branching ratio which is
explicitly given for any individual channel in this file is **always
used** regardless of any matrix-element value.

.. _Decay channel files:

Decay channel files
-------------------

A decay channel file contains various information about that specific decay
channel. There are different sections, some of which are optional:

*
  .. code-block:: xml

     <Options>
         AlwaysIntegrate = 0
         CPAsymmetryC = 0.0
         CPAsymmetryS = 0.0
     </Options>

  * ``AlwaysIntegrate = [0,1]`` For each decay channel, one needs an
    integration result for unweighting the kinematics (see
    below). This result is stored in the decay channel file, such that
    the integration is not needed for each run. The AlwaysIntegrate
    option allows to bypass the stored integration result, and do the
    integration nonetheless (same effect as deleting the integration
    result).

  * ``CPAsymmetryC/CPAsymmetryS`` If one wants to include time
    dependent CP asymmetries through interference between mixing and
    decay one can set the coefficients of the cos and sin terms
    respectively.  HADRONS++ will then respect these asymmetries
    between particle and anti-particle in the choice of decay
    channels.

*
  .. code-block:: xml

     <Phasespace>
       1.0 MyIntegrator1
       0.5 MyIntegrator2
     </Phasespace>

  Specifies the phase-space mappings and their weight.

*
  .. code-block:: xml

     <ME>
       1.0 0.0 my_matrix_element[X,X,X,X,X,...]
       1.0 0.0 my_current1[X,X,...] my_current2[X,X,X,...]
     </ME>

  Specifies the matrix elements or currents used for the kinematics,
  their respective weights, and the order in which the particles
  (momenta) enter them. For more details, the reader is referred to
  :cite:`Krauss2010xx`.

*
  .. code-block:: text

     <my_matrix_element[X,X,X,X,X,...]>
       parameter1 = value1
       parameter2 = value2
       ...
     </my_matrix_element[X,X,X,X,X,...]>

  Each matrix element or current may have an additional section where
  one can specify needed parameters, e.g. which form factor model to
  choose.  Each parameter has to be specified on a new line as shown
  above. Available parameters are listed in :cite:`Krauss2010xx`.
  Parameters not specified get a default value, which might not make
  sense in specific decay channels. One may also specify often needed
  parameters in ``HadronConstants.dat``, but they will get overwritten
  by channel specific parameters, should these exist.

*
  .. code-block:: xml

     <Result>
       3.554e-11 6.956e-14 1.388e-09;
     </Result>

  These last three lines have quite an important meaning. If they are
  missing, HADRONS++ integrates this channel during the initialization
  and adds the result lines.  If this section exists though, and
  ``AlwaysIntegrate`` is off (the default value, see above) then
  HADRONS++ reads in the maximum for the kinematics unweighting.

  Consequently, if some parameters are changed (also masses of
  incoming and outgoing particles) the maximum might change such that
  a new integration is needed in order to obtain correct kinematical
  distributions. There are two ways to enforce the integration: either
  by deleting the last three lines or by setting ``AlwaysIntegrate``
  to 1. When a channel is re-integrated, HADRONS++ copies the old
  decay channel file into ``.<filename>.dat.old``.


HadronConstants.dat
--------------------


``HadronConstants.dat`` may contain some globally needed parameters
(e.g.  for neutral meson mixing, see :cite:`Krauss2010xx`) and also
fall-back values for all matrix-element parameters which one specifies
in decay channel files. Here, the ``Interference_X = 1`` switch would
enable rate asymmetries due to CP violation in the interference
between mixing and decay (cf. :ref:`Decay channel files`), and setting
``Mixing_X = 1`` enables explicit mixing in the event record according
to the time evolution of the flavour states. By default, all mixing
effects are turned off.

Mixing parameters with some example values

.. code-block:: python

   x_K = 0.946
   y_K = -0.9965
   qoverp2_K = 1.0
   Interference_K = 0
   Mixing_K = 0

   x_D = 0.0
   y_D = 0.0
   qoverp2_D = 1.0
   Interference_D = 0
   Mixing_D = 0

   x_B = 0.776
   y_B = 0.0
   qoverp2_B = 1.0
   Interference_B = 1
   Mixing_B = 0

   x_B(s) = 30.0
   y_B(s) = 0.155
   qoverp2_B(s) = 1.0
   Interference_B(s) = 0
   Mixing_B(s) = 0

Further remarks
---------------

.. index:: SOFT_SPIN_CORRELATIONS
.. index:: HARD_SPIN_CORRELATIONS

**Spin correlations:** a spin correlation algorithm is implemented. It
can be switched on through the setting
:option:`SOFT_SPIN_CORRELATIONS: 1`.

If spin correlations for tau leptons produced in the hard scattering
process are supposed to be taken into account, one needs to specify
:option:`HARD_SPIN_CORRELATIONS: 1` as well. If using AMEGIC++ as ME
generator, note that the Process libraries have to be re-created if
this is changed.

**Adding new channels:** if new channels are added to HADRONS++
(choosing isotropic decay kinematics) a new decay table must be
defined and the corresponding hadron must be added to
``HadronDecays.dat``.  The decay table merely needs to consist of the
outgoing particles and branching ratios, i.e. the last column (the one
with the decay channel file name) can safely be dropped. By running
Sherpa it will automatically produce the decay channel files and write
their names in the decay table.

**Some details on tau decays:** :math:`\tau` decays are treated within
the HADRONS++ framework, even though the :math:`\tau` is not a
hadron. As for many hadron decays, the hadronic tau decays have form
factor models implemented, for details the reader is referred to
:cite:`Krauss2010xx`.
