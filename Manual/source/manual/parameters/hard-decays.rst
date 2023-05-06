.. _Hard decays:

***********
Hard decays
***********

.. index:: HARD_DECAYS
.. index:: PARTICLE_DATA_Stable

The handler for decays of particles produced in the hard scattering
process (e.g. W, Z, top, Higgs) can be enabled and configured using
the :option:`HARD_DECAYS` collection of settings (and a small number
of other other top-level settings).  Which (anti)particles IDs should
be treated as unstable is determined by the
:option:`PARTICLE_DATA:<id>:Stable` switch described in :ref:`Models`.

The syntax to configure :option:`HARD_DECAYS` sub-settings is:

.. code-block:: yaml

   HARD_DECAYS:
     <sub-setting>: <value>
     # more sub-settings ...
     Channels:
       <channel id>:
         <channel sub-setting>: <value>
         # more sub-settings for <channel>
       # more channels ...

The central setting to enable the hard decays is

.. code-block:: yaml

   HARD_DECAYS:
     Enabled: true

The channel ID codes are of the form ``a,b,c,...``, where ``a`` is the
PDG ID of the decaying particle and ``b,c,...`` are the decay products.
The IDs for the decay channels can also be found in the decay table printed to
screen during the run.

This decay module can also be used on top of NLO matrix elements, but
it does not include any NLO corrections in the decay matrix elements
themselves.

Note that the decay handler is an afterburner at the event generation
level.  It does not affect the calculation and integration of the hard
scattering matrix elements. The cross section is thus unaffected
during integration, and the branching ratios (if any decay channels
have been disabled) are only taken into account for the event weights
and cross section output at the end of event generation (if not
disabled with the :option:`HARD_DECAYS:Apply_Branching_Ratios` option,
cf. below).  Furthermore any cuts or scale definitions are not
affected by decays and operate only on the inclusively produced
particles before decays.

.. contents::
   :local:

.. _Status:

Status
======

.. index:: Status

This sub-setting to each channel defined in :option:`HARD_DECAYS:Channels`
allows to explicitly force or disable a decay channel. The status can take the
following values:

:option:`Status: -1`
  Decay channel is disabled and does not contribute to total width.

:option:`Status: 0`
  Decay channel is disabled but contributes to total width.

:option:`Status: 1 (default)`
  Decay channel is enabled.

:option:`Status: 2`
  Decay channel is forced.

For example, to disable the hadronic decay channels of the W boson one would use:

.. code-block:: yaml

   HARD_DECAYS:
     Channels:
       24,2,-1:  { Status: 0 }
       24,4,-3:  { Status: 0 }
       -24,-2,1: { Status: 0 }
       -24,-4,3: { Status: 0 }

In the same way, the bottom decay mode of the Higgs could be forced using:

.. code-block:: yaml

   25,5,-5:  { Status: 2 }

Note that the ordering of the decay products in :option:`<channel id>` is
important and has to be identical to the ordering in the decay table
printed to screen.  It is also possible to request multiple forced
decay channels (:option:`Status: 2`) for the same particle, all other
channels will then automatically be disabled.

.. _Width:

Width
=====

.. index:: Width

This option allows to overwrite the calculated partial width (in GeV)
of a given decay channel, and even to add new inactive channels which
contribute to the total width. This is useful to adjust the branching
ratios, which are used for the relative contributions of different
channels and also influence the cross section during event generation,
as well as the total width which is used for the lineshape of the
resonance.

An example to set (/add) the partial widths of the ``H->ff``,
``H->gg`` and ``H->yy`` channels can be seen in the following. The
values have been taken from `LHC Higgs WG
<https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWG/Higgs_XSBR_YR4_update.xlsx>`_):

.. code-block:: yaml

   PARTICLE_DATA:
     25:
       Mass: 125.09
       Width: 0.0041

   HARD_DECAYS:
     Enabled: true
     Channels:
       25,5,-5:    { Width: 2.382E-03 }
       25,15,-15:  { Width: 2.565E-04 }
       25,13,-13:  { Width: 8.901E-07 }
       25,4,-4:    { Width: 1.182E-04 }
       25,3,-3:    { Width: 1E-06 }
       25,21,21:   { Width: 3.354E-04 }
       25,22,22:   { Width: 9.307E-06 }
       25,23,22:   { Width: 6.318E-06 }

Another example, setting the leptonic and hadronic decay channels of W
and Z bosons to the PDG values, would be specified as follows:

.. code-block:: yaml

   HARD_DECAYS:
     Enabled: true
     Channels:
       24,2,-1:    { Width: 0.7041 }
       24,4,-3:    { Width: 0.7041 }
       24,12,-11:  { Width: 0.2256 }
       24,14,-13:  { Width: 0.2256 }
       24,16,-15:  { Width: 0.2256 }
       -24,-2,1:   { Width: 0.7041 }
       -24,-4,3:   { Width: 0.7041 }
       -24,-12,11: { Width: 0.2256 }
       -24,-14,13: { Width: 0.2256 }
       -24,-16,15: { Width: 0.2256 }
       23,1,-1:    { Width: 0.3828 }
       23,2,-2:    { Width: 0.2980 }
       23,3,-3:    { Width: 0.3828 }
       23,4,-4:    { Width: 0.2980 }
       23,5,-5:    { Width: 0.3828 }
       23,11,-11:  { Width: 0.0840 }
       23,12,-12:  { Width: 0.1663 }
       23,13,-13:  { Width: 0.0840 }
       23,14,-14:  { Width: 0.1663 }
       23,15,-15:  { Width: 0.0840 }
       23,16,-16:  { Width: 0.1663 }
       6,24,5:     { Width: 1.32 }
       -6,-24,-5:  { Width: 1.32 }

See also :option:`Use_HO_SM_Widths` below for a global automatic switch to set these values.

.. _Use_HO_SM_Widths:

Use_HO_SM_Widths
================

.. index:: Use_HO_SM_Widths

The partial decay widths (and thus BRs) calculated and used by the decay
handler are only LO accurate. For SM setups, we provide pre-defined decay
widths taking higher-order corrections into account. By default
(:option:`HARD_DECAYS: { Use_HO_SM_Widths: true }`) these will overwrite
the LO widths with the values given in the :option:`Width` example above.


.. _HARD_DECAYS_Spin_Correlations:

Spin_Correlations
=================

.. index:: Spin_Correlations

Spin correlations between the hard scattering process and the
following decay processes are enabled by default. If you want to
disable them, e.g. for spin correlation studies, you can specify the
option :option:`Spin_Correlations: 0`.

.. _Store_Results:

Store_Results
=============

.. index:: Store_Results

The decay table and partial widths are calculated on-the-fly during
the initialization phase of Sherpa from the given model and its
particles and interaction vertices. To store these results in the
``Results/Decays`` directory, one has to specify :option:`HARD_DECAYS:
{ Store_Results: 1 }`.  In case existing decay tables are to be read
in the same configuration should be done. Please note, that Sherpa
will delete decay channels present in the read in results but not in
the present model with present parameters by default. To prevent
Sherpa from updating the decay table files accordingly specify
:option:`HARD_DECAYS: { Store_Results: 2 }`.

.. _hard_Result_Directory:

Result_Directory
================

.. index:: Result_Directory

Specifies the name of the directory where the decay results are to be
stored. Defaults to the value of the top-level setting
:ref:`RESULT_DIRECTORY`.

.. _Set_Widths:

Set_Widths
==========

.. index:: Set_Widths
.. index:: PARTICLE_DATA_Width

The decay handler computes LO partial and total decay widths and
generates decays with corresponding branching fractions, independently
from the particle widths specified by
:option:`PARTICLE_DATA:<id>:Width`. The latter are relevant only for
the core process and should be set to zero for all unstable particles
appearing in the core-process final state. This guarantees
on-shellness and gauge invariance of the core process, and subsequent
decays can be handled by the afterburner.  In constrast,
:option:`PARTICLE_DATA:<id>:Width` should be set to the physical width
when unstable particles appear (only) as intermediate states in the
core process, i.e. when production and decay are handled as a full
process or using ``Decay``/``DecayOS``.  In this case, the option
:option:`HARD_DECAYS: { Set_Widths: true }` permits to overwrite the
:option:`PARTICLE_DATA:<id>:Width` values of unstable particles by the
LO widths computed by the decay handler.

.. _Apply_Branching_Ratios:

Apply_Branching_Ratios
======================

.. index:: Apply_Branching_Ratios

By default (:option:`HARD_DECAYS: { Apply_Branching_Ratios: true }`),
weights for events which involve a hard decay are multiplied with the
corresponding branching ratios (if decay channels have been
disabled). This also means that the total cross section at the end of
the event generation run already includes the appropriate BR
factors. If you want to disable that, e.g. because you want to
multiply with your own modified BR, you can set the option
:option:`{HARD_DECAYS: { Apply_Branching_Ratios: false }`.

.. _Mass_Smearing:

Mass_Smearing
=============

.. index:: Mass_Smearing

With the default of :option:`HARD_DECAYS: { Mass_Smearing: 1 }` the
kinematic mass of the unstable propagator is distributed according to
a Breit-Wigner shape a posteriori. All matrix elements are still
calculated in the narrow-width approximation with onshell
particles. Only the kinematics are affected.  To keep all intermediate
particles onshell :option:`{HARD_DECAYS: { Mass_Smearing: 0 }`.

.. _Resolve_Decays:

Resolve_Decays
==============

.. index:: Resolve_Decays
.. index:: Min_Prop_Width

There are different options how to decide when a 1->2 process should
be replaced by the respective 1->3 processes built from its decaying
daughter particles.

:option:`Resolve_Decays: Threshold`
  (default)
  Only when the sum of decay product masses exceeds the decayer mass.

:option:`Resolve_Decays: ByWidth`
  As soon as the sum of 1->3 partial widths exceeds the 1->2 partial width.

:option:`Resolve_Decays: None`
  No 1->3 decays are taken into account.

In all cases, one can exclude the replacement of a particle below a
given width threshold using :option:`Min_Prop_Width: <threshold>`
(default 0.0).  Both settings are sub-settings of
:option:`HARD_DECAYS`:

.. code-block:: yaml

   HARD_DECAYS:
     Resolve_Decays: <mode>
     Min_Prop_Width: <threshold>

.. _Decay_Tau:

Decay_Tau
=========

.. index:: Decay_Tau

By default, the tau lepton is decayed by the hadron decay module,
:ref:`Hadron decays`, which includes not only the leptonic decay
channels but also the hadronic modes. If :option:`Decay_Tau: true` is
specified, the tau lepton will be decayed in the hard decay handler,
which only takes leptonic and partonic decay modes into account. Note,
that in this case the tau needs to also be set massive:

.. code-block:: yaml

   PARTICLE_DATA:
     15:
       Massive: true
   HARD_DECAYS:
     Decay_Tau: true

.. _Decay table integration settings:

Decay table integration settings
================================

.. index:: Int_Accuracy
.. index:: Int_Target_Mode
.. index:: Int_NIter

Three parameters can be used to steer the accuracy and time
consumption of the calculation of the partial widths in the decay
table: :option:`Int_Accuracy: 0.01` specifies a relative accuracy for
the integration. The corresponding target reference is either the
given total width of the decaying particle (:option:`Int_Target_Mode:
0`, default) or the calculated partial decay width
(:option:`Int_Target_Mode: 1`). The option :option:`Int_NIter: 2500`
can be used to change the number of points per integration iteration,
and thus also the minimal number of points to be used in an
integration.  All decay table integration settings are sub-settings of
:option:`HARD_DECAYS`.

.. _Simulation of polarized cross sections for intermediate particles:


Simulation of polarized cross sections for intermediate particles
=================================================================

.. index:: Simulation of polarized cross sections for intermediate particles

In this chapter it is described how Sherpa can be used to simulate
polarized cross sections for unstable intermediate state particles.
At the moment, only the simulation of polarized
cross sections for massive vector bosons is supported.
Sherpa can simulate all polarized cross sections in one simulation
run. The polarized cross sections are handled during event 
generation and printed out as additional event weights similar to 
variation weights.   
By default, the cross sections for all polarization combinations of 
the intermediate particles are printed out. For massive vector 
bosons also all transverse weights are calculated automatically. 
Beside this, user-specified weights can be determined which is described
in section :ref:`Custom polarized cross sections`.
Weight names for automatically provided weights have the form
``PolWeight_ReferenceSystem.particle1.helicity1_particle2.helicity2...``
e.g. for two W+-bosons PolWeight.W+.+_W+.-. Ordering of the particles
in the weight name corresponds to Sherpa's internal particle ordering
which can be read off from the ordering in the process printed out
when Sherpa starts running. 
More details about the definition of polarization for intermediate 
vector bosons and the implementation can be found in     .

.. contents::
   :local:

.. _General procedure:

General procedure
------------------

.. index:: General procedure

The definition of polarization for particles in intermediate states is only 
possible for processes which can be factorized into a production matrix
element and decay matrix elements of them. 
To neglect possible non-resonant processes, for which this factorization
is not possible, Sherpa applies an extended narrow-width approximation. All intermediate
particles are considered as on shell but all spin correlations are preserved. 
The production process is specified in the :ref:`Processes` part 
whereas the possible decays are characterized in the :ref:`Hard decays` section.
Details about :option:`PROCESSES` and :option:`HARD_DECAYS` definition are described
in the corresponding chapters of this manual. The following example shows 
:option:`PROCESSES` and :option:`HARD_DECAYS` definition of the same-sign :math:`W^+ W^+`-scattering
with the W decaying to electrons or muons.

.. code-block:: yaml

   PARTICLE_DATA:
     24: 
       Width: 0

   HARD_DECAYS:
     Enabled: true
     Mass_Smearing: 1
     Channels:
      24,12,-11: {Status: 2}
      24,14,-13: {Status: 2}

   PROCESSES:
   - 93 93 -> 24 24 93 93:
      Order: {QCD: 0, EW: 4} 

Things to notice:

* In :option:`PARTICLE_DATA` the :option:`Width` of the intermediate particles must
  be set to zero since they are handled as stable for the hard process
  matrix element calculation. The particles are then decayed by the internal
  (hard) decay module.

* Spin correlations must be enabled during the hard decays (which is the default).

The central setting to enable the calculation of polarized cross sections is:

.. code-block:: yaml

   HARD_DECAYS:
     Pol_Cross_Section: 
       Enabled: true

The polarization vectors of massive vector bosons are implemented according to
:cite:`Dittmaier1998nn`, equation (3.19). Specifically, the polarization vectors 
are expressed in terms of Weyl spinors. For that, an arbitrary light-like vector
needs to be chosen. The definition of vector boson polarization
is not unambiguous. It can be specified by the options described in the following 
sections: 
:option:`Pol_Cross_Section:Spin_Basis` and 
:option:`Pol_Cross_Section:Reference_System`.  

.. _Spin basis:

Spin basis
----------

.. index:: Spin basis

For massive vector bosons the choice of a light-like vector for their description 
in the Weyl spinor formalism is not really arbitrary since it characterizes the 
spin axis chosen to define the polarization.
By default, the reference vector is selected such that polarization vectors are expressed 
in the helicity basis since this is the common choice for vector boson polarization.
The polarization vectors are then eigenvectors of the helicity operator and have the same
form as in (3.15) in :cite:`Dittmaier1998nn` after transformation from spinor to vector
representation. Sherpa provides several gauge choices for the Weyl spinors. To really get
the polarization vectors in this form, the following spinor gauge choice must be chosen:

.. code-block:: yaml

   COMIX_DEFAULT_GAUGE: 0

.. code-block:: yaml

   HARD_DECAYS:
     Pol_Cross_Section: 
       Enabled: true
       Spin_Basis: Helicity

If :option:`Spin_Basis: ComixDefault` is selected the COMIX default reference vector 
specified by :option:`COMIX_DEFAULT_GAUGE` (default 1 which corresponds to
(1.0, 0.0, 1/:math:`\sqrt{2}`, 1/:math:`\sqrt{2}`)) is used. 
Furthermore, it is possible to hand over any constant reference vector:

.. code-block:: yaml

   HARD_DECAYS:
     Pol_Cross_Section:
       Enabled: true
       Spin_Basis: 1.0, 0.0, 1.0, 0.0

.. _Reference system:

Reference system
-----------------

.. index:: Reference system

The helicity of a massive particle is not Lorentz invariant. Therefore a reference system
needs to be chosen to define vector boson polarization unambiguous. Sherpa supports the 
following options: 

:option:`Reference_System: Lab` (default)
  Vector boson polarization is defined in the laboratory frame. 

:option:`Reference_System: COM`
  Vector boson polarization is defined in the center of mass system of all 
  hard-decaying particles.

:option:`Reference_System: PPFr`
  Vector boson polarization is defined in the center of mass system of the two interacting
  partons.

:option:`Reference_System: RestFrames`
  Vector boson polarization for each hard decaying particle is defined in its own
  center of mass system. If the helicity basis is selected as spin basis the spin axis lies along the flight
  direction of the hard decaying particle in the laboratory frame.

Sherpa allows the calculation of polarized cross sections for different polarization 
definitions specified by different reference systems in one simulation run:

.. code-block:: yaml

   HARD_DECAYS:
     Pol_Cross_Section: 
       Enabled: true
       Reference_System: [Lab, COM]

Additionally to the options explained above, each reference system defined by one or several
hard process initial or final state particles can be used. This can be specified by the 
particle numbers of the desired particles according to the Sherpa numbering scheme.
Distinct particle numbers should only be separated by a single white space, at least if more 
than one reference system is specified. The second reference frame in the following example
is the parton-parton rest frame.

.. code-block:: yaml

   HARD_DECAYS:
     Pol_Cross_Section: 
       Enabled: true
       Reference_System: [Lab, 0 1]

In the Sherpa event output, polarized cross sections of vector bosons defined
in different frames are distinguished by adding the reference frame to the weight names, 
e.g. PolWeight_Lab.W+.+_W+.-. For reference systems defined by particle numbers,  
``refsystemn`` is added to avoid commas in weight names. n is the place in the reference 
system list specified in the YAML-File starting at 0. For the example above this means
e.g. PolWeight_refsystem1.W+.+_W+.-.

.. _Custom polarized cross sections:

Custom polarized cross sections
--------------------------------

.. index:: Custom polarized cross sections

Sherpa provides the calculation of two different types of custom polarized cross sections.
On the one hand, it is possible to hand over a comma separated list of weight names from the 
automatically calculated cross sections.
These cross sections are then added by Sherpa and printed out as additional event weight.
On the other hand, partially unpolarized cross sections can be calculated. 
These can be specified by the numbers of the particles which should be considered 
as unpolarized. Also here the numbering of the particles is according to the Sherpa numbering
scheme.  
Custom weights are generally specified by :option:`Weight`. By adding numbers to 
:option:`Weight` e.g. :option:`Weight1`, :option:`Weight2` ... more than one custom cross 
section can be calculated. The number is limited to :option:`Weight10` by default but can
be increased by using :option:`Number_Of_Custom_Weights`.
In the following example the W- is considered as unpolarized: 

.. code-block:: yaml

   HARD_DECAYS:
    Enabled: true
    Mass_Smearing: 1
    Channels:
     24,12,-11: {Status: 2}
     -24,-14,13: {Status: 2}
    Pol_Cross_Section: 
      Enabled: true
      Weight1: W+.+_W-.+, W+.-_W-.+
      Weight2: 3

   PROCESSES:
    - 93 93 -> 24 -24 93 93:
      Order: {QCD: 0, EW: 4} 

In the case of partially unpolarized cross sections the helicity of the unpolarized particle
is set to ``u`` in the weight names, e.g. PolWeight_Lab.W+.+_W-.u. For custom cross
sections specified by weight names PolWeight_refsystem. ``Weightn`` is used instead to avoid
long weight names. Hereby, :option:`Weightn` corresponds to the corresponding setting in the 
YAML-File.  
