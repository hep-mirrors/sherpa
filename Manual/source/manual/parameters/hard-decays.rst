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
       "<channel id>":
         <channel sub-setting>: <value>
         # more sub-settings for <channel>
       # more channels ...

The channel ID codes are of the form ``a -> b c ...``.  The particle
IDs for the decay channels can be found in the decay table printed to
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
       "24 -> 2 -1":  { Status: 0 }
       "24 -> 4 -3":  { Status: 0 }
       "-24 -> -2 1": { Status: 0 }
       "-24 -> -4 3": { Status: 0 }

In the same way, the bottom decay mode of the Higgs could be forced using:

.. code-block:: yaml

   "25 -> 5 -5":  { Status: 2 }

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
values have been taken from `LHCHXSWG
<https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3>`_):

.. code-block:: yaml

   PARTICLE_DATA:
     25:
       Mass: 125
       Width: 0.00407

   HARD_DECAYS:
     Enabled: true
     Channels:
       "25 -> 5 -5":   { Width: 2.35e-3 }
       "25 -> 15 -15": { Width: 2.57e-4 }
       "25 -> 13 -13": { Width: 8.91e-7 }
       "25 -> 4 -4":   { Width: 1.18e-4 }
       "25 -> 3 -3":   { Width: 1.00e-6 }
       "25 -> 21 21":  { Width: 3.49e-4 }
       "25 -> 22 22":  { Width: 9.28e-6 }

Another example, setting the leptonic and hadronic decay channels of W
and Z bosons to the PDG values, would be specified as follows:

.. code-block:: yaml

   HARD_DECAYS:
     Enabled: true
     Channels:
       "24 -> 2 -1":    { Width: 0.7041 }
       "24 -> 4 -3":    { Width: 0.7041 }
       "24 -> 12 -11":  { Width: 0.2256 }
       "24 -> 14 -13":  { Width: 0.2256 }
       "24 -> 16 -15":  { Width: 0.2256 }
       "-24 -> -2 1":   { Width: 0.7041 }
       "-24 -> -4 3":   { Width: 0.7041 }
       "-24 -> -12 11": { Width: 0.2256 }
       "-24 -> -14 13": { Width: 0.2256 }
       "-24 -> -16 15": { Width: 0.2256 }
       "23 -> 1 -1":    { Width: 0.3828 }
       "23 -> 2 -2":    { Width: 0.2980 }
       "23 -> 3 -3":    { Width: 0.3828 }
       "23 -> 4 -4":    { Width: 0.2980 }
       "23 -> 5 -5":    { Width: 0.3828 }
       "23 -> 11 -11":  { Width: 0.0840 }
       "23 -> 12 -12":  { Width: 0.1663 }
       "23 -> 13 -13":  { Width: 0.0840 }
       "23 -> 14 -14":  { Width: 0.1663 }
       "23 -> 15 -15":  { Width: 0.0840 }
       "23 -> 16 -16":  { Width: 0.1663 }

.. _HARD_SPIN_CORRELATIONS:

HARD_SPIN_CORRELATIONS
======================

.. index:: HARD_SPIN_CORRELATIONS

Spin correlations between the hard scattering process and the
following decay processes are enabled by default. If you want to
disable them, e.g. for spin correlation studies, you can specify the
option :option:`HARD_SPIN_CORRELATIONS: 0`. This is a top-level
setting as opposed to the other ``HARD_DECAYS``-related settings.

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
