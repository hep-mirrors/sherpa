.. _Shower Parameters:

**************
Parton showers
**************

The following parameters are used to steer the shower setup.

.. contents::
   :local:

.. _SHOWER_GENERATOR:

SHOWER_GENERATOR
================

There are two shower generators in Sherpa, :option:`Dire` (default)
and :option:`CSS`. See the module summaries in :ref:`Basic structure`
for details about these showers.

Other shower modules are in principle supported and more choices will
be provided by Sherpa in the near future.  To list all available
shower modules, the tag ``SHOW_SHOWER_GENERATORS: 1`` can be specified
on the command line.

``SHOWER_GENERATOR: None`` switches parton showering off completely.
However, even in the case of strict fixed order calculations, this
might not be the desired behaviour as, for example, then neither the
METS scale setter, cf. :ref:`SCALES`, nor Sudakov rejection weights
can be employed.  To circumvent when using the Dire or CS Shower see
:ref:`Sherpa Shower options`.

.. _JET_CRITERION:

JET_CRITERION
=============

This option uses the value for ``SHOWER_GENERATOR`` as its default.
Correspondingly, the only natively supported options in Sherpa are
:option:`CSS` and :option:`Dire`. The corresponding jet criterion is
described in :cite:`Hoeche2009rj`. A custom jet criterion, tailored to
a specific experimental analysis, can be supplied using Sherpa's
plugin mechanism.

.. _MASSIVE_PS:

MASSIVE_PS
==========

This option instructs Sherpa to treat certain partons as massive in
the shower, which have been considered massless by the matrix element.
The argument is a list of parton flavours, for example
:option:`MASSIVE_PS: [4, 5]`, if both c- and b-quarks are to be
treated as massive.

.. _MASSLESS_PS:

MASSLESS_PS
===========

When hard decays are used, Sherpa treats all flavours as massive in
the parton shower. This option instructs Sherpa to treat certain
partons as massless in the shower nonetheless. The argument is a list
of parton flavours, for example :option:`MASSLESS_PS: [1, 2, 3]`, if
u-, d- and s-quarks are to be treated as massless.

.. _Sherpa Shower options:

Sherpa Shower options
=====================

.. index:: KIN_SCHEME
.. index:: IS_PT2MIN
.. index:: FS_PT2MIN
.. index:: IS_AS_FAC
.. index:: FS_AS_FAC
.. index:: MAXEM

Sherpa's default shower module is based on :cite:`Schumann2007mg`.  A
new ordering parameter for initial state splitters was introduced in
:cite:`Hoeche2009rj` and a novel recoil strategy for initial state
splittings was proposed in :cite:`Hoeche2009xc`.  While the ordering
variable is fixed, the recoil strategy for dipoles with initial-state
emitter and final-state spectator can be changed for systematics
studies. Setting :option:`SHOWER:KIN_SCHEME: 0` corresponds to using the
recoil scheme proposed in :cite:`Hoeche2009xc`, while
:option:`SHOWER:KIN_SCHEME: 1` (default) enables the original recoil
strategy.  The lower cutoff of the shower evolution can be set via
:option:`SHOWER:FS_PT2MIN` and :option:`SHOWER:IS_PT2MIN` for final and
initial state shower, respectively.  Note that this value is specified
in GeV^2. Scale factors for the evaluation of the strong coupling in
the parton shower are given by :option:`SHOWER:FS_AS_FAC` and
:option:`SHOWER:IS_AS_FAC`. They multiply the ordering parameter, which
is given in units of GeV^2.

Setting :option:`SHOWER:MAXEM: <N>` forces the CS Shower to truncate its
evolution at the Nth emission. Note that in this case not all of the
Sudakov weights might be computed correctly. On the other hand, the
use of CS Shower in the METS scale setter is not affected,
cf. :ref:`SCALES`.

The parton shower coupling scales, PDF scales and PDF themselves
can be varied on-the-fly, along with the on-the-fly variations
of the corresponding matrix element parameters.
See :ref:`On-the-fly event weight variations`
to find out how specify the variations
and enable them in the shower.

Most parton showers available in Sherpa allow the same options.
These options are specified as follows:

.. code-block:: yaml

   SHOWER:
     KIN_SCHEME: <scheme>
     IS_AS_FAC: <factor>
     # other shower settings ...

When the parton shower is used for MC@NLO matching, the options
can be set differently. They are then specified as follows:

.. code-block:: yaml

   MC@NLO:
     IS_AS_FAC: <factor>
     # other shower settings ...

.. _CS Shower options:

CS Shower options
=================

.. index:: EW_MODE
.. index:: MASS_THRESHOLD
.. index:: EVOLUTION_SCHEME
.. index:: SCALE_SCHEME

By default, only QCD splitting functions are enabled in the CS shower.
If you also want to allow for photon splittings, you can enable them
by using :option:`SHOWER:EW_MODE: true`. Note, that if you have leptons
in your matrix-element final state, they are by default treated by a
soft photon resummation as explained in :ref:`QED Corrections`. To
avoid double counting, this has to be disabled as explained in that
section.

The evolution variable of the CS shower can be changed using
:option:`SHOWER:EVOLUTION_SCHEME`. Several options are currently implemented:

:option:`0`
  transverse momentum ordering

:option:`1`
  modified transverse momentum ordering.

:option:`2`
  like `0` but parton masses taken into account

:option:`3`
  like `1` but parton masses taken into account

:option:`20`
  like `0` but parton masses taken into account only for g->QQ

:option:`30`
  like `1` but parton masses taken into account only for g->QQ

The scale can be set differently for final- and initial-state shower.
The two values are combined as `FS+100*IS`, where `FS` is the choice
for the final state, and `IS` is the choice for the initial state.
The scale at which the strong coupling for shower splittings
is evaluated can be chosen with
:option:`SHOWER:SCALE_SCHEME`:

The default is to evaluate the strong coupling at the transverse momentum
in the parton splitting. Gluon splittings into quarks in the final state
are evaluated at the virtuality of the gluon, as are branchings into
a soft t-channel gluon in the initial state. Options are additive.

:option:`0`
  default

:option:`1`
  evaluate final-state gluon splitting into quarks at the transverse momentum

:option:`2`
  evaluate initial-state quark to gluon splittings at the transverse momentum

:option:`20`
  evaluate initial-state gluon splitting into soft t-channel gluons at the transverse momentum

Additionally,
the CS shower allows to disable splittings at scales below the
on-shell mass of heavy quarks. The upper limit for the corresponding
heavy quark mass is set using :option:`SHOWER:MASS_THRESHOLD`.
