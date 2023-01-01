.. _Selectors:

*********
Selectors
*********

.. index:: SHOW_SELECTOR_SYNTAX

Sherpa provides the following selectors that set up cuts at the matrix
element level:

.. contents::
   :local:

Some selectors modify the momenta and flavours of the set of final
state particles. These selectors also take a list of subselectors
which then act on the modified flavours and momenta.  Details are
explained in the respective selectors' description.

.. _Inclusive selectors:

Inclusive selectors
===================

The selectors listed here implement cuts on the matrix element level,
based on event properties.  The corresponding syntax is

.. code-block:: yaml

   SELECTORS:
     - [<keyword>, <parameter 1>, <parameter 2>, ...]
     - # other selectors ...

Parameters that accept numbers can also be given in a form that is
understood by the internal algebra interpreter, see
:ref:`Interpreter`.  The selectors act on *all* particles in the
event.  Their respective keywords are

:option:`[N, <kf>, <min value>, <max value>]`
  Minimum and maximum multiplicity of flavour ``<kf>`` in the
  final state.

:option:`[PTmis, <min value>, <max value>]`
  Missing transverse momentum cut (at the moment only neutrinos are
  considered invisible)

:option:`[ETmis, <min value>, <max value>]`
  Missing transverse energy cut (at the moment only neutrinos are
  considered invisible)

:option:`[IsolationCut, <kf>, <dR>, <exponent>, <epsilon>, <optional: mass_max>]`
  Smooth cone isolation :cite:`Frixione1998jh`, the parameters given
  are the flavour ``<kf>`` to be isolated against massless partons
  and the isolation cone parameters.

:option:`[NJ, <N>, <algo>, <min value>, <max value>]`
  NJettiness from :cite:`Stewart2010tn`, where ``<algo>`` specifies
  the jet finding algorithm to determine the hard jet directions and
  ``<N>`` is their multiplicity.
  ``algo=kt|antikt|cambridge|siscone,PT:<ptmin>,R:<dR>[,[ETA:<etamax>,Y:<ymax>]]``

.. _One particle selectors:

One particle selectors
======================

The selectors listed here implement cuts on the
matrix element level, based on single particle kinematics.
The corresponding  syntax is

.. code-block:: yaml

   SELECTORS:
     - [<keyword>, <flavour code>, <min value>, <max value>]
     - # other selectors ...

:option:`<min value>` and :option:`<max value>` are floating point
numbers, which can also be given in a form that is understood by the
internal algebra interpreter, see :ref:`Interpreter`.  The selectors
act on *all* possible particles with the given flavour. Their
respective keywords are

:option:`PT`
  transverse momentum cut

:option:`ET`
  transverse energy cut

:option:`Y`
  rapidity cut

:option:`Eta`
  pseudorapidity cut

:option:`PZIN`
  cut on the z-component of the momentum, acts on initial-state flavours only
  (commonly used in DIS analyses)

:option:`HT`
  Visible transverse energy cut


.. _Two particle selectors:

Two particle selectors
======================

The selectors listed here implement cuts on the matrix element level,
based on two particle kinematics.  The corresponding is

.. code-block:: yaml

   SELECTORS:
     - [<keyword>, <flavour1 code>, <flavour2 code>, <min value>, <max value>]
     - # other selectors ...

:option:`<min value>` and :option:`<max value>` are floating point
numbers, which can also be given in a form that is understood by the
internal algebra interpreter, see :ref:`Interpreter`.  The selectors
act on *all* possible particles with the given flavour. Their
respective keywords are

:option:`Mass`
  invariant mass

:option:`Q2`
  DIS-like virtuality

:option:`PT2`
  pair transverse momentum

:option:`MT2`
  pair transverse mass

:option:`DY`
  rapidity separation

:option:`DEta`
  pseudorapidity separation

:option:`DPhi`
  azimuthal separation

:option:`DR`
  angular separation (build from eta and phi)

:option:`DR(y)`
  angular separation (build from y and phi)

:option:`INEL`
  inelasticity, one of the flavours must be in the initial-state (commonly used
  in DIS analyses)

.. _Decay selectors:

Decay selectors
===============

The selectors listed here implement cuts on the matrix element level,
based on particle decays, see :ref:`Decay` and :ref:`DecayOS`.

:option:`DecayMass`
  Invariant mass of a decaying particle. The syntax is

  .. code-block:: yaml

     - [DecayMass, <flavour code>, <min value>, <max value>]

:option:`Decay`
  Any kinematic variable of a decaying particle. The syntax is

  .. code-block:: yaml

     - [Decay(<expression>), <flavour code>, <min value>, <max value>]

  where ``<expression>`` is an expression handled by the internal
  interpreter, see :ref:`Interpreter`.

:option:`Decay2`
  Any kinematic variable of a pair of decaying particles. The syntax is

  .. code-block:: yaml

     - [Decay2(<expression>), <flavour1 code>, <flavour2 code>, <min value>, <max value>]

  where ``<expression>`` is an expression handled by the internal
  interpreter, see :ref:`Interpreter`.

Particles are identified by flavour, i.e. the cut is applied on
*all* decaying particles that match :option:`<flavour code>`.
:option:`<min value>` and :option:`<max value>` are floating point
numbers, which can also be given in a format that is understood by the
internal algebra interpreter, see :ref:`Interpreter`.

.. _Particle dressers:

Particle dressers
=================

.. _Jet selectors:

Jet selectors
=============

There are two different types of jet finders

:option:`NJetFinder`
  k_T-type algorithm to select on a given number of jets

:option:`FastjetFinder`
  Select on a given number of jets using FastJet algorithms

Their respective syntax and defaults are

.. code-block:: yaml

   SELECTORS:
   - NJetFinder:
       N: 0
       PTMin: 0.0
       ETMin: 0.0
       R: 0.4
       Exp: 1
       EtaMax: None
       YMax: None
       MassMax: 0.0
   - FastjetFinder:
       Algorithm: kt
       N: 0
       PTMin: 0.0
       ETMin: 0.0
       DR: 0.4
       f: 0.75        # Siscone f parameter
       EtaMax: None
       YMax: None
       Nb: -1
       Nb2: -1

Note that all parameters are optional. If they are not specified,
their respective default values as indicated in the above snippet is
used.  However, at the very least the number of jets, :option:`N`,
should be specified to require a non-zero number of jets.

The :option:`NJetFinder` allows to select for kinematic configurations
with at least :option:`<N>` jets that satisfy both, the
:option:`<PTMin>` and the :option:`<ETMin>` minimum requirements and
that are in a pseudo-rapidity region :option:`|eta|<EtaMax>`. The
:option:`<Exp>` (exponent) allows to apply a kt-algorithm
(1) or an anti-kt algorithm (-1). As only massless partons are clustered by
default, the :option:`<MassMax>` allows to also include partons with a mass
up to the specified values. This is useful e.g. in calculations with massive
b-quarks which shall nonetheless satisfy jet criteria.

The second option :option:`FastjetFinder` allows to use the `FastJet
<http://www.fastjet.fr>`_ plugin, through fjcore.
It takes the following arguments: ``<Algorithm>`` can take the values
``kt,antikt,cambridge,siscone,eecambridge,jade``, ``<N>`` is the
minimum number of jets to be found, ``<PTMin>`` and ``<ETMin>`` are
the minimum transverse momentum and/or energy, ``<DR>`` is the radial
parameter.  Optional arguments are: ``<f>`` (default 0.75, only
relevant for the Siscone algorithm), ``<EtaMax>`` and ``<YMax>`` as
maximal absolute (pseudo-)rapidity, ``<Nb>`` and ``<Nb2>`` set the
number of required b-jets, where for the former both b and anti-b
quarks are counted equally towards b-jets, while for the latter they
are added with a relative sign as constituents, i.e. a jet containing
b and anti-b is not tagged (default: -1, i.e. no b jets are required).
Note that only ``<Algorithm>``, ``<N>`` and ``<PTMin>`` are relevant
for the lepton-lepton collider algorithms.

The selector :option:`FastjetVeto` allows to use the `FastJet
<http://www.fastjet.fr>`_ plugin to apply jet veto cuts. Its syntax is
identical to :option:`FastjetFinder`.

The momenta and
nodal values of the jets found with FastJet can also be used to
calculate more elaborate selector criteria. The syntax of this
selector is

.. code-block:: yaml

   - FastjetSelector:
       Expression: <expression>
       Algorithm: kt
       N: 0
       PTMin: 0.0
       ETMin: 0.0
       DR: 0.4
       f: 0.75
       EtaMax: None
       YMax: None
       BMode: 0

wherein ``Algorithm`` can take the values
``kt,antikt,cambridge,siscone,eecambridge,jade``.  In the algebraic
``<expression>``, ``MU_n2`` (n=2..njet+1) signify the nodal values of
the jets found and ``p[i]`` are their momenta. For details see
:ref:`Scale setters`. For example, in lepton pair production in
association with jets

.. code-block:: yaml

   - FastjetSelector:
       Expression: Mass(p[4]+p[5])>100
       Algorithm: antikt
       N: 2
       PTMin: 40
       ETMin: 0
       DR: 0.5

selects all phase space points where two anti-kt jets with at least 40
GeV of transverse momentum and an invariant mass of at least 100 GeV
are found. The expression must calculate a boolean value.  The
``BMode`` parameter, if specified different from its default 0, allows
to use b-tagged jets only, based on the parton-level constituents of
the jets.  There are two options: With ``BMode: 1`` both b and anti-b
quarks are counted equally towards b-jets, while for ``BMode: 2`` they
are added with a relative sign as constituents, i.e. a jet containing
b and anti-b is not tagged.  Note that only ``<epression>``,
``<algorithm>``, ``<n>`` and ``<ptmin>`` are relevant when using the
lepton-lepton collider algorithms.

.. _Isolation selector:

Isolation selector
==================

Instead of the simple ``IsolationCut`` (:ref:`Inclusive selectors`), you may
also use the more flexible ``Isolation_Selector`` to require photons (or other
particles) with a smooth cone isolation and additionally apply further criteria
to them. Example:

.. code-block:: yaml

   SELECTORS:
   - Isolation_Selector:
       Isolation_Particle: 22
       Rejection_Particles: [93]
       Isolation_Parameters:
         R: 0.1
         EMAX: 0.1
         EXP: 2
         PT: 0.
         Y: 2.7
       NMin: 2
       Remove_Nonisolated: true
       Subselectors:
       - VariableSelector:
           Variable: PT
           Flavs: [22]
           Ranges:
           - [20, E_CMS]
           - [18, E_CMS]
           Ordering: [PT_UP]
       - [DR, 22, 22, 0.2, 10000.0 ]
       #for integration efficiency: m_yy >= sqrt(2 pTmin1 pTmin2 (1-cos dR))
       - [Mass, 22, 22, 3.7, E_CMS]

.. _Universal selector:

Universal selector
==================

.. index:: SHOW_VARIABLE_SYNTAX

The universal selector is intended to implement non-standard cuts on
the matrix element level. Its syntax is

.. code-block:: yaml

   SELECTORS:
   - VariableSelector:
       Variable: <variable>
       Flavs: [<kf1>, ..., <kfn>]
       Ranges:
       - [<min1>, <max1>]
       - ...
       - [<minn>, <maxn>]
       Ordering: [<order1>, ..., <orderm>]

The ``Variable`` parameter defines the name of the variable to cut
on. The keywords for available predefined can be figured by running
Sherpa :option:`SHOW_VARIABLE_SYNTAX: true`.  Or alternatively, an
arbitrary cut variable can be constructed using the internal
interpreter, see :ref:`Interpreter`. This is invoked with the command
``Calc(...)``. In the formula specified there you have to use place
holders for the momenta of the particles: ``p[0]`` ... ``p[n]`` hold
the momenta of the respective particles ``kf1`` ... ``kfn``. A list of
available vector functions and operators can be found here
:ref:`Interpreter`.

:option:`<kf1>,..,<kfn>` specify the PDG codes of the particles the
variable has to be calculated from.  The ranges :option:`[<min>,
<max>]` then define the cut regions.

If the ``Ordering`` parameter is not given, the order of cuts is
determined internally, according to Sherpa's process classification
scheme.  This then has to be matched if you want to have different
cuts on certain different particles in the matrix element. To do this,
you should put enough (for the possible number of combinations of your
particles) arbitrary ranges at first and run Sherpa with debugging
output for the universal selector: ``Sherpa 'FUNCTION_OUTPUT:
{"Variable_Selector::Trigger": 15}'``.  This will start to produce
lots of output during integration, at which point you can interrupt
the run (Ctrl-c). In the ``Variable_Selector::Trigger(): {...}``
output you can see, which particle combinations have been found and
which cut range your selector has held for them (vs. the arbitrary
range you specified). From that you should get an idea, in which order
the cuts have to be specified.

If the fourth argument is given, particles are ordered before the cut
is applied. Possible orderings are :option:`PT_UP`, :option:`ET_UP`,
:option:`E_UP`, :option:`ETA_UP` and :option:`ETA_DOWN`, (increasing
p_T, E_T, E, eta, and decreasing eta). They have to be specified for
each of the particles, separated by commas.

Examples

.. code-block:: yaml

   SELECTORS:
   # two-body transverse mass
   - VariableSelector:
       Variable: mT
       Flavs: [11, -12]
       Ranges:
       - [50, E_CMS]

   # cut on the pT of only the hardest lepton in the event
   - VariableSelector:
       Variable: PT
       Flavs: 90
       Ranges:
       - [50, E_CMS]
       Ordering: [PT_UP]

   # using bool operations to restrict eta of the electron to |eta| < 1.1 or
   # 1.5 < |eta| < 2.5
   - VariableSelector:
       Variable: Calc(abs(Eta(p[0]))<1.1||(abs(Eta(p[0]))>1.5&&abs(Eta(p[0]))<2.5))
       Flavs: 11
       Ranges:
       - [1, 1]  # NOTE: this means true for bool operations

   # requesting opposite side tag jets in VBF
   - VariableSelector:
       Variable: Calc(Eta(p[0])*Eta(p[1]))
       Flavs: [93, 93]
       Ranges:
       - [-100, 0]
       Ordering: [PT_UP, PT_UP]

   # restricting electron+photon mass to be outside of [87.0,97.0]
   - VariableSelector:
       Variable: Calc(Mass(p[0]+p[1])<87.0||Mass(p[0]+p[1])>97.0)
       Flavs: [11, 22]
       Ranges:
       - [1, 1]

   # in ``Z[lepton lepton] Z[lepton lepton]``, cut on mass of lepton-pairs
   # produced from Z's
   - VariableSelector:
       Variable: m
       Flavs: [90, 90]
       # here we use knowledge about the internal ordering to cut only on the
       # correct lepton pairs
       Ranges:
       - [80, 100]
       - [0, E_CMS]
       - [0, E_CMS]
       - [0, E_CMS]
       - [0, E_CMS]
       - [80, 100]


.. _Minimum selector:

Minimum selector
================

This selector can combine several selectors to pass an event if either
those passes the event. It is mainly designed to generate more
inclusive samples that, for instance, include several jet finders and
that allows a specification later. The syntax is

.. code-block:: yaml

   SELECTORS:
   - MinSelector:
       Subselectors:
       - <selector 1>
       - <selector 2>
       ...

The :ref:`Minimum selector` can be used if constructed with other
selectors mentioned in this section
