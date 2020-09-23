.. _Processes:

*********
Processes
*********

The process setup takes the following general form

.. code-block:: yaml

   PROCESSES:
   # Process 1:
   - <process declaration>:
       <parameter>: <value>
       <multiplicities-to-be-applied-to>:
         <parameter>: <value>
         ...
   # Process 2:
   - <process declaration>
       ...

I.e. ``PROCESSES`` followed by a list of process definitions.  Each
definition consists of a key-value mapping with a single key. The key
defines the initial- and final-state particles, and its value gives the
optional parameters for this process steering its setup.
Most of these parameters can be grouped under a multiplicity key, which can
either be a single or a range of multiplicities, e.g. ``2->2-4: {
<settings-that-affect-only-2->2,2->3 and 2->4 processes> }``.

The following parameters are used to steer the process setup:

.. contents::
   :local:
   :depth: 1

.. _Process:

Process
=======

Each process definition starts with the specification of the
(core) process itself. The initial and final state particles are
specified by their PDG codes, or by particle containers, see
:ref:`Particle containers`. Examples are

``- 93 93 -> 11 -11``
  Sets up a Drell-Yan process group with light quarks
  in the initial state.

``- 11 -11 -> 93 93 93{3}``
  Sets up jet production in e+e- collisions with up to three
  additional jets.

The syntax for specifying processes is explained in the following
sections:

.. contents::
   :local:

.. _PDG codes:

PDG codes
---------

Initial and final state particles are specified using their PDG codes
(cf. `PDG
<http://pdg.lbl.gov/2009/mcdata/mc_particle_id_contents.html>`_).  A
list of particles with their codes, and some of their properties, is
printed at the start of each Sherpa run, when the :ref:`OUTPUT` is set
at level :option:`2`.

.. _Particle containers:

Particle containers
-------------------

Sherpa contains a set of containers that collect particles with
similar properties, namely

* lepton (carrying number ``90``),

* neutrino (carrying number ``91``),

* fermion (carrying number ``92``),

* jet (carrying number ``93``),

* quark (carrying number ``94``).


These containers hold all massless particles and anti-particles of the
denoted type and allow for a more efficient definition of initial and
final states to be considered. The jet container consists of the gluon
and all massless quarks, as set by

.. code-block:: yaml

   PARTICLE_DATA:
     <id>:
       Mass: 0
       # ... and/or ...
       Massive: false

A list of particle containers is printed at the start of each Sherpa
run, when the :ref:`OUTPUT` is set at level :option:`2`.

.. index:: PARTICLE_CONTAINER

It is also possible to define a custom particle container using the
keyword ``PARTICLE_CONTAINER``. The container must be given an
unassigned particle ID (kf-code) and its name (freely chosen by you)
and the flavour content must be specified.  An example would be the
collection of all down-type quarks using the unassigned ID 98, which
could be declared as

.. code-block:: yaml

   PARTICLE_CONTAINER:
     98:
       Name: downs
       Flavours: [1, -1, 3, -3, 5, -5]

Note that, if wanted, you have to add both particles and
anti-particles.

.. _Parentheses:

Parentheses
-----------

The parenthesis notation allows to group a list of processes with
different flavor content but similar structure. This is most useful in
the context of simulations containing heavy quarks.  In a setup with
massive b-quarks, for example, the b-quark will not be part of the
jets container. In order to include b-associated processes easily, the
following can be used:

.. code-block:: yaml

   PARTICLE_DATA:
     5: {Massive: true}
   PARTICLE_CONTAINER:
     98: {Name: B, Flavours: [5, -5]}
   PROCESSES:
   - 11 -11 -> (93,98) (93,98):
     ...

.. _Curly brackets:

Curly brackets
--------------

The curly bracket notation when specifying a process allows up to a
certain number of jets to be included in the final state. This is
easily seen from an example, ``11 -11 -> 93 93 93{3}`` sets
up jet production in e+e- collisions. The matix element final state
may be 2, 3, 4 or 5 light partons or gluons.

.. _Decay:

Decay
=====

Specifies the exclusive decay of a particle produced in the matrix
element. The virtuality of the decaying particle is sampled according
to a Breit-Wigner distribution. In practice this amouts to selecting
only those diagrams containing s-channels of the specified flavour
while the phase space is kept general. Consequently, all spin
correlations are preserved.  An example would be

.. code-block:: yaml

    - 11 -11 -> 6[a] -6[b]:
       Decay:
       - 6[a] -> 5 24[c]
       - -6[b] -> -5 -24[d]
       - 24[c] -> -13 14
       - -24[d] -> 94 94


.. _DecayOS:

DecayOS
=======

Specifies the exclusive decay of a particle produced in the matrix
element. The decaying particle is on mass-shell, i.e.  a strict
narrow-width approximation is used. This tag can be specified
alternatively as :option:`DecayOS`. In practice this amouts to
selecting only those diagrams containing s-channels of the specified
flavour and the phase space is factorised as well. Nonetheless, all
spin correlations are preserved.  An example would be

.. code-block:: yaml

   - 11 -11 -> 6[a] -6[b]:
       DecayOS:
       - 6[a] -> 5 24[c]
       - -6[b] -> -5 -24[d]
       - 24[c] -> -13 14
       - -24[d] -> 94 94

.. _No_Decay:

No_Decay
========

Remove all diagrams associated with the decay/s-channel of the given
flavours.  Serves to avoid resonant contributions in processes like
W-associated single-top production. Note that this method breaks gauge
invariance!  At the moment this flag can only be set for Comix.  An
example would be

.. code-block:: yaml

   - 93 93 -> 6[a] -24[b] 93{1}:
       Decay: 6[a] -> 5 24[c]
       DecayOS:
       - 24[c] -> -13 14
       - -24[b] -> 11 -12
       No_Decay: -6

.. _proc_Scales:

Scales
======

Sets a process-specific scale.  For the corresponding syntax see
:ref:`SCALES`.

.. _proc_Couplings:

Couplings
=========

Sets process-specific couplings.  For the corresponding syntax see
:ref:`COUPLINGS`.

.. _CKKW:

CKKW
====

Sets up multijet merging according to :cite:`Hoeche2009rj`.  The
additional argument specifies the parton separation criterion
("merging cut") Q_{cut} in GeV.  It can be given in any form which is
understood by the internal interpreter, see
:ref:`Interpreter`. Examples are


* Hadronic collider: ``CKKW: 20``

* Leptonic collider: ``CKKW: pow(10,-2.5/2.0)*E_CMS``

* DIS: ``CKKW: $(QCUT)/sqrt(1.0+sqr($(QCUT)/$(SDIS))/Abs2(p[2]-p[0]))``

.. _param_Process_Selectors:

Process_Selectors
=================

Using ``Selectors: [<selector 1>, <selector 2>]`` in a process
definition sets up process-specific selectors. They use the same
syntax as describes in :ref:`Selectors`.

.. _Order:

Order
=====

Sets a process-specific coupling order.  Orders are counted at the
amplitude level.  For example, the process 1 -1 -> 2 -2 would have
orders ``{QCD: 2, EW: 0``}, ``{QCD: 1, EW: 1}`` and ``{QCD: 0,
EW: 2}``. There can also be a third entry that is model specific
(e.g. for HEFT couplings). Half-integer orders are so far supported
only by Comix.  The word "Any" can be used as a wildcard.

Note that for decay chains this setting applies to the full process,
see :ref:`Decay` and :ref:`DecayOS`.


.. _Max_Order:

Max_Order
=========

Sets a process-specific maximum coupling order.  See :ref:`Order` for
the syntax and additional information.

.. _Min_Order:

Min_Order
=========

Sets a process-specific minimum coupling order.  See :ref:`Order` for
the syntax and additional information.

.. _Min_N_Quarks:

Min_N_Quarks
============

Limits the minimum number of quarks in the process to the given value.

.. _Max_N_Quarks:

Max_N_Quarks
============

Limits the maximum number of quarks in the process to the given value.

.. _Min_N_TChannels:

Min_N_TChannels
===============

Limits the minimum number of t-channel propagators in the process to
the given value.

.. _Max_N_TChannels:

Max_N_TChannels
===============

Limits the maximum number of t-channel propagators in the process to
the given value.

.. _Print_Graphs:

Print_Graphs
============

Writes out Feynman graphs in LaTeX format. The parameter specifies a
directory name in which the diagram information is stored. This
directory is created automatically by Sherpa. The LaTeX source files
can be compiled using the command

.. code-block:: shell-session

   $ ./plot_graphs <graphs directory>

which creates an html page in the graphs directory that can be viewed
in a web browser.

.. _Name_Suffix:

Name_Suffix
===========

Defines a unique name suffix for the process.

.. _Integration_Error:

Integration_Error
=================

Sets a process-specific relative integration error target.
An example to specify an error target of 2% for
2->3 and 2->4 processes would be:

.. code-block:: yaml

   - 93 93 -> 93 93 93{2}:
       2->3-4:
         Integration_Error: 0.02

.. _Max_Epsilon:

Max_Epsilon
===========

Sets epsilon for maximum weight reduction.  The key idea is to allow
weights larger than the maximum during event generation, as long as
the fraction of the cross section represented by corresponding events
is at most the epsilon factor times the total cross section. In other
words, the relative contribution of overweighted events to the
inclusive cross section is at most epsilon.

.. _Enhance_Factor:

Enhance_Factor
==============

Sets a process specific enhance factor.

.. _RS_Enhance_Factor:

RS_Enhance_Factor
=================

Sets an enhance factor for the RS-piece of an MC\@NLO process.

.. _Enhance_Function:

Enhance_Function
================

Sets a process specific enhance function.

.. note::

   This feature can only be used when generating weighted events.

Note that the convergence of the Monte Carlo integration can be worse
if enhance functions are employed and therefore the integration can
take significantly longer. The reason is that the default phase space
mapping, which is constructed according to diagrammatic information
from hard matrix elements, is not suited for event generation
including enhancement. It must first be adapted, which, depending on
the enhance function and the final state multiplicity, can be an
intricate task.

*If Sherpa cannot achieve an integration error target due to the use
of enhance functions, it might be appropriate to locally redefine this
error target*, see :ref:`Integration_Error`.

.. _Enhance_Observable:

Enhance_Observable
==================


Allows for the specification of a ME-level observable in which the event
generation should be flattened. Of course, this induces an appropriate weight
for each event. This option is available for both weighted and unweighted event
generation, but for the latter as mentioned above the weight stemming from the
enhancement is introduced.

An example would be:

.. code-block:: yaml

   - 93 93 -> 11 -11 93{1}:
       2->3:
         Enhance_Observable: VAR{log10(PPerp(p[2]+p[3]))}|1|3

Here, the 1-jet process is flattened with respect to the logarithmic
transverse momentum of the lepton pair in the limits 1.0 (10 GeV) to
3.0 (1 TeV).  For the calculation of the observable one can use any
function available in the algebra interpreter (see
:ref:`Interpreter`).

Note that the convergence of the Monte Carlo integration can be worse
if enhance observables are employed and therefore the integration can
take significantly longer. The reason is that the default phase space
mapping, which is constructed according to diagrammatic information
from hard matrix elements, is not suited for event generation
including enhancement. It must first be adapted, which, depending on
the enhance function and the final state multiplicity, can be an
intricate task.

*If Sherpa cannot achieve an integration error target due to the use
of enhance functions, it might be appropriate to locally redefine this
error target*, see :ref:`Integration_Error`.

.. _NLO_Mode:

NLO_Mode
========

This setting specifies whether and in which mode an NLO calculation
should be performed. Possible values are:

``None``
  perform a leading-order calculation (this is the default)

``Fixed_Order``
  perform a fixed-order next-to-leading order calculation

``MC@NLO``
  perform an MC\@NLO-type matching of a fixed-order next-to-leading order
  calculation to the resummation of the parton shower

The usual multiplicity identifier applies to this switch as well.
Note that using a value other than ``None`` implies ``NLO_Part: BVIRS`` for
the relevant multiplicities.
For fixed-order NLO calculations (``NLO_Mode: Fixed_Order``), this can be
overridden by setting ``NLO_Part`` explicitly, see :ref:`NLO_Part`.

Note that Sherpa includes only a very limited selection of one-loop
corrections. For processes not included external codes can be
interfaced, see :ref:`External one-loop ME`

.. _NLO_Part:

NLO_Part
========

In case of fixed-order NLO calculations this switch specifies which
pieces of a NLO calculation are computed, also see :ref:`NLO_Mode`.
Possible choices are

``B``
  born term

``V``
  virtual (one-loop) correction

``I``
  integrated subtraction terms

``RS``
  real correction, regularized using Catani-Seymour subtraction terms

Different pieces can be combined in one processes setup. Only pieces
with the same number of final state particles and the same order in
alpha_S and alpha can be treated as one process, otherwise they will
be automatically split up.

.. _NLO_Order:

NLO_Order
=========

Specifies the relative order of the NLO correction wrt. the considered
Born process. For example, ``NLO_Order: {QCD: 1, EW: 0}`` specifies
a QCD correction while ``NLO_Order: {QCD: 0, EW: 1}`` specifies an
EW correction.

.. _Subdivide_Virtual:

Subdivide_Virtual
=================

Allows to split the virtual contribution to the total cross section
into pieces.  Currently supported options when run with
`BlackHat <https://projects.hepforge.org/blackhat>`_ are
:option:`LeadingColor` and :option:`FullMinusLeadingColor`. For
high-multiplicity calculations these settings allow to adjust the
relative number of points in the sampling to reduce the overall
computation time.

.. _ME_Generator:

ME_Generator
============

Set a process specific nametag for the desired tree-ME generator, see
:ref:`ME_GENERATORS`.

.. _RS_ME_Generator:

RS_ME_Generator
===============

Set a process specific nametag for the desired ME generator used for
the real minus subtraction part of NLO calculations. See also
:ref:`ME_GENERATORS`.

.. _Loop_Generator:

Loop_Generator
==============

Set a process specific nametag for the desired loop-ME generator. The
only Sherpa-native option is ``Internal`` with a few hard coded loop
matrix elements. Other loop matrix elements are provided by external
libraries.

.. _Integrator:

Integrator
==========

Sets a process-specific integrator, see :ref:`int_INTEGRATOR`.

.. _PSI_ItMin:

PSI_ItMin
=========

Sets the number of points per optimization step, see :ref:`PSI`.

.. _RS_PSI_ItMin:

RS_PSI_ItMin
============

Sets the number of points per optimization step in real-minus-subtraction
parts of fixed-order and MC\@NLO calculations, see :ref:`PSI`.

.. _Special Group:

Special Group
=============

Allows to split up individual flavour processes within a process group for
integrating them separately. This can help improve the integration/unweighting
efficiency. Note: Only works with Comix so far.
Example for usage:

.. code-block:: yaml

   Process 93 93 -> 11 -11 93
   Special Group(0-1,4)
   [...]
   End process
   Process 93 93 -> 11 -11 93
   Special Group(2-3,5-7)
   [...]
   End process

The numbers for each individual process can be found using a script in
the AddOns directory: :file:`AddOns/ShowProcessIds.sh Process/Comix.zip`

.. _End process:

End process
===========

Completes the setup of a process or a list of processes with common
properties.
