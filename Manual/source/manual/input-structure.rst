.. _Input structure:

###############
Input structure
###############


A Sherpa setup is steered by various parameters, associated with the
different components of event generation.

.. index:: PATH
.. index:: RUNDATA

These have to be specified in a configuration file which by default is
named :file:`Sherpa.yaml` in the current working directory.  If you
want to use a different setup directory for your Sherpa run, you have
to specify it on the command line as :option:`-p <dir>` or ``'PATH:
<dir>'`` (including the quotes).  To read parameters from a
configuration file with a different name, you may specify :option:`-f
<file>` or ``'RUNDATA: <file>'``.

Sherpa's configuration files are writtin in the `YAML <yaml.org>`_ format.
Most settings are just written as the settings' name followed by its value,
like this:

.. code-block:: yaml

   EVENTS: 100M
   BEAMS: 2212
   BEAM_ENERGIES: 7000
   ...

Others use a more nested structure:

.. code-block:: yaml

   HARD_DECAYS:
       Enabled: true
       Apply_Branching_Ratios: false

where ``Enabled`` and ``Apply_Branching_Ratios`` are sub-settings of
the top-level ``HARD_DECAYS`` setting, which is denoted by indentation
(here two additional spaces).

The different settings and their structure are described in detail in
another chapter of this manual, see :ref:`Parameters`.

All parameters can be overwritten on the command line, i.e.
command-line input has the highest priority.
Each argument is parsed as a single YAML line. This usually means that you have
to quote each argument:

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa 'KEYWORD1: value1' 'KEYWORD2: value2' ...

Because each argument is parsed as YAML, you can also specify nested settings,
e.g. to disable hard decays (even if it is enabled in the config file) you can
write:

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa 'HARD_DECAYS: {Enabled: false}'

Or you can specify the list of matrix-element generators writing:

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa 'ME_GENERATORS: [Comix, Amegic]'

All over Sherpa, particles are defined by the particle code proposed
by the PDG. These codes and the particle properties will be listed
during each run with ``OUTPUT: 2`` for the elementary particles and
``OUTPUT: 4`` for the hadrons.  In both cases, antiparticles are
characterized by a minus sign in front of their code, e.g. a mu- has
code ``13``, while a mu+ has ``-13``.

All quantities have to be specified in units of GeV and
millimeter. The same units apply to all numbers in the event output
(momenta, vertex positions).  Scattering cross sections are denoted in
pico-barn in the output.

There are a few extra features for an easier handling of the parameter
file(s), namely global tag replacement, see `Tags`_, and algebra
interpretation, see `Interpreter`_.


.. contents::
   :local:

.. _Interpreter:

***********
Interpreter
***********

Sherpa has a built-in interpreter for algebraic expressions, like
``cos(5/180*M_PI)``.  This interpreter is employed when reading
integer and floating point numbers from input files, such that certain
parameters can be written in a more convenient fashion.  For example
it is possible to specify the factorisation scale as ``sqr(91.188)``.

There are predefined tags to alleviate the handling

``M_PI``
  Ludolph's Number to a precision of 12 digits.

``M_C``
  The speed of light in the vacuum.

``E_CMS``
  The total centre of mass energy of the collision.

The expression syntax is in general C-like, except for the extra
function ``sqr``, which gives the square of its argument. Operator
precedence is the same as in C.  The interpreter can handle functions
with an arbitrary list of parameters, such as ``min`` and ``max``.

The interpreter can be employed to construct arbitrary variables from
four momenta, like e.g. in the context of a parton level selector, see
:ref:`Selectors`.  The corresponding functions are

:samp:`Mass({v})`
  The invariant mass of :samp:`{v}` in GeV.

:samp:`Abs2({v})`
  The invariant mass squared of :samp:`{v}` in GeV^2.

:samp:`PPerp({v})`
  The transverse momentum of :samp:`{v}` in GeV.

:samp:`PPerp2({v})`
  The transverse momentum squared of :samp:`{v}` in GeV^2.

:samp:`MPerp({v})`
  The transverse mass of :samp:`{v}` in GeV.

:samp:`MPerp2({v})`
  The transverse mass squared of :samp:`{v}` in GeV^2.

:samp:`Theta({v})`
  The polar angle of :samp:`{v}` in radians.

:samp:`Eta({v})`
  The pseudorapidity of :samp:`{v}`.

:samp:`Y({v})`
  The rapidity of :samp:`{v}`.

:samp:`Phi({v})`
  The azimuthal angle of :samp:`{v}` in radians.

:samp:`Comp({v},{i})` The :samp:`{i}`'th component of the vector
  :samp:`{v}`. :samp:`{i}` = 0 is the energy/time component,
  :samp:`{i}` = 1, 2, and 3 are the x, y, and z components.

:samp:`PPerpR({v1},{v2})`
  The relative transverse momentum between :samp:`{v1}` and :samp:`{v2}` in GeV.

:samp:`ThetaR({v1},{v2})`
  The relative angle between :samp:`{v1}` and :samp:`{v2}` in radians.

:samp:`DEta({v1},{v2})`
  The pseudo-rapidity difference between :samp:`{v1}` and :samp:`{v2}`.

:samp:`DY({v1},{v2})`
  The rapidity difference between :samp:`{v1}` and :samp:`{v2}`.

:samp:`DPhi({v1},{v2})`
  The relative polar angle between :samp:`{v1}` and :samp:`{v2}` in radians.

.. _Tags:

****
Tags
****

Tag replacement in Sherpa is performed through the data reading
routines, which means that it can be performed for virtually all
inputs.  Specifying a tag on the command line or in the configuration
file using the syntax ``TAGS: {<Tag>: <Value>}`` will replace every
occurrence of ``<Tag>`` in all files during read-in. An example
tag definition could read

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa 'TAGS: {QCUT: 20, NJET: 3}'

and then be used in the configuration file like:

.. code-block:: yaml

   RESULT_DIRECTORY: Result_$(QCUT)
   PROCESSES:
   - 93 93 -> 11 -11 93{$(NJET)}:
       Order: {QCD: Any, EW: 2}
       CKKW: $(QCUT)
