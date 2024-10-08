.. _Input structure:
.. _PATH:
.. _RUNDATA:

###############
Input structure
###############


A Sherpa setup is steered by various parameters, associated with the
different components of event generation.

.. index:: PATH
.. index:: RUNDATA

These have to be specified in a configuration file which by default is
named :file:`Sherpa.yaml` residing in the current working directory.  If you
want to use a different setup directory for your Sherpa run, you have
to specify it on the command line as :option:`-p \<dir\>` or
``'PATH: <dir>'`` (including the quotes).

To read parameters from a configuration file with a different name,
you may give the file name as a positional argument on the command line
like this: ``Sherpa <file>``. Note that you can also pass more than
one file like this: ``Sherpa <file1> <file2> ...`` In this case, settings
in files to the right take precedence. This can be useful to reduce
duplication in the case that you have several setups that share a common
set of settings.

Note that you can also pass filenames using the legacy syntax
:option:`-f \<file\>` or ``'RUNDATA: [<file1>, <file2>]'``.
However, this is deprecated.
Use positional arguments instead. Mixing this legacy syntax and positional
arguments for specifying configuration files yields undefined behaviour.

Sherpa’s configuration files are written in the `YAML <https://yaml.org>`_
format.
Most settings are just written as the settings’ name followed by its value,
like this:

.. code-block:: yaml

   EVENTS: 100M
   BEAMS: 2212
   BEAM_ENERGIES: 7000
   ...

In other words, they are key-value pairs of the top-level mapping.
For some settings, the value is itself a mapping.
Hence, we get a nested structure, for example:

.. code-block:: yaml

   HARD_DECAYS:
     Enabled: true
     Apply_Branching_Ratios: false

where ``Enabled`` and ``Apply_Branching_Ratios`` are sub-settings of
the top-level ``HARD_DECAYS`` setting.
The hierarchy is denoted by indentation here.
In YAML, this is called block style and relies on proper formatting
(i.e. each element must be on a separate line, and indentation must be consistent).
Alternatively, one can use flow style, using indicators such as braces
instead of whitespace to indicate structure.
For the previous example, the inner mapping can be written with curly braces and commas:

.. code-block:: yaml

   HARD_DECAYS: { Enabled: true, Apply_Branching_Ratios: false }

Other settings are sequences of elements.
An example would be a sequence of two scale variations:

.. code-block:: yaml

   SCALE_VARIATIONS:
   - 0.25
   - 4.00

In block style, each sequential item is prepended with a single dash.
Equivalently, the snippet can be rewritten in flow style using square brackets
(line breaks are completely optional then, and are omitted here):

.. code-block:: yaml

   SCALE_VARIATIONS: [0.25, 4.00]

Each ``SCALE_VARIATIONS`` item can itself be a sequence (to specify different
variations for the factorisation and renormalisation scale).
Block and flow style can be freely mixed in the different levels:

.. code-block:: yaml

   SCALE_VARIATIONS:
   - 0.25
   - [0.25, 1.00]
   - [1.00, 0.25]

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

Note that we have used flow style here,
because block style would require line breaks,
which are difficult to deal with on the command line.

For scalar (i.e. single-valued) settings, you can use a more convenient
syntax on the command line, where the levels are separated with a colon:

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa KEYWORD1:value1 KEYWORD2:value2 ...
   $ <prefix>/bin/Sherpa HARD_DECAYS:Enabled:false

As this syntax needs no space after the colon,
you can normally suppress quotation marks as we did here.
For non-nested scalar settings, there is yet another possibility,
using an equal sign instead of a colon:

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa KEYWORD1=value1 KEYWORD2=value2 ...

All over Sherpa, particles are defined by the particle code proposed
by the PDG. These codes and the particle properties will be listed
during each run with ``OUTPUT: 2`` for the elementary particles and
``OUTPUT: 4`` for the hadrons.  In both cases, antiparticles are
characterized by a minus sign in front of their code, e.g. a mu- has
code ``13``, while a mu+ has ``-13``.

All dimensionful quantities need to be specified in units of GeV and
millimeter. The same units apply to all numbers in the event output
(momenta, vertex positions).  Scattering cross sections are quoted in
pico-barn in the output.

There are a few extra features for an easier handling of the parameter
file(s), namely algebra interpretation, see `Interpreter`_,
and global tag replacement, see `Tags`_.
In addition, Sherpa compiles a report after each run which lists
all settings for the run, see `Settings Report`_.
Here, the user can review the value for each setting,
whether the setting has actually been used,
and whether it has been customised in some way.

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
occurrence of ``$(<Tag>)`` in all files during read-in. An example
tag definition could read

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa 'TAGS: {QCUT: 20, NJET: 3}'

and then be used in the configuration file like:

.. code-block:: yaml

   RESULT_DIRECTORY: Result_$(QCUT)
   PROCESSES:
   - 93 93 -> 11 -11 93{$(NJET)}:
       Order: {QCD: 0, EW: 2}
       CKKW: $(QCUT)

.. _Settings Report:

***************
Settings Report
***************

Sherpa monitors the use and customization of user settings
during its execution.

This can in particular be useful to alert the user
that some setting has not actually been used at all,
which can be unexpected and perhaps due to a typo on the command line
or in an input file.

Therefore, Sherpa prints a list of top-level settings at the end of each run,
which are unused themselves or contain a subsetting that has not been used.
Thus, the user can decide to inspect any of these settings further,
if this is indeed unexpected.

For further inspection, Sherpa writes the directory ``Settings_Report``
to the working directory towards the end of a run.
It includes a markdown file that contains up to three sections.
The top section lists all settings that have been customized by the user,
but which have not been used by Sherpa during its execution.
If there is no such setting, this section is omitted.
As explained above, this can be consulted to ensure that a given setting
has indeed been used by Sherpa.
If a setting is unexpectedly listed as unused,
it might be due to a typo in the setting name.

The middle section lists all user-customized settings
along with their default values,
the (perhaps multiple) customizations across all configuration files
and the command line, as well as the final value used by Sherpa.
This section can be consulted to double-check
which setting value from what input method
has eventually been used by Sherpa.

The last sections lists all settings that have been used by Sherpa,
but which have not been customized by the user,
along with their default values.
It can be used to understand what other (uncustomized) settings
are relevant for a given Sherpa run,
and to discover further customization points.

The ``Settings_Report`` directory also contains a ``Makefile``
to compile a HTML page from the markdown file,
by simply executing `make` within the directory.
This requires the tool ``pandoc`` to be installed and available
in the ``PATH`` environment variable.
Alternatively, any markdown conversion tool can be used instead.
The compiled HTML page, ``Settings_Report.html``,
can be opened in any standard web browser.
