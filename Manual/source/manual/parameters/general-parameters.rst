.. _General Parameters:

******************
General parameters
******************


The following parameters describe general run information.  See
:ref:`Input structure` for how to use them in a configuration file or
on the command line.

.. contents::
   :local:

.. _param_EVENTS:

EVENTS
======

.. index:: EVENTS

This parameter specifies the number of events to be generated.

It can alternatively be set on the command line through option
:option:`-e`, see :ref:`Command line`.

.. _EVENT_TYPE:

EVENT_TYPE
==========

.. index:: EVENT_TYPE

This parameter specifies the kind of events to be generated.  It can
alternatively be set on the command line through option :option:`-t`,
see :ref:`Command line`.

* The default event type is ``StandardPerturbative``, which will
  generate a hard event through exact matrix elements matched and/or
  merged with the paerton shower, eventually including hadronization,
  hadron decays, etc..

Alternatively there are two more specialised modes, namely:

* ``MinimumBias``, which generates minimum bias events through the
  SHRIMPS model implemented in Sherpa, see :ref:`Minimum Bias`

* ``HadronDecay``, which allows to simulate the decays of a specific
  hadron.

.. _SHERPA_VERSION:

SHERPA_VERSION
==============

.. index:: SHERPA_VERSION

This parameter ties a config file to a specific Sherpa version, e.g.
``SHERPA_VERSION: 2.2.0``. If two parameters are given they are
interpreted as a range of Sherpa versions: ``SHERPA_VERSION: [2.2.0,
2.2.5]`` specifies that this config file can be used with any Sherpa
version between (and including) 2.2.0 and 2.2.5.

.. _TUNE:

TUNE
====

.. index:: TUNE

.. warning::

   This parameter is currently not supported.

..
   This parameter specifies which tune is to be used. Setting different
   tunes using this parameter ensures, that consistent settings are
   employed. This affects mostly :ref:`MPI Parameters` and
   :ref:`Intrinsic Transverse Momentum` parameters. Possible values are
   (for Sherpa 2.1.1):

   * ``CT10`` MPI tune for the Sherpa's default PDF, CT10. This is the default tune.

   * ``CT10_UEup`` Upward variation of MPI activity, variation of the CT10 tune to
     assess MPI uncertainties.

   * ``CT10_UEdown`` Downward variation of MPI activity, variation of the CT10 tune to
     assess MPI uncertainties.


.. _OUTPUT:

OUTPUT
======

.. index:: OUTPUT
.. index:: OUTPUT_PRECISION
.. index:: EVT_OUTPUT
.. index:: EVT_OUTPUT_START
.. index:: FUNCTION_OUTPUT

This parameter specifies the screen output level (verbosity) of the
program.  If you are looking for event file output options please
refer to section :ref:`Event output formats`.

It can alternatively be set on the command line through option
:option:`-O`, see :ref:`Command line`. A different output level can be
specified for the event generation step through :option:`EVT_OUTPUT`
or command line option :option:`-o`, see :ref:`Command line`

The value can be any sum of the following:

* 0: Error messages (-> always displayed).
* 1: Event display.
* 2: Informational messages during the run.
* 4: Tracking messages (lots of output).
* 8: Debugging messages (even more output).

E.g. :option:`OUTPUT=3` would display information, events and
errors. Use :option:`OUTPUT_PRECISION` to set the default output
precision (default ``6``).  Note: this may be overriden in specific
functions' output.

For expert users: The output level can be overriden for individual
functions, e.g. like this

.. code-block:: yaml

   FUNCTION_OUTPUT:
     "void SHERPA::Matrix_Element_Handler::BuildProcesses()": 8
     ...

where the function signature is given by the value of
``__PRETTY_FUNCTION__`` in the function block.  Another expert
parameter is :option:`EVT_OUTPUT_START`, with which the first event
affected by :option:`EVT_OUTPUT` can be specified. This can be useful
to generate debugging output only for events affected by a some issue.

.. _LOG_FILE:

LOG_FILE
========

.. index:: LOG_FILE

This parameter specifies the log file. If set, the standard output
from Sherpa is written to the specified file, but output from child
processes is not redirected. This option is particularly useful to
produce clean log files when running the code in MPI mode, see
:ref:`MPI parallelization`.  A file name can alternatively be
specified on the command line through option :option:`-l`, see
:ref:`Command line`.

.. _RANDOM_SEED:

RANDOM_SEED
===========

.. index:: RANDOM_SEED

Sherpa uses different random-number generators. The default is the
Ran3 generator described in :cite:`NumRec2007`.  Alternatively, a
combination of George Marsaglias KISS and SWB :cite:`marsaglia1991`
can be employed, see `this
<http://groups.google.co.uk/group/sci.stat.math/msg/edcb117233979602>`_
`website
<http://groups.google.co.uk/group/sci.math.num-analysis/msg/eb4ddde782b17051>`_.
The integer-valued seeds of the generators are specified by
:option:`RANDOM_SEED: [A, .., D]`. They can also be set individually
using :option:`RANDOM_SEED1: A` through :option:`RANDOM_SEED4: D`. The
Ran3 generator takes only one argument (in this case, you can simply
use :option:`RANDOM_SEED: A`). This value can also be set using the
command line option :option:`-R`, see :ref:`Command line`.

.. _EVENT_SEED_MODE:

EVENT_SEED_MODE
===============

The tag :option:`EVENT_SEED_MODE` can be used to enforce the same
seeds in different runs of the generator. When set to 1, existing
random seed files are read and the seed is set to the next available
value in the file before each event. When set to 2, seed files are
written to disk.  These files are gzip compressed, if Sherpa was
compiled with option :option:`-DSHERPA_ENABLE_GZIP=ON`.  When set to 3, Sherpa
uses an internal bookkeeping mechanism to advance to the next
predefined seed.  No seed files are written out or read in.

.. _ANALYSIS:

ANALYSIS
========

.. index:: ANALYSIS

Analysis routines can be switched on or off using the ANALYSIS
parameter.  The default is no analysis.  This parameter can also be
specified on the command line using option :option:`-a`, see
:ref:`Command line`.

The following analysis handlers are currently available

:option:`Internal`
  | Sherpa's internal analysis handler.
  | To use this option, the package must be configured with option
  | :option:`-DSHERPA_ENABLE_ANALYSIS=ON`. An output directory can
  | be specified using :ref:`ANALYSIS_OUTPUT`.

:option:`Rivet`
  | The Rivet package, see `Rivet Website <http://projects.hepforge.org/rivet/>`_.
  | To enable it, Rivet and HepMC have to be installed and Sherpa must be configured
  | as described in :ref:`Rivet analyses`.

Multiple options can also be specified, e.g. ``ANALYSIS: [Internal,
Rivet]``.

.. _ANALYSIS_OUTPUT:

ANALYSIS_OUTPUT
===============

.. index:: ANALYSIS_OUTPUT

Name of the directory for histogram files when using the internal
analysis and name of the Yoda file when using Rivet, see
:ref:`ANALYSIS`.  The directory/file will be created w.r.t. the
working directory. The default value is ``Analysis/``. This parameter
can also be specified on the command line using option :option:`-A`,
see :ref:`Command line`.

.. _TIMEOUT:

TIMEOUT
=======

.. index:: TIMEOUT

A run time limitation can be given in user CPU seconds through
:option:`TIMEOUT`. This option is of some relevance when running
SHERPA on a batch system. Since in many cases jobs are just
terminated, this allows to interrupt a run, to store all relevant
information and to restart it without any loss. This is particularly
useful when carrying out long integrations.  Alternatively, setting
the :option:`TIMEOUT` variable to -1, which is the default setting,
translates into having no run time limitation at all. The unit is
seconds.

.. _RLIMIT_AS:

RLIMIT_AS
=========

.. index:: RLIMIT_AS
.. index:: RLIMIT_BY_CPU
.. index:: MEMLEAK_WARNING_THRESHOLD

A memory limitation can be given to prevent Sherpa to crash the system
it is running on as it continues to build up matrix elements and loads
additional libraries at run time. Per default the maximum RAM of the
system is determined and set as the memory limit. This can be changed
by giving :option:`RLIMIT_AS: <size>` where the size is given as
e.g. ``500 MB``, ``4 GB``, or ``10 %``.  When running with :ref:`MPI
parallelization` it might be necessary to divide the total maximum by
the number of cores. This can be done by setting ``RLIMIT_BY_CPU:
true``.

Sherpa checks for memory leaks during integration and event
generation.  If the allocated memory after start of integration or
event generation exceeds the parameter
:option:`MEMLEAK_WARNING_THRESHOLD`, a warning is printed.  Like
:option:`RLIMIT_AS`, :option:`MEMLEAK_WARNING_THRESHOLD` can be set
using units.  The warning threshold defaults to ``16MB``.

.. _BATCH_MODE:

BATCH_MODE
==========

.. index:: BATCH_MODE
.. index:: EVENT_DISPLAY_INTERVAL

Whether or not to run Sherpa in batch mode. The default is ``1``,
meaning Sherpa does not attempt to save runtime information when
catching a signal or an exception. On the contrary, if option ``0`` is
used, Sherpa will store potential integration information and analysis
results, once the run is terminated abnormally. All possible settings
are:

:samp:`{0}`
      Sherpa attempts to write out integration and analysis
      results when catching an exception.

:samp:`{1}`
      Sherpa does not attempt to write out integration and
      analysis results when catching an exception.

:samp:`{2}`
      Sherpa outputs the event counter continously, instead of
      overwriting the previous one (default when using
      :ref:`LOG_FILE`).

:samp:`{4}`
      Sherpa increases the on-screen event counter in constant
      steps of 100 instead of an increase relative to the current
      event number. The interval length can be adjusted with
      ``EVENT_DISPLAY_INTERVAL``.

:samp:`{8}`
      Sherpa prints the name of the hard process for the
      last event at each print out.

:samp:`{16}`
      Sherpa prints the elapsed time and time left in
      seconds only.

The settings are additive such that multiple settings can be employed
at the same time.

.. note::

   When running the code on a cluster or in a grid environment,
   BATCH_MODE should always contain setting 1
   (i.e. ``BATCH_MODE: 1`` or ``3`` or ``5`` etc.).

   The command line option :option:`-b` should therefore not be used
   in this case, see :ref:`Command line`.

.. _NUM_ACCURACY:

NUM_ACCURACY
============

.. index:: NUM_ACCURACY

The targeted numerical accuracy can be specified through
:option:`NUM_ACCURACY`, e.g. for comparing two numbers. This might
have to be reduced if gauge tests fail for numerical reasons.  The
default is ``1E-10``.

.. _SHERPA_CPP_PATH:

SHERPA_CPP_PATH
===============

.. index:: SHERPA_CPP_PATH

The path in which Sherpa will eventually store dynamically created C++
source code.  If not specified otherwise, sets
:option:`SHERPA_LIB_PATH` to ``$SHERPA_CPP_PATH/Process/lib``. This
value can also be set using the command line option :option:`-L`, see
:ref:`Command line`. Both settings can also be set using environment
variables.

.. _SHERPA_LIB_PATH:

SHERPA_LIB_PATH
===============

.. index:: SHERPA_LIB_PATH

The path in which Sherpa looks for dynamically linked libraries from
previously created C++ source code, cf. :ref:`SHERPA_CPP_PATH`.

.. _Event output formats:

Event output formats
====================

.. index:: HepMC3_GenEvent
.. index:: HepMC3_Short
.. index:: LHEF
.. index:: Root
.. index:: FILE_SIZE
.. index:: EVENT_FILE_PATH
.. index:: EVENT_OUTPUT_PRECISION
.. index:: EVENT_OUTPUT
.. index:: EVENT_INPUT

Sherpa provides the possibility to output events in various formats,
e.g. the HepMC format.  The
authors of Sherpa assume that the user is sufficiently acquainted with
these formats when selecting them.

If the events are to be written to file, the parameter
:option:`EVENT_OUTPUT` must be specified together with a file name. An
example would be ``EVENT_OUTPUT: HepMC3_GenEvent[MyFile]``, where
``MyFile`` stands for the desired file base name. More than one output
can also be specified:

.. code-block:: yaml

   EVENT_OUTPUT:
     - HepMC3_GenEvent[MyFile]
     - Root[MyFile]

The following formats are currently available:

:option:`HepMC3_GenEvent`
  Generates output using HepMC3 library. The format of the output is
  set with ``HEPMC3_IO_TYPE: <0|1|2|3|4>`` tag.  The default value is
  0 and corresponds to ASCII GenEvent. Other available options are 1:
  HepEvt 2: ROOT file with every event written as an object of class
  GenEvent. 3: ROOT file with GenEvent objects writen into TTree.
  Otherwise similar to ``HepMC3_GenEvent``.

  The HepMC::GenEvent::m_weights weight vector stores the
  following items: ``[0]`` event weight, ``[1]`` combined matrix
  element and PDF weight (missing only phase space weight information,
  thus directly suitable for evaluating the matrix element value of
  the given configuration), ``[2]`` event weight normalisation (in
  case of unweighted events event weights of ~ +/-1 can be obtained by
  (event weight)/(event weight normalisation)), and ``[3]`` number of
  trials. The total cross section of the simulated event sample can be
  computed as the sum of event weights divided by the sum of the
  number of trials.  This value must agree with the total cross
  section quoted by Sherpa at the end of the event generation run, and
  it can serve as a cross-check on the consistency of the HepMC event
  file.  Note that Sherpa conforms to the Les Houches 2013 suggestion
  (http://phystev.in2p3.fr/wiki/2013:groups:tools:hepmc) of indicating
  interaction types through the GenVertex type-flag.  Multiple event
  weights can also be used, cf.
  :ref:`On-the-fly event weight variations`. The following additional
  customisations can be used.

  ``HEPMC_USE_NAMED_WEIGHTS: <true|false>`` Enable filling weights
  with an associated name. The nominal event weight has the key
  ``Weight``. ``MEWeight``, ``WeightNormalisation`` and ``NTrials``
  provide additional information for each event as described
  above. The default value is ``true``.

  ``HEPMC_EXTENDED_WEIGHTS: <false|true>`` Write additional event
  weight information needed for a posteriori reweighting into the
  WeightContainer, cf. :ref:`A posteriori scale and PDF variations
  using the HepMC GenEvent Output`. Necessitates the use of
  ``HEPMC_USE_NAMED_WEIGHTS``. The default value is ``false``.

  ``HEPMC_TREE_LIKE: <false|true>`` Force the event record to be
  stricly tree-like. Please note that this removes some information
  from the matrix-element-parton-shower interplay which would be
  otherwise stored.
  Requires ``-DHepMC3_DIR=/path/to/hepmc3``
  (or ``-DSHERPA_ENABLE_HEPMC3=ON``, if HepMC3 is installed in a
  standard location). The default value is ``false``.

:option:`HepMC3_Short`

  Generates output in HepMC::IO_GenEvent format, however, only
  incoming beams and outgoing particles are stored. Intermediate and
  decayed particles are not listed. The event weights stored as the
  same as above, and ``HEPMC_USE_NAMED_WEIGHTS`` and
  ``HEPMC_EXTENDED_WEIGHTS`` can be used to customise.

  Requires ``-DHepMC3_DIR=/path/to/hepmc3``
  (or ``-DSHERPA_ENABLE_HEPMC3=ON``, if HepMC3 is installed in a
  standard location).

:option:`LHEF`
  Generates output in Les Houches Event File format. This output
  format is intended for output of **matrix element configurations
  only**. Since the format requires PDF information to be written out
  in the outdated PDFLIB/LHAGLUE enumeration format this is only
  available automatically if LHAPDF is used, the identification
  numbers otherwise have to be given explicitly via
  ``LHEF_PDF_NUMBER`` (``LHEF_PDF_NUMBER_1`` and ``LHEF_PDF_NUMBER_2``
  if both beams carry different structure functions).  This format
  currently outputs matrix element information only, no information
  about the large-Nc colour flow is given as the LHEF output format is
  not suited to communicate enough information for meaningful parton
  showering on top of multiparton final states.

:option:`Root`
  Generates output in ROOT ntuple format **for NLO event generation
  only**.  For details on the ntuple format, see :ref:`A posteriori
  scale and PDF variations using the ROOT NTuple Output <A posteriori
  scale and PDF variations using the ROOT NTuple Output>`. ROOT ntuples can be
  read back into Sherpa and analyzed using the option
  :option:`EVENT_INPUT`. This feature is described in :ref:`NTuple production`.

  Requires ``-DROOT_DIR=/path/to/root``
  (or ``-DSHERPA_ENABLE_ROOT=ON``, if ROOT is installed in a
  standard location).

The output can be further customized using the following options:

:option:`FILE_SIZE`
  Number of events per file (default: unlimited).

:option:`EVENT_FILE_PATH`
  Directory where the files will be stored.

:option:`EVENT_OUTPUT_PRECISION`
  Steers the precision of all numbers written to file (default: 12).

For all output formats except ROOT, events can be written
directly to gzipped files instead of plain text. The option
:option:`-DSHERPA_ENABLE_GZIP=ON` must be given during installation to enable
this feature.

.. _On-the-fly event weight variations:

On-the-fly event weight variations
==================================

Sherpa can compute alternative event weights on-the-fly, resulting in
alternative weights for the generated event.
An important example is the variation of QCD scales and input PDF.
There are also on-the-fly variations for approximate electroweak corrections,
this is discussed in its own section, :ref:`Approximate Electroweak
Corrections`.

Specifying variations
---------------------

There are two ways to specify scale and PDF variations.
Either using the unified ``VARIATIONS`` list,
and/or by using the specialised ``SCALE_VARIATIONS``
and ``PDF_VARIATIONS``, and ``QCUT_VARIATIONS`` lists.
Only the ``VARIATIONS`` list allows to specify
correlated variations (i.e. varying both scales and PDFs at the same time),
but it is more verbose and therefore harder to remember.
Therefore, we suggest to use the more specialised variants
whenever uncorrelated variations are required.

They are evoked using the following syntax:

.. _SCALE_VARIATIONS:
.. _PDF_VARIATIONS:
.. _QCUT_VARIATIONS:

.. index:: SCALE_VARIATIONS
.. index:: PDF_VARIATIONS
.. index:: QCUT_VARIATIONS

.. code-block:: yaml

   SCALE_VARIATIONS:
   - [<muF2-fac-1>, <muR2-fac-1>]
   - [<muF2-fac-2>, <muR2-fac-2>]
   - <mu2-fac-3>

   PDF_VARIATIONS:
   - <PDF-1>
   - <PDF-2>

   QCUT_VARIATIONS:
   - <qcut-fac-1>
   - <qcut-fac-2>

This example specifies a total of seven on-the-fly variations.

Scale factors in ``SCALE_VARIATIONS`` can be given
as a list of two numbers, or as a single number.
When two numbers are given, they are applied to the factorisation and the renomalisation scale, respectively.
If only a single number is given, it is applied to both scales at the same time.
The factors for the renormalisation and factorisation scales
must be given in their quadratic form, i.e. a "4.0" in the settings means that the
(unsquared) scale is to be multiplied by a factor of 2.0.

For the ``PDF_VARIATIONS``, any set present in any of the PDF library
interfaces loaded through ``PDF_LIBRARY`` can be used. If no PDF set is given
it defaults to the nominal one. Specific PDF members can be specified by
appending the PDF set name with ``/<member-id>``.

It can be painful to write every variation explicitly, e.g. for 7-point scale
factor variations or if one want variations for all members of a PDF set.
Therefore an asterisk can be appended to some values, which results in an
*expansion*.  For PDF sets, this means that the variation is repeated for each
member of that set.  For scale factors, ``4.0*`` is expanded to itself, unity,
and its inverse: ``1.0/4.0, 1.0, 4.0``.  A special meaning is reserved for
specifying a single number ``4.0*`` as a ``SCALE_VARIATIONS`` list item,
which expands to a 7-point scale variation:

.. code-block:: yaml

   SCALE_VARIATIONS:
   - 4.0*

is therefore equivalent to

.. code-block:: yaml

   SCALE_VARIATIONS:
   - [0.25, 0.25]
   - [0.25, 1.00]
   - [1.00, 0.25]
   - [1.00, 1.00]
   - [4.00, 1.00]
   - [1.00, 4.00]
   - [4.00, 4.00]

Equivalently, one can even just write ``SCALE_VARIATIONS: 4.0*``,
because a single scalar on the right-hand side will automatically
be interpreted as the first item of a list when the setting
expects a list.

Such expansions may include trivial scale variations and the central
PDF set, resulting
in the specification of a completely trivial variation,
which would just repeat the nominal calculation.
Per default, these trivial variations are automically omitted during the
calculation, since the nominal calculation is anyway included in the Sherpa
output. If required (e.g. for debugging), this filtering
can be explicitly disabled using
``VARIATIONS_INCLUDE_CV: true``.

We now discuss the alternative ``VARIATIONS`` syntax.
The following snippet
specifies two on-the-fly variations,
where scales and PDFs are varied
simultaneously:

.. _VARIATIONS:

.. index:: VARIATIONS

.. code-block:: yaml

   VARIATIONS:
   - ScaleFactors:
       MuR2: <muR2-fac-1>
       MuF2: <muF2-fac-1>
       QCUT: <qcut-fac-1>
     PDF: <PDF-1>
   - ScaleFactors:
       MuR2: <muR2-fac-2>
       MuF2: <muF2-fac-2>
       QCUT: <qcut-fac-2>
     PDF: <PDF-2>
   ...

The key word ``VARIATIONS`` takes a list of variations.  Each variation is
specified by a set of scale factors, and a PDF choice (or AlphaS(MZ) choice,
see below).

Scale factors can be given for the renormalisation, factorisation and for the
merging scale.  The corresponding keys are ``MuR2``, ``MuF2`` and ``QCUT``,
respectively.
The factors for the renormalisation and factorisation scales
must be given in their quadratic form, i.e. a ``MUR2: 4.0`` means that the
(unsquared) renormalisation scale is to be multiplied by a factor of 2.0.
All scale factors can be omitted
(they default to 1.0). Instead of ``MuR2`` and ``MuF2``, one can also use the
keyword ``Mu2``. In this case, the given factor is applied to both the
renormalisation and the factorisation scale.

Instead of using ``PDF: <PDF>`` (which consistently also varies the strong
coupling if the PDF has a different specification of it!), one can also specify
a pure AlphaS variation by giving its value at the Z mass scale: ``AlphaS(MZ):
<alphas(mz)-value>``. This can be useful e.g. for leptonic productions,
and is currently exclusive to the ``VARIATIONS`` syntax.

Also ``VARIATIONS`` can expand values using the star syntax:

.. code-block:: yaml

   VARIATIONS:
     - ScaleFactors:
         Mu2: 4.0*

is therefore equivalent to

.. code-block:: yaml

   VARIATIONS:
     - ScaleFactors:
         MuF2: 0.25
         MuR2: 0.25
     - ScaleFactors:
         MuF2: 1.0
         MuR2: 0.25
     - ScaleFactors:
         MuF2: 0.25
         MuR2: 1.0
     - ScaleFactors:
         MuF2: 1.0
         MuR2: 1.0
     - ScaleFactors:
         MuF2: 4.0
         MuR2: 1.0
     - ScaleFactors:
         MuF2: 1.0
         MuR2: 4.0
     - ScaleFactors:
         MuF2: 4.0
         MuR2: 4.0

As another example, a complete variation using the PDF4LHC convention would
read

.. code-block:: yaml

   VARIATIONS:
     - ScaleFactors:
         Mu2: 4.0*
     - PDF: CT10nlo*
     - PDF: MMHT2014nlo68cl*
     - PDF: NNPDF30_nlo_as_0118*

Please note, this syntax will create :math:`6+52+50+100=208` additional weights
for each event. Even though reweighting is used to reduce the amount of
additional calculation as far as possible, this can still necessitate a
considerable amount of additional CPU hours, in particular when parton-shower
reweighting is enabled (see below).

The rest of this section applies to both the combined ``VARIATIONS``
and the individual ``SCALE_VARIATIONS`` etc. syntaxes.

Variation output
----------------

.. index:: OUTPUT_ME_ONLY_VARIATIONS

The total cross section for all variations along with the nominal cross section
are written to the standard output after the event generation has finalized.
Additionally, some event output (see :ref:`Event output formats`) and analysis methods
(see :ref:`ANALYSIS`) are able to process alternate event weights.
Currently, the only supported event output method is
``HepMC3_GenEvent`` (requires configuration with HepMC version 3 or later).
The supported analysis methods are ``Rivet`` and ``Internal``.

The alternative event weight names follow the MC naming convention, i.e. they
are named ``MUR=<fac>__MUF=<fac>__LHAPDF=<id>``.  When using Sherpa's
interface to Rivet 2, :ref:`Rivet analyses`, separate instances of
Rivet, one for each alternative event weight in addition to the
nominal one, are instantiated leading to one set of histograms each.
They are again named using the ``MUR=<fac>__MUF=<fac>__LHAPDF=<id>``
convention.
For Rivet 3, the internal multi-weight handling capabilities are used instead,
such that there are no alternate histogram files, just one containing
histograms for all variations.
Extending the naming convention, for pure strong coupling variations, an additional
tag ``ASMZ=<val>`` is appended.
If shower scale variations are disabled (either implicitly, because ``SHOWER_GENERATOR: None``,
or explicitly, see below),
you will find ``ME.MUR``/``ME.MUF`` tags instead of the simple ones
to make explicit that the parton-shower scales are not varied with the ME scales.

If parton-shower variations are enabled, ``CSS_REWEIGHT: true``
(the default if parton showering is enabled),
then pure ME-only variations are included along with the full variations in the
HepMC/Rivet output by default. This can be disabled using
``OUTPUT_ME_ONLY_VARIATIONS: false``.
All weight names of ME-only variations
include a "ME" as part of the keys to indicate that
only the ME part of the calculation has been varied, e.g.
``ME:MUR=<fac>__ME:MUF=<fac>__ME:LHAPDF=<id>``.

The user must also be aware that, of course, the cross section of the
event sample, changes when using an alternative event weight as
compared to the nominal one. Any histogramming therefore has to account
for this and recompute the total cross section as the sum of weights
divided by the number of trials, cf. :ref:`Cross section
determination`.
For HepMC 3, Sherpa writes alternate cross sections directly to the
GenCrossSection entry of the event record, such that no manual intervention is
required (as long as the correct cross section variation is picked in
downstream processing steps).

Varying the PDFs of a single beam
---------------------------------

.. _PDF_VARIATION_BEAMS:
.. _PDF_VARIATION_ALPHAS_BEAM:

.. index:: PDF_VARIATION_BEAMS
.. index:: PDF_VARIATION_ALPHAS_BEAM

The ``PDF_VARIATION_BEAMS`` setting can be used to restrict for which beams
a PDF variation is applied. Its default is ``[1, 2]``, i.e. both beams
will undergo a given PDF variation. Use ``1`` or ``2`` to only apply it to a
single beam. This is a global setting for all PDF variations, i.e. it is
currently not possible to do this on the basis of a single PDF
variation.

When using ``PDF_VARIATION_BEAMS``,
there is an ambiguity which beam's PDF
should be used to evaluate the strong coupling.
For that, the setting ``PDF_VARIATION_ALPHAS_BEAM``
can be used.  Its default is ``0``, which means that the first
available beam's PDF is used. Use ``1`` or ``2`` to select a specific beam's
PDF instead.

Having different PDFs for each beam will be reflected in the
:ref:`Variation output`.
Consider the following example:
``MUR=1__MUF=1__LHAPDF.BEAM1=93300__LHAPDF.BEAM2=93301``,
where the beams' LHAPDF IDs are specified individually.

Variations for different event generation modes
-----------------------------------------------

The on-the-fly reweighting works for all event generation modes
(weighted or (partially) unweighted) and all calculation types (LO,
LOPS, NLO, NLOPS, NNLO, NNLOPS, MEPS\@LO, MEPS\@NLO and MENLOPS).

NLO calculations
````````````````

.. index:: NLO_MUR_COEFFICIENT_FROM_VIRTUAL

For NLO calculations, note that some loop providers (e.g. Recola)
do not provide the pole coefficients, while others do (e.g. OpenLoops).
For the former, Sherpa will automatically exclude the IR pole
coefficients from the scale variation.
One can also manually exclude them using
``NLO_MUR_COEFFICIENT_FROM_VIRTUAL: false``.
If they are excluded,
then IR pole cancellation is assumed and, thus,
only the UV renormalisation term pole coefficient is considered in the scale variation.

Parton shower emissions
```````````````````````

.. index:: CSS_REWEIGHT
.. index:: CSS_REWEIGHT_SCALE_CUTOFF
.. index:: CSS_MAX_REWEIGHT_FACTOR

By default, the reweighting of parton shower emissions is included in the variations.
It can be disabled explicitly,
using :option:`CSS_REWEIGHT: false`.  This should work out of the box for all
types of variations. However, parton-shower reweighting (even though formally
exact), tends to be numerically less stable than the reweighting of the hard
process. If numerical issues are encountered, one can try to
increase :option:`CSS_REWEIGHT_SCALE_CUTOFF` (default: 5, measured in GeV).
This disables shower variations for emissions at scales below the value.
An additional safeguard against rare spuriously large shower variation
weights is implemented as :option:`CSS_MAX_REWEIGHT_FACTOR` (default: 1e3).
Any variation weights accumulated during an event and larger than this factor
will be ignored and reset to 1.

.. _MPI parallelization:

MPI parallelization
===================

MPI parallelization in Sherpa can be enabled using the configuration
option :option:`-DSHERPA_ENABLE_MPI=ON`. Sherpa supports `OpenMPI
<http://www.open-mpi.org/>`_ and `MPICH2
<http://www.mcs.anl.gov/research/projects/mpich2/>`_ . For detailed
instructions on how to run a parallel program, please refer to the
documentation of your local cluster resources or the many excellent
introductions on the internet. MPI parallelization is mainly intended
to speed up the integration process, as event generation can be
parallelized trivially by starting multiple instances of Sherpa with
different random seed, cf.  :ref:`RANDOM_SEED`. However, both the
internal analysis module and the Root NTuple writeout can be used with
MPI. Note that these require substantial data transfer.

Please note that the process information contained in the ``Process``
directory for both Amegic and Comix needs to be generated without MPI
parallelization first. Therefore, first run

.. code-block:: shell-session

   $ Sherpa INIT_ONLY=1 <Sherpa.yaml>

and, in case of using Amegic, compile the libraries. Then start your
parallized integration, e.g.

.. code-block:: shell-session

   $ mpirun -n <n> Sherpa -e 0 <Sherpa.yaml>

After the integration has finished, you can submit individual jobs to generate
event samples (with a different random seed for each job).  Upon completion,
the results can be merged.
