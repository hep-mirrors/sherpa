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
compiled with option :option:`--enable-gzip`.  When set to 3, Sherpa
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
  | To use this option, the package must be configured with option :option:`--enable-analysis`.
  | An output directory can be specified using :ref:`ANALYSIS_OUTPUT`.

:option:`Rivet`
  | The Rivet package, see `Rivet Website <http://projects.hepforge.org/rivet/>`_.
  | To enable it, Rivet and HepMC have to be installed and Sherpa must be configured
  | as described in :ref:`Rivet analyses`.

:option:`HZTool`
  | The HZTool package, see `HZTool Website <http://hztool.hepforge.org/>`_.
  | To enable it, HZTool and CERNLIB have to be installed and Sherpa must be configured
  | as described in :ref:`HZTool analyses`.


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

The settings are additive such that multiple settings can be employed
at the same time.

.. note::

   When running the code on a cluster or in a grid environment,
   BATCH_MODE should always contain setting 1
   (i.e. ``BATCH_MODE=[1|3|5|7]``).

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

.. index:: HepMC_GenEvent
.. index:: HepMC_Short
.. index:: HEPEVT
.. index:: LHEF
.. index:: Root
.. index:: Delphes
.. index:: FILE_SIZE
.. index:: EVENT_FILE_PATH
.. index:: EVENT_OUTPUT_PRECISION
.. index:: EVENT_OUTPUT
.. index:: EVENT_INPUT

Sherpa provides the possibility to output events in various formats,
e.g. the HepEVT common block structure or the HepMC format.  The
authors of Sherpa assume that the user is sufficiently acquainted with
these formats when selecting them.

If the events are to be written to file, the parameter
:option:`EVENT_OUTPUT` must be specified together with a file name. An
example would be ``EVENT_OUTPUT: HepMC_GenEvent[MyFile]``, where
``MyFile`` stands for the desired file base name. More than one output
can also be specified:

.. code-block:: yaml

   EVENT_OUTPUT:
     - HepMC_GenEvent[MyFile]
     - Root[MyFile]

The following formats are currently available:

:option:`HepMC_GenEvent`
  Generates output in HepMC::IO_GenEvent
  format. The HepMC::GenEvent::m_weights weight vector stores the
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
  weights can also be enabled with HepMC versions >=2.06, cf.
  :ref:`Scale and PDF variations`. The following additional
  customisations can be used

  ``HEPMC_USE_NAMED_WEIGHTS: <false|true>`` Enable filling weights
  with an associated name. The nominal event weight has the key
  ``Weight``. ``MEWeight``, ``WeightNormalisation`` and ``NTrials``
  provide additional information for each event as described
  above. Needs HepMC version >=2.06.

  ``HEPMC_EXTENDED_WEIGHTS: <false|true>`` Write additional event
  weight information needed for a posteriori reweighting into the
  WeightContainer, cf. :ref:`A posteriori scale and PDF variations
  using the HepMC GenEvent Output`. Necessitates the use of
  ``HEPMC_USE_NAMED_WEIGHTS``.

  ``HEPMC_TREE_LIKE: <false|true>`` Force the event record to be
  stricly tree-like. Please note that this removes some information
  from the matrix-element-parton-shower interplay which would be
  otherwise stored.

:option:`HepMC_Short`

  Generates output in HepMC::IO_GenEvent format, however, only
  incoming beams and outgoing particles are stored. Intermediate and
  decayed particles are not listed. The event weights stored as the
  same as above, and ``HEPMC_USE_NAMED_WEIGHTS`` and
  ``HEPMC_EXTENDED_WEIGHTS`` can be used to customise.

:option:`HepMC3_GenEvent`
  Generates output using HepMC3 library. The format of the output is
  set with ``HEPMC3_IO_TYPE: <0|1|2|3|4>`` tag.  The default value is
  0 and corresponds to ASCII GenEvent. Other available options are 1:
  HepEvt 2: ROOT file with every event written as an object of class
  GenEvent. 3: ROOT file with GenEvent objects writen into TTree.
  Otherwise similar to ``HepMC_GenEvent``.

:option:`Delphes_GenEvent`
  Generates output in `Root <http://root.cern.ch>`_ format, which can
  be passed to `Delphes <http://cp3.irmp.ucl.ac.be/projects/delphes>`_
  for analyses.  Input events are taken from the HepMC
  interface. Storage space can be reduced by up to 50% compared to
  gzip compressed HepMC. This output format is available only if
  Sherpa was configured and installed with options
  :option:`--enable-root` and
  :option:`--enable-delphes=/path/to/delphes`.

:option:`Delphes_Short`
  Generates output in `Root <http://root.cern.ch>`_ format, which can
  be passed to `Delphes <http://cp3.irmp.ucl.ac.be/projects/delphes>`_
  for analyses.  Only incoming beams and outgoing particles are
  stored.

:option:`PGS`
  Generates output in `StdHEP <http://cepa.fnal.gov/psm/stdhep>`_
  format, which can be passed to `PGS
  <http://www.physics.ucdavis.edu/~conway/research/software/pgs/pgs4-general.htm>`_
  for analyses. This output format is available only if Sherpa was
  configured and installed with options
  :option:`--enable-hepevtsize=4000` and
  :option:`--enable-pgs=/path/to/pgs`.  Please refer to the PGS
  documentation for how to pass StdHEP event files on to PGS.  If you
  are using the LHC olympics executeable, you may run
  ``./olympics --stdhep events.lhe <other options>``.

:option:`PGS_Weighted`
  Generates output in `StdHEP <http://cepa.fnal.gov/psm/stdhep>`_
  format, which can be passed to `PGS
  <http://www.physics.ucdavis.edu/~conway/research/software/pgs/pgs4-general.htm>`_
  for analyses. Event weights in the HEPEV4 common block are stored in
  the event file.

:option:`HEPEVT`
  Generates output in HepEvt format.

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
  scale and PDF variations using the ROOT NTuple Output>`.  This
  output option is available only if Sherpa was linked to ROOT during
  installation by using the configure option
  ``--enable-root=/path/to/root``.  ROOT ntuples can be read back into
  Sherpa and analyzed using the option :option:`EVENT_INPUT`. This
  feature is described in :ref:`NTuple production`.

The output can be further customized using the following options:

:option:`FILE_SIZE`
  Number of events per file (default: unlimited).

:option:`EVENT_FILE_PATH`
  Directory where the files will be stored.

:option:`EVENT_OUTPUT_PRECISION`
  Steers the precision of all numbers written to file (default: 12).

For all output formats except ROOT and Delphes, events can be written
directly to gzipped files instead of plain text. The option
:option:`--enable-gzip` must be given during installation to enable
this feature.

.. _Scale and PDF variations:

Scale and PDF variations
========================

Sherpa can compute alternative event weights for different scale and
PDF choices on-the-fly, resulting in alternative weights for the generated
event. The can be evoked with the following syntax

.. code-block:: yaml

   VARIATIONS:
   - "<muR-fac-1>,<muF-fac-1>,<PDF-1>"
   - "<muR-fac-2>,<muF-fac-2>,<PDF-2>"
   ...

The key word ``VARIATIONS`` takes a list of
variation factors for the nominal renormalisation and factorisation scale
and an associated PDF set. If the scale factors are omitted, they default to 1.
Any set present in any of the PDF library interfaces
loaded through ``PDF_LIBRARY`` can be used. If no PDF set is given it
defaults to the nominal one. Specific PDF members can be
specified by appending the PDF set name with ``/<member-id>``.  Enclosing
the pdf set with square brackets will expand to variations over all members of
that set. This only works with LHAPDF6 sets or the internal default sets.
Please note that scales are, as always in Sherpa, given in their quadratic
form.  Thus, a variation of factor 4 of the squared scale [GeV^2] means a
variation of factor 2 on the scale itself [GeV].
Scales also support square bracket expansion, e.g. ``[4]`` expands to
``1/4, 1, 4``. Enclosing both scale tags with square brackets expands to the
7-point scale variation:

.. code-block:: yaml

   VARIATIONS:
     - "[4,4]"

   # is equivalent to
   VARIATIONS:
     - 0.25,0.25
     - 1,0.25
     - 0.25,1
     - 1,1
     - 4,1
     - 1,4
     - 4,4

Thus, a complete variation using the PDF4LHC convention would read

.. code-block:: yaml

   VARIATIONS:
     - "[4,4]"
     - "[CT10nlo]"
     - "[MMHT2014nlo68cl]"
     - "[NNPDF30_nlo_as_0118]"

Please note, this syntax will create :math:`7+53+51+101=212`
additional weights for each event.

Note that the square bracket expansion includes trivial scale
variations and the central PDF set. This can be disabled with
``VARIATIONS_INCLUDE_CV: false``.

The additional event weights can then be written into the event
output.  However, this is currently only supported for
``HepMC_GenEvent`` and ``HepMC_Short`` with versions >=2.06 and
``HEPMC_USE_NAMED_WEIGHTS: true``.  The alternative event weights
follow the Les Houches naming convention for such variations, ie. they
are named ``MUR<fac>_MUF<fac>_PDF<id>``.  When using Sherpa's
interface to Rivet, :ref:`Rivet analyses`, separate instances of
Rivet, one for each alternative event weight in addition to the
nominal one, are instantiated leading to one set of histograms each.
They are again named using the ``MUR<fac>_MUF<fac>_PDF<id>``
convention.

The user must also be aware that, of course, the cross section of the
event sample, changes when using an alternative event weight as
compared to the nominal one. Any histograming therefore has to account
for this and recompute the total cross section as the sum of weights
devided by the number of trials, cf. :ref:`Cross section
determination`.

The on-the-fly reweighting works for all event generation modes
(weighted or (partially) unweighted) and all calculation types (LO,
LOPS, NLO, NLOPS, MEPS\@LO, MEPS\@NLO and MENLOPS).
However, the reweighting of parton shower emissions has to be enabled explicitly,
using :option:`CSS_REWEIGHT: true`.  This should work out of the box for both scale
and PDF variations.  If numerical issues are encountered, one can try to
increase :option:`CSS_REWEIGHT_SCALE_CUTOFF` (default: 5, measured in GeV).
This disables shower variations for emissions at scales below the value.
An additional safeguard against rare spuriously large shower variation
weights is implemented as @code{CSS_MAX_REWEIGHT_FACTOR} (default: 1e3).
Any variation weights accumulated during an event and larger than this factor
will be ignored and reset to 1.

To include the ME-only variations along with the full variations in the
HepMC/Rivet output, you can use ``HEPMC_INCLUDE_ME_ONLY_VARIATIONS:
true`` and ``RIVET: @{ INCLUDE_HEPMC_ME_ONLY_VARIATIONS: true @``},
respectively.

.. _MPI parallelization:

MPI parallelization
===================

MPI parallelization in Sherpa can be enabled using the configuration
option :option:`--enable-mpi`. Sherpa supports `OpenMPI
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

   $ Sherpa -f <run-card> INIT_ONLY=1

and, in case of using Amegic, compile the libraries. Then start your
parallized integration, e.g.

.. code-block:: shell-session

   $ mpirun -n <n> Sherpa -f <run-card> -e 0

After the integration has finished, you can submit individual jobs to generate
event samples (with a different random seed for each job).  Upon completion,
the results can be merged.
