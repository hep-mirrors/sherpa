.. _Getting started:

###############
Getting started
###############


.. contents::
   :local:

.. _Installation:

************
Installation
************


Sherpa is distributed as a tarred and gzipped file named
:samp:`sherpa-{<VERSION>}.tar.gz`, and can be unpacked in the
current working directory with

.. code-block:: shell-session

   $ tar -zxf sherpa-<VERSION>.tar.gz

Alternatively, it can also be accessed via Git through the location
specified on the download page.

To guarantee successful installation, the following tools should be
available on the system:

  * C++ compiler
  * cmake
  * make or ninja


Recommended:
  * Fortran compiler
  * LHAPDF  (including devel packages). If not available, use the `-DSHERPA_ENABLE_INSTALL_LHAPDF=ON` cmake option to install LHAPDF on-the-fly during the Sherpa installation (internet connection required).
  * libzip  (including devel packages). If not available, use the `-DSHERPA_ENABLE_INSTALL_LIBZIP=ON` cmake option to install libzip on-the-fly during the Sherpa installation (internet connection required).

Compilation and installation proceed through the following commands if
you use the distribution tarball:

.. code-block:: shell-session

   $ cd sherpa-<VERSION>/
   $ cmake -S . -B <builddir> [+ optional configuration options described below]
   $ cmake --build <builddir> [other build options, e.g. -j 8]
   $ cmake --install <builddir>

where `<builddir>` has to be replaced with the (temporary) directory in which intermediate files are stored during the build process.
You can simply use the current working directory, i.e. `cmake -S . -B .` to compile in-source if you want to keep everything Sherpa-related in one directory.

Note that re-running ``cmake`` with different configuration options is `not the same <https://gitlab.kitware.com/cmake/cmake/-/issues/19622>`_ as running it in a fresh working directory. Use ``ccmake .`` instead to check/change the current configuration. To start afresh, e.g. to pick up a different version of a dependency, you can use the ``cmake --fresh [...]`` option in recent versions of cmake, or just delete the cache (``rm -rf CMakeCache.txt CMakeFiles``).

If not specified differently, the directory structure after
installation is organized as follows


``$(prefix)/bin``
  Sherpa executable and scripts

``$(prefix)/include``
  headers for process library compilation

``$(prefix)/lib``
  basic libraries

``$(prefix)/share``
  PDFs, Decaydata, fallback run cards


The installation directory ``$(prefix)`` can be specified by using the
``-DCMAKE_INSTALL_PREFIX=/path/to/installation/target`` directive and
defaults to the current working directory (`.`).

If Sherpa has to be moved to a different directory after the
installation, one has to set the following environment variables for
each run:

  * ``SHERPA_INCLUDE_PATH=$newprefix/include/SHERPA-MC``
  * ``SHERPA_SHARE_PATH=$newprefix/share/SHERPA-MC``
  * ``SHERPA_LIBRARY_PATH=$newprefix/lib/SHERPA-MC``
  * ``LD_LIBRARY_PATH=$SHERPA_LIBRARY_PATH:$LD_LIBRARY_PATH``


Sherpa can be interfaced with various external packages,
e.g. `HepMC <http://hepmc.web.cern.ch/hepmc/>`_ for event output,
`LHAPDF <https://lhapdf.hepforge.org/>`_ for PDF sets,
or `Rivet <https://rivet.hepforge.org/>`_ for analysis.
For this to work, the user has to add the corresponding options to the
cmake configuration, e.g. for Rivet:

.. code-block:: shell-session

   $  cmake [...] -DSHERPA_ENABLE_RIVET=ON

If your Rivet installation is not in a standard directory, you instead have to
point cmake to the path where Rivet is installed:

.. code-block:: shell-session

   $  cmake [...] -DRIVET_DIR=/my/rivet/install/dir

The provided path has to point to the top level installation
directory of the external package, i.e. the one containing the
``lib/``, ``share/``, ... subdirectories.

Other external packages are activated using equivalent configuration options,
i.e. either using `-DSHERPA_ENABLE_<PACKAGENAME>=ON`
or using `-D<PackageName>_DIR=/my/package/install/dir`
(or both, but enabling a package is already implied
if its directory is given).
Note that the package name in the `SHERPA_ENABLE_<PACKAGENAME>`
is always capitalised,
while the capitalisation can differ
in `<PackageName>_DIR`,
as defined by the third-party package.
For a complete list of possible configuration options
(and their correct capitalisation),
run ``cmake -LA``.

The Sherpa package has successfully been compiled, installed and tested
on various Linux distributions
(Arch Linux, SUSE Linux, RHEL, Scientific Linux, Debian, Ubuntu) and on macOS,
using the GNU Compiler Collection (GCC), Clang and Intel OneAPI 2022.

If you have multiple compilers installed on your system, you can use
shell environment variables to specify which of these are to be
used. A list of the available variables is printed with

.. code-block:: shell-session

   $ -DCMAKE_CXX_COMPILER=myc++compiler

in the Sherpa top level directory and looking at the last
lines. Depending on the shell you are using, you can set these
variables e.g. with export (bash) or setenv (csh). For example,
to use a specific versioned GCC installation, you could run
before calling ``cmake``:

.. code-block:: bash

   export CXX=g++-13
   export CC=gcc-13
   export CPP=cpp-13


Installation on Cray XE6 / XK7
==============================

Sherpa has been installed successfully on Cray XE6 and Cray XK7.  The
following configure command should be used:

.. code-block:: shell-session

   $ cmake -DSHERPA_ENABLE_MPI=ON <your options>

Sherpa can then be run with

.. code-block:: shell-session

   $ aprun -n <nofcores> <prefix>/bin/Sherpa -lrun.log

The modularity of the code requires setting the environment variable
:option:`CRAY_ROOTFS`, cf. the Cray system documentation.

Installation on IBM BlueGene/Q
==============================

Sherpa has been installed successfully on an IBM BlueGene/Q system.
The following cmake command should be used

.. code-block:: shell-session

   $ cmake <your options> -DSHERPA_ENABLE_MPI=ON -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_Fortran_COMPILER=mpif90

Sherpa can then be run with

.. code-block:: shell-session

   $ qsub -A <account> -n <nofcores> -t 60 --mode c16 <prefix>/bin/Sherpa -lrun.log

macOS Installation
==================

The installation on macOS works analogously to an installation
on a GNU/Linux system.
You might need to ensure first that the Xcode Command Line Tools are installed.
Other missing command line tools
can be installed through a package manager like `Homebrew <https://brew.sh/>`_
or `MacPorts <https://www.macports.org>`_.

.. _Running Sherpa:

**************
Running Sherpa
**************

The ``Sherpa`` executable resides in the directory ``<prefix>/bin/``
where ``<prefix>`` denotes the path to the Sherpa installation
directory. The way a particular simulation will be accomplished is
defined by several parameters, which can all be listed in a common
file, or data card (Parameters can be alternatively specified on the
command line; more details are given in :ref:`Input structure`).  This
steering file is called ``Sherpa.yaml`` and some example setups
(i.e. ``Sherpa.yaml`` files) are distributed with Sherpa.
They can be found in the directory ``<prefix>/share/SHERPA-MC/Examples/``,
and descriptions of some of their key features can be found
in the section :ref:`Examples`.

.. note:: It is not in general possible to reuse steering files from
   previous Sherpa versions. Often there are small changes in the
   parameter syntax of the files from one version to the next.  These
   changes are documented in our manuals. In addition, update any
   custom Decaydata directories you may have used (and reapply any
   changes which you might have applied to the old ones), see
   :ref:`Hadron decays`.

The very first step in running Sherpa is therefore to adjust all
parameters to the needs of the desired simulation. The details for
doing this properly are given in :ref:`Parameters`. In this section,
the focus is on the main issues for a successful operation of
Sherpa. This is illustrated by discussing and referring to the
parameter settings that come in the example steering file
``./Examples/V_plus_Jets/LHC_ZJets/Sherpa.LO.yaml``,
cf. :ref:`LHC_ZJets`.  This is a simple configuration created to show
the basics of how to operate Sherpa. **It should be stressed
that this steering file relies on many of Sherpa's default settings,
and, as such, you should understand those settings before using it to
look at physics.** For more information on the settings and parameters
in Sherpa, see :ref:`Parameters`, and for more examples see the
:ref:`Examples` section.

.. _Process selection and initialization:

Process selection and initialization
====================================

Central to any Monte Carlo simulation is the choice of the hard
processes that initiate the events. These hard processes are described
by matrix elements. In Sherpa, the selection of processes happens in
the ``PROCESSES`` part of the steering file.  Only a few 2->2
reactions have been hard-coded. They are available in the EXTRA_XS
module.  The more usual way to compute matrix elements is to employ
one of Sherpa's automated tree-level generators, AMEGIC++ and Comix,
see :ref:`Basic structure`.  If no matrix-element generator is
selected, using the :ref:`ME_GENERATORS` tag, then Sherpa will use
whichever generator is capable of calculating the process, checking
Comix first, then AMEGIC++ and then EXTRA_XS. Therefore, for some
processes, several of the options are used. In this example, however,
all processes will be calculated by Comix.

To begin with the example, the Sherpa run has to be started by
changing into the
``<prefix>/share/SHERPA-MC/Examples/V_plus_Jets/LHC_ZJets/`` directory
and executing

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa

You may also run from an arbitrary directory, employing
``<prefix>/bin/Sherpa --path=<prefix>/share/SHERPA-MC/Examples/V_plus_Jets/LHC_ZJets``.
In the example, an absolute path is passed to the optional argument
--path.  It may also be specified relative to the current working
directory. If it is not specified at all, the current working
directory is understood.

For good book-keeping, it is highly recommended to reserve different
subdirectories for different simulations as is demonstrated with
the example setups.

If AMEGIC++ is used, Sherpa requires an initialization run, where C++
source code is written to disk. This code must be compiled into
dynamic libraries by the user by running the makelibs script in the
working directory.  After this step Sherpa is run again for the
actual cross section integrations and event generation.  For more
information on and examples of how to run Sherpa using AMEGIC++, see
:ref:`Running Sherpa with AMEGIC++`.

If the internal hard-coded matrix elements or Comix are used, and
AMEGIC++ is not, an initialization run is not needed, and Sherpa will
calculate the cross sections and generate events during the first run
already.

As the cross sections get integrated, the integration over phase space
is optimized to arrive at an efficient event generation.  Subsequently
events are generated if a number of events is passed to the optional
argument :option:`--events` or set in the :file:`Sherpa.yaml` file with the
:ref:`param_EVENTS` parameters.

The generated events are not stored into a file by default; for
details on how to store the events see :ref:`Event output
formats`. Note that the computational effort to go through this
procedure of generating, compiling and integrating the matrix elements
of the hard processes depends on the complexity of the parton-level
final states. For low multiplicities (2->2,3,4), of course, it can be
followed instantly.

.. _Results directory:

Usually more than one generation run is wanted. As long as the
parameters that affect the matrix-element integration are not changed,
it is advantageous to store the cross sections obtained during the
generation run for later use. This saves CPU time especially for large
final-state multiplicities of the matrix elements. Per default, Sherpa
stores these integration results in a directory called :file:`Results/`.
The name of the output directory can be customised via
`Results directory`_


.. code-block:: shell-session

   <prefix>/bin/Sherpa -r <result>/

or with ``RESULT_DIRECTORY: <result>/`` in the steering file, see
:ref:`RESULT_DIRECTORY`. The storage of the integration results can be
prevented by either using

.. code-block:: shell-session

   <prefix>/bin/Sherpa -g

or by specifying ``GENERATE_RESULT_DIRECTORY: false`` in the steering
file.

If physics parameters change, the cross sections have to be
recomputed.  The new results should either be stored in a new
directory or the ``<result>`` directory may be re-used once it has
been emptied.  Parameters which require a recomputation are any
parameters affecting the :ref:`Models`, :ref:`Matrix Elements` or
:ref:`Selectors`.  Standard examples are changing the magnitude of
couplings, renormalisation or factorisation scales, changing the PDF
or centre-of-mass energy, or, applying different cuts at the parton
level. If unsure whether a recomputation is required, a simple test is
to temporarily use a different value for the ``RESULT_DIRECTORY``
option and check whether the new integration numbers (statistically)
comply with the stored ones.

A warning on the validity of the process libraries is in order here:
it is absolutely mandatory to generate new library files, whenever the
physics model is altered, i.e. particles are added or removed and
hence new or existing diagrams may or may not anymore contribute to
the same final states.  Also, when particle masses are switched on or
off, new library files must be generated (however, masses may be
changed between non-zero values keeping the same process
libraries). The best recipe is to create a new and separate setup
directory in such cases. Otherwise the ``Process`` and ``Results``
directories have to be erased:

.. code-block:: shell-session

   $ rm -rf Process/ Results/

In either case one has to start over with the whole initialization
procedure to prepare for the generation of events.



The example set-up: Z+Jets at the LHC
=====================================

The setup file (:file:`Sherpa.yaml`) provided in
``./Examples/V_plus_Jets/LHC_ZJets/`` can be considered as a standard
example to illustrate the generation of fully hadronised events in
Sherpa, cf. :ref:`LHC_ZJets`. Such events will include effects from
parton showering, hadronisation into primary hadrons and their
subsequent decays into stable hadrons. Moreover, the example chosen
here nicely demonstrates how Sherpa is used in the context of merging
matrix elements and parton showers :cite:`Hoeche2009rj`. In addition
to the aforementioned corrections, this simulation of inclusive
Drell-Yan production (electron-positron channel) will then include
higher-order jet corrections at the tree level. As a result the
transverse-momentum distribution of the Drell-Yan pair and the
individual jet multiplicities as measured by the ATLAS and CMS
collaborations at the LHC can be well described.

Before event generation, the initialization procedure as described in
:ref:`Process selection and initialization` has to be completed. The
matrix-element processes included in the setup are the following: ::

  proton proton -> parton parton -> electron positron + up to five partons


In the ``PROCESSES`` list of the steering file this translates into

.. code-block:: yaml

   PROCESSES:
   - 93 93 -> 11 -11 93{5}:
       Order: {QCD: 0, EW: 2}
       CKKW: 20
     [...]

Fixing the order of electroweak
couplings to :option:`2`, matrix elements of all partonic subprocesses
for Drell-Yan production without any and with up to five extra QCD
parton emissions will be generated.  Proton--proton collisions are
considered at beam energies of 6.5 TeV.
Model parameters and couplings can all be defined in
the :file:`Sherpa.yaml` file as you will see in the rest of this manual.

The QCD radiation matrix elements have to be regularised to obtain
meaningful cross sections. This is achieved by specifying ``CKKW: 20``
when defining the process in :file:`Sherpa.yaml`. Simultaneously, this
tag initiates the ME-PS merging procedure.  To eventually obtain fully
hadronised events, the ``FRAGMENTATION`` setting has been left on it's
default value :option:`Ahadic` (and therefore been omitted from the
steering file), which will run Sherpa's cluster hadronisation, and the
``DECAYMODEL`` setting has it's default value :option:`Hadrons`, which
will run Sherpa's hadron decays. Additionally corrections owing to
photon emissions are taken into account.

For a first example run with this setup, we suggest to simplify the run card
significantly and only later, for physics studies, going back to the
full-featured run card. So replace the full process listing with
a short and simple

.. code-block:: yaml

   PROCESSES:
   - 93 93 -> 11 -11 93{1}:
       Order: {QCD: 0, EW: 2}
       CKKW: 20

for now.
In addition, remove the :option:`ASSOCIATED_CONTRIBUTIONS_VARIATIONS` listing,
since these can no longer be calculated after modifying :option:`PROCESSES`.
You might also need to remove (or install) PDFs that are not available
on your system from the :option:`PDF_VARIATIONS` listing.
Then you can go ahead and start Sherpa for the first time by running the

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa

command as described in :ref:`Running Sherpa`. Sherpa displays some
output as it runs. At the start of the run, Sherpa initializes the
relevant model, and displays a table of particles, with their
:ref:`PDG codes` and some properties. It also displays the
:ref:`Particle containers`, and their contents. The other relevant
parts of Sherpa are initialized, including the matrix element
generator(s). The Sherpa output will look like:

.. code-block:: console

   Welcome to Sherpa, <user name> on <host name>.
   Initialization of framework underway ...
   [...]
   Seed: 1234
   [...]
   Initializing beam spectra ...
     Type: Collider Setup
     Beam 1: P+ (enabled = 0, momentum = (6500,0,0,6500))
     Beam 2: P+ (enabled = 0, momentum = (6500,0,0,-6500))
   Initializing PDFs ...
     Hard scattering:    PDF4LHC21_40_pdfas + PDF4LHC21_40_pdfas
     MPI:                PDF4LHC21_40_pdfas + PDF4LHC21_40_pdfas
   [...]
   Fixed electroweak parameters
     Input scheme: alpha(mZ)-mZ-sin(theta_W) scheme, input: 1/\alphaQED(m_Z), sin^2(theta_W), m_Z, m_h, widths
     Ren. scheme:  alphamZsW
     Parameters:   sin^2(\theta_W) = 0.23113
                   vev              = 246.16 - 3.36725 i
   Set \alpha according to EW scheme
     1/\alpha(0)   = 128.802
     1/\alpha(def) = 128.802
   Particle data:
   [...]
   Initializing showers ...
   Initializing matrix elements for the hard processes ...
   Building processes (3 ME generators, 1 process blocks) ...
   Setting up processes ........ done (59 MB, 0s/0s)
   Performing tests ........ done (60 MB, 0s/0s)
   [...]
   Initializing hadron particle information ...
   Initialized fragmentation
   Initialized hadron decays (model = HADRONS++)
   Initialized soft photons
   [...]

Then Sherpa will start to integrate the cross sections. The output
will look like:

.. code-block:: console

   Calculating xs for '2_2__j__j__e-__e+' (Comix) ...
   Integration parameters: n_{min} = 5000, N_{opt} = 10, N_{max} = 1, exponent = 0.5
   Starting the calculation at 09:43:52. Lean back and enjoy ... .
   1630.48 pb +- ( 54.9451 pb = 3.36988 % ) 5000 ( 5029 -> 99.4 % )
   full optimization:  ( 0s elapsed / 5s left ) [09:43:52]
   1692.37 pb +- ( 46.7001 pb = 2.75945 % ) 12071 ( 12133 -> 99.5 % )
   full optimization:  ( 0s elapsed / 5s left ) [09:43:52]
   [...]

The first line here displays the process which is being calculated. In
this example, the integration is for the :math:`2 \to 2` process, parton, parton
:math:`\to` electron, positron. The matrix element generator used is displayed
after the process.  As the integration progresses, summary lines are
displayed, like the ones shown above. The current estimate of the cross
section is displayed, along with its statistical error estimate. The
number of phase space points calculated is displayed after this,
and the efficiency is displayed after that. On the line below, the time
elapsed is shown, and an estimate of the total time till the optimisation
is complete.  In square brackets is an output of the system clock.

When the integration is complete, the output will look like:

.. code-block:: console

   [...]
   1677.75 pb +- ( 1.74538 pb = 0.104031 % ) 374118 ( 374198 -> 100 % )
   full optimization:  ( 4s elapsed / 1s left ) [09:43:56]
   1677.01 pb +- ( 1.36991 pb = 0.0816873 % ) 534076 ( 534157 -> 99.9 % )
   integration time:   ( 5s elapsed / 0s left ) [09:43:58]
   2_2__j__j__e-__e+ : 1677.01 pb +- ( 1.36991 pb = 0.0816873 % )  exp. eff: 20.6675 %
     reduce max for 2_2__j__j__e-__e+ to 1 ( eps = 0.001 -> exp. eff 0.206675 ) 

with the final cross section result and its statistical error displayed.

Sherpa will then move on to integrate the other processes specified in the
run card.

When the integration is complete, the event generation will start.  As
the events are being generated, Sherpa will display a summary line
stating how many events have been generated, and an estimate of how
long it will take.  When the event generation is complete, Sherpa's
output looks like:

.. code-block:: console

   [...]
     Event 100 ( 1 s total ) = 1.12208e+07 evts/day                    
   Summarizing the run may take some time ...
   +----------------------------------------------------------------------------+
   | Nominal or variation name     XS [pb]      RelDev  AbsErr [pb]      RelErr |
   +----------------------------------------------------------------------------+
   | Nominal                       1739.79         0 %      171.304      9.84 % |
   | ME & PS: MUR=0.5 MUF=0.5      1635.61     -5.98 %      187.894     11.48 % |
   | ME & PS: MUR=2 MUF=2          2261.57     29.99 %      387.031     17.11 % |
   | [<results for other variations>]                                           |
   +----------------------------------------------------------------------------+

A summary of the number of events generated is displayed, with the
total cross section for the process and possible systematic variations
:cite:`Bothmann2016nao` and :ref:`Input structure`.

The generated events are not stored into a file by default; for
details on how to store the events see :ref:`Event output formats`.


.. _Parton-level event generation with Sherpa:

Parton-level event generation with Sherpa
=========================================

Sherpa has its own tree-level matrix-element generators called
AMEGIC++ and Comix.  Furthermore, with the module PHASIC++,
sophisticated and robust tools for phase-space integration are
provided. Therefore Sherpa obviously can be used as a cross-section
integrator. Because of the way Monte Carlo integration is
accomplished, this immediately allows for parton-level event
generation. Taking the ``LHC_ZJets`` setup, users have to modify just
a few settings in ``Sherpa.yaml`` and would arrive at a parton-level
generation for the process gluon down-quark to electron positron and
down-quark, to name an example. When, for instance, the options
"``EVENTS: 0``" and "``OUTPUT: 2``" are added to the steering file, a
pure cross-section integration for that process would be obtained with
the results plus integration errors written to the screen.

For the example, the process definition in ``PROCESSES`` simplifies to

.. code-block:: yaml

   - 21 1 -> 11 -11 1:
       Order: {QCD: 1, EW: 2}

with all other settings in the process block removed.  And under the
assumption to start afresh, the initialization procedure has to be
followed as before.  Picking the same collider environment as in the
previous example there are only two more changes before the
:file:`Sherpa.yaml` file is ready for the calculation of the hadronic
cross section of the process g d to e- e+ d at LHC and subsequent
parton-level event generation with Sherpa. These changes read
``SHOWER_GENERATOR: None``, to switch off parton showering,
``FRAGMENTATION: None``, to do so for the hadronisation effects,
``MI_HANDLER: None``, to switch off multiparton interactions, and
``ME_QED: {ENABLED: false}``, to switch off resummed QED corrections
onto the :math:`Z \rightarrow e^- e^+` decay. Additionally, the
non-perturbative intrinsic transverse momentum may be wished to not be
taken into account, therefore set ``BEAM_REMNANTS: false``.

.. _Multijet merged event generation with Sherpa:

Multijet merged event generation with Sherpa
============================================

For a large fraction of LHC final states, the application of
reconstruction algorithms leads to the identification of several hard
jets. Calculations therefore need to describe as accurately as
possible both the hard jet production as well as the subsequent
evolution and the interplay of multiple such topologies. Several
scales determine the evolution of the event.

Various such merging schemes have been proposed: :cite:`Catani2001cc`,
:cite:`Lonnblad2001iq`, :cite:`Mangano2001xp`, :cite:`Krauss2002up`,
:cite:`Mangano2006rw`, :cite:`Lavesson2008ah`, :cite:`Hoeche2009rj`,
:cite:`Hamilton2009ne`, :cite:`Hamilton2010wh`, :cite:`Hoeche2010kg`,
:cite:`Lonnblad2011xx`, :cite:`Hoeche2012yf`, :cite:`Gehrmann2012yg`,
:cite:`Lonnblad2012ng`, :cite:`Lonnblad2012ix`.  Comparisons of the
older approaches can be found e.g. in :cite:`Hoche2006ph`,
:cite:`Alwall2007fs`. The currently most advanced treatment at
tree-level, detailed in :cite:`Hoeche2009rj`, :cite:`Hoeche2009xc`,
:cite:`Carli2009cg`, is implemented in Sherpa.

How to setup a multijet merged calculation is detailed in most
:ref:`Examples`, e.g. :ref:`LHC_WJets`, :ref:`LHC_ZJets` or
:ref:`TopsJets`.



.. _Running Sherpa with AMEGIC++:

Running Sherpa with AMEGIC++
============================

When Sherpa is run using the matrix element generator AMEGIC++, it is
necessary to run it twice. During the first run (the initialization
run) Feynman diagrams for the hard processes are constructed and
translated into helicity amplitudes. Furthermore suitable phase-space
mappings are produced. The amplitudes and corresponding integration
channels are written to disk as C++ source code, placed in a
subdirectory called ``Process``. The initialization run is started
using the standard Sherpa executable, as described in :ref:`Running
Sherpa`. The relevant command is

.. code-block:: shell-session

   $ <prefix>/bin/Sherpa

The initialization run stops with the message "New libraries
created. Please compile.", which is nothing but the request to carry
out the compilation and linking procedure for the generated
matrix-element libraries. The ``makelibs`` script, provided for this
purpose and created in the working directory, must be invoked by the
user (see ``./makelibs -h`` for help):

.. code-block:: shell-session

   $ ./makelibs

Note that the ``cmake`` tool has to be available for this step

.. index:: AMEGIC_LIBRARY_MODE

Another option is :kbd:`./makelibs -m`, which creates one library per
subprocess. This can be useful for very complex processes, in
particular if the default combined library generation fails due to a
limit on the number of command line arguments.  Note that this option
requires that Sherpa is run with ``AMEGIC_LIBRARY_MODE: 0`` (default:
1).

Afterwards Sherpa can be restarted using the same command as
before. In this run (the generation run) the cross sections of the
hard processes are evaluated. Simultaneously the integration over
phase space is optimized to arrive at an efficient event generation.

.. _Cross section determination:

***************************
Cross section determination
***************************

To determine the total cross section, in particular in the context of
multijet merging with Sherpa, the final output of the event generation
run should be used, e.g.

.. code-block:: console

   +----------------------------------------------------------------------------+
   | Nominal or variation name     XS [pb]      RelDev  AbsErr [pb]      RelErr |
   +----------------------------------------------------------------------------+
   | Nominal                       1739.79         0 %      171.304      9.84 % |
   | ME & PS: MUR=0.5 MUF=0.5      1635.61     -5.98 %      187.894     11.48 % |
   | ME & PS: MUR=2 MUF=2          2261.57     29.99 %      387.031     17.11 % |
   | [<results for other variations>]                                           |
   +----------------------------------------------------------------------------+

Note that the Monte Carlo error quoted for the total cross section is
determined during event generation. It, therefore, might differ
substantially from the errors quoted during the integration step, and
it can be reduced simply by generating more events.

In contrast to plain fixed order results, Sherpa's total cross section
in multijet merging setups (MEPS, MENLOPS, MEPS\@NLO) is composed of
values from various fixed order processes, namely those which are
combined by applying the multijet merging, see :ref:`Multijet merged
event generation with Sherpa`. In this context, it is important to
note:

**The higher multiplicity tree-level  cross sections determined during
the integration step are meaningless by themselves, only the inclusive
cross section printed at the end of  the event generation run is to be
used.**

**Sherpa total  cross sections  have leading  order accuracy  when the
generator is run  in LO merging mode (MEPS), in  NLO merging (MENLOPS,
MEPS\@NLO) mode they have NLO accuracy.**


Differential cross sections from single events
==============================================


To calculate the expectation value of an observable defined through a
series of cuts and requirements each event produced by Sherpa has to
be evaluated whether it meets the required criteria. The expectation
value is then given by

.. math::
   \langle O\rangle = \frac{1}{N_\text{trial}} \cdot \sum_i^n {w_i(\Phi_i) O(\Phi_i)}.

Therein the :math:`w_i(\Phi_i)` are the weight of the event with the
phase space configuration :math:`\Phi_i` and :math:`O(\Phi_i)` is the
value of the observable at this point. :math:`N_\text{trial} =
\sum_i^n n_{\text{trial,i}}` is the sum of number of trials
:math:`n_\text{trial,i}` of all events. A good cross check is to
reproduce the inclusive cross section as quoted by Sherpa (see above).

In case of unweighted events one might want to rescale the uniform
event weight to unity using ``w_norm``. The above equation then reads

.. math::
   \langle O \rangle = \frac{w_\text{norm}}{N_\text{trial}} \cdot \sum_i^n{\frac{w_i(\Phi_i)}{w_\text{norm} O(\Phi_i)}}

wherein :math:`\frac{w_i(\Phi_i)}{w_\text{norm}} = 1`, i.e. the sum simply
counts how many events pass the selection criteria of the
observable. If however, ``PartiallyUnweighted`` event weights or
``Enhance_Factor`` or ``Enhance_Observable`` are used, this is no
longer the case and the full form needs to be used.

All required quantities, :math:`w_i`, :math:`w_\text{norm}` and
:math:`n_{\text{trial},i}`, accompany each event and are written
e.g. into the HepMC output (cf. :ref:`Event output formats`).
