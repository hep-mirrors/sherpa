.. _Tips and Tricks:

###############
Tips and tricks
###############

.. contents::
   :local:

.. _Shell completion:

****************
Shell completion
****************


Sherpa will install a file named
``$prefix/share/SHERPA-MC/sherpa-completion`` which contains tab completion
functionality for the bash shell. You simply have to
source it in your active
shell session by running

.. code-block:: shell-session

   $ .  $prefix/share/SHERPA-MC/sherpa-completion

and you will be able to tab-complete any parameters on a Sherpa
command line.

To permanently enable this feature in your bash shell, you'll have to add the
source command above to your :file:`~/.bashrc`.

.. _Rivet analyses:

**************
Rivet analyses
**************

.. index:: ANALYSIS_OUTPUT

Sherpa is equipped with an interface to the analysis tool `Rivet
<http://projects.hepforge.org/rivet/>`_. To enable it, Rivet and
`HepMC <http://lcgapp.cern.ch/project/simu/HepMC/>`_ have to be
installed (e.g. using the Rivet bootstrap script) and your Sherpa
compilation has to be configured with the following options:

.. code-block:: shell-session

   $ cmake -DSHERPA_ENABLE_HEPMC3=ON -DHepMC3_DIR=/path/to/hepmc3 \
           -DSHERPA_ENABLE_RIVET=ON -DRIVET_DIR=/path/to/rivet

(Note: Both paths are equal if you used the Rivet bootstrap script.)

To use the interface, you need to enable it using the
:option:`ANALYSIS` option and to configure it it using the
:option:`RIVET` settings group as follows:

.. code-block:: yaml

   ANALYSIS: Rivet
   RIVET:
     --analyses:
       - D0_2008_S7662670
       - CDF_2007_S7057202
       - D0_2004_S5992206
       - CDF_2008_S7828950

The analyses list specifies which Rivet analyses to run and the
histogram output file can be changed with the normal ``ANALYSIS_OUTPUT``
switch.

Further Rivet options (especially for Rivet v3) can be passed through
the interface. The following ones are currently implemented:

.. code-block:: yaml

   ANALYSIS: Rivet
   RIVET:
     --analyses:
       - MC_ZINC
     --ignore-beams: 1
     --skip-weights: 0
     --match_weights: ".*MUR.*"
     --unmatch-weights: "NTrials"
     --nominal-weight: "Weight"
     --weight-cap: 100.0
     --nlo-smearing: 0.1

You can also use ``rivet-mkhtml`` (distributed with Rivet) to create
plot webpages from Rivet's output files:

.. code-block:: shell-session

   $ source /path/to/rivetenv.sh   # see below
   $ rivet-mkhtml -o output/ file1.yoda [file2.yoda, ...]
   $ firefox output/index.html &

If your Rivet installation is not in a standard location, the bootstrap script
should have created a :file:`rivetenv.sh` which you have to source before running
the ``rivet-mkhtml`` script.

.. _HZTool analyses:

***************
HZTool analyses
***************

.. index:: ANALYSIS_OUTPUT

Sherpa is equipped with an interface to the analysis tool `HZTool
<http://projects.hepforge.org/hztool/>`_. To enable it, HZTool and
`CERNLIB <http://cernlib.web.cern.ch/>`_ have to be installed and your
Sherpa compilation has to be configured with the following options:

.. code-block:: shell-session

   $ cmake -DSHERPA_ENABLE_HZTOOL=ON -DHZTOOL_DIR=/path/to/hztool \
           -DCERNLIB_DIR=/path/to/cernlib -DHEPEVT_CB_SIZE=4000

Note that an example CERNLIB installation bootstrap script is provided
in ``AddOns/HZTool/start_cern_64bit``. Note that this script is only
provided for convenience, we will not provide support if it is not
working as expected.

To use the interface, enable it using the :option:`ANALYSIS` and
configure it using the :option:`HZTool` settings group:

.. code-block:: yaml

   ANALYSIS: HZTool
   HZTOOL:
     HISTO_NAME: output.hbook
     HZ_ENABLE:
     - hz00145
     - hz01073
     - hz02079
     - hz03160

The ``HZ_ENABLE`` list specifies which HZTool analyses to run.  The
histogram output directory can be changed using the
``ANALYSIS_OUTPUT`` switch, while ``HZTOOL:HISTO_NAME`` specifies the
hbook output file.

.. _MCFM interface:

**************
MCFM interface
**************

.. index:: Loop_Generator

Sherpa is equipped with an interface to the NLO library of `MCFM
<http://mcfm.fnal.gov/>`_ for decdicated processes.  To enable it,
MCFM has to be installed and compiled into a single library,
libMCFM.a. To this end, an installation script is provided in
``AddOns/MCFM/install_mcfm.sh``. Please note, due to some process
specific changes that are made by the installation script to the MCFM
code, only few selected processes of MCFM-6.3 are available through
the interface.

Finally, your Sherpa compilation has to be configured with the
following options:

.. code-block:: yaml

   $ cmake -DSHERPA_ENABLE_MCFM=ON -DMCFM_DIR=/path/to/MCFM

To use the interface, specify

.. code-block:: yaml

   Loop_Generator: MCFM

in the process section of the run card and add it to the list of
generators in :ref:`ME_GENERATORS`. Of course, MCFM's process.DAT file
has to be copied to the current run directory.

.. _Debugging a crashing/stalled event:

**********************************
Debugging a crashing/stalled event
**********************************

Crashing events
===============

If an event crashes, Sherpa tries to obtain all the information needed to
reproduce that event and writes it out into a directory named

.. code-block:: text

  Status__<date>_<time>

If you are a Sherpa user and want to report this crash to the Sherpa
team, please attach a tarball of this directory to your email. This
allows us to reproduce your crashed event and debug it.

To debug it yourself, you can follow these steps (Only do this if you
are a Sherpa developer, or want to debug a problem in an addon library
created by yourself):

* Copy the random seed out of the status directory into your run path:

  .. code-block:: shell-session

     $ cp  Status__<date>_<time>/random.dat  ./

* Run your normal Sherpa commandline with an additional parameter:

  .. code-block:: shell-session

     $ Sherpa [...] 'STATUS_PATH: ./'

  Sherpa will then read in your random seed from "./random.dat" and
  generate events from it.

* Ideally, the first event will lead to the crash you saw earlier, and
  you can now turn on debugging output to find out more about the
  details of that event and test code changes to fix it:

  .. code-block:: shell-session

     $ Sherpa [...] --output 15 'STATUS_PATH: ./'

Stalled events
==============

If event generation seems to stall, you first have to find out
the number of the current event. For that you would terminate the stalled
Sherpa process (using Ctrl-c) and check in its final output for the number
of generated events.
Now you can request Sherpa to write out the random seed for the event before the
stalled one:

.. code-block:: shell-session

   $ Sherpa [...] --events <#events - 1> 'SAVE_STATUS: Status/'

(Replace ``<#events - 1>`` using the number you figured out earlier.)

The created status directory can either be sent to the Sherpa
developers, or be used in the same steps as above to reproduce that
event and debug it.

.. _Versioned installation:

**********************
Versioned installation
**********************

If you want to install different Sherpa versions into the same prefix
(e.g. `/usr/local`), you have to enable versioning of the installed
directories by using the configure option ``-DSHERPA_ENABLE_VERSIONING=ON``.
Optionally you can even pass an argument to this parameter of what you
want the version tag to look like.

.. _NLO calculations:

****************
NLO calculations
****************

.. contents::
   :local:

.. _Choosing DIPOLES ALPHA:

Choosing DIPOLES ALPHA
======================

A variation of the parameter ``DIPOLES:ALPHA`` (see :ref:`Dipole
subtraction`) changes the contribution from the real (subtracted)
piece (``RS``) and the integrated subtraction terms (``I``), keeping
their sum constant.  Varying this parameter provides a nice check of
the consistency of the subtraction procedure and it allows to optimize
the integration performance of the real correction. This piece has the
most complicated momentum phase space and is often the most time
consuming part of the NLO calculation.  The optimal choice depends on
the specific setup and can be determined best by trial.

Hints to find a good value:

* The smaller ``DIPOLES:ALPHA`` is the less dipole term have to be
  calculated, thus the less time the evaluation/phase space point
  takes.

* Too small choices lead to large cancellations between the ``RS``
  and the ``I`` parts and thus to large statistical errors.

* For very simple processes (with only a total of two partons in the
  initial and the final state of the born process) the best choice is
  typically ``DIPOLES: {ALPHA: 1``}.  The more complicated a process
  is the smaller ``DIPOLES:ALPHA`` should be (e.g. with 5 partons the
  best choice is typically around 0.01).

* A good choice is typically such that the cross section from the
  ``RS`` piece is significantly positive but not much larger than
  the born cross section.

.. _Integrating complicated Loop-ME:

Integrating complicated Loop-ME
===============================

For complicated processes the evaluation of one-loop matrix elements
can be very time consuming. The generation time of a fully optimized
integration grid can become prohibitively long. Rather than using a
poorly optimized grid in this case it is more advisable to use a grid
optimized with either the born matrix elements or the born matrix
elements and the finite part of the integrated subtraction terms only,
working under the assumption that the distributions in phase space are
rather similar.

This can be done by one of the following methods:

#. Employ a dummy virtual (requires no computing time, returns a
   finite value as its result) to optimise the grid. This only works
   if ``V`` is not the only ``NLO_Part`` specified.

   #. During integration set the ``Loop_Generator`` to ``Dummy``. The
      grid will then be optimised to the phase space distribution of
      the sum of the Born matrix element and the finite part of the
      integrated subtraction term, plus a finite value from ``Dummy``.

      .. note::

         The cross section displayed during integration will also
         correspond to these contributions.

   #. During event generation reset ``Loop_Generator`` to your
      generator supplying the virtual correction. The events generated
      then carry the correct event weight.

#. Suppress the evaluation of the virtual and/or the integrated
   subtraction terms. This only works if Amegic is used as the matrix
   element generator for the ``BVI`` pieces and ``V`` is not the only
   ``NLO_Part`` specified.


   #. During integration add ``AMEGIC: { NLO_BVI_MODE: <num> }`` to
      your configuration. ``<num>`` takes the following values:
      ``1``-``B``, ``2``-``I``, and ``4``-``V``. The values are
      additive, i.e.  ``3``-``BI``.


      .. note::

         The cross section displayed during integration will match the parts
         selected by ``NLO_BVI_MODE``.

   #. During event generation remove the switch again and the events
      will carry the correct weight.


.. note::

   this will not work for the ``RS`` piece!

.. _Avoiding misbinning effects:

Avoiding misbinning effects
===========================

Close to the infrared limit, the real emission matrix element and
corresponding subtraction events exhibit large cancellations. If the
(minor) kinematics difference of the events happens to cross a
parton-level cut or analysis histogram bin boundary, then large
spurious spikes can appear.

These can be smoothed to some extend by shifting the weight from the
subtraction kinematic to the real-emission kinematic if the dipole
measure alpha is below a given threshold. The fraction of the shifted
weight is inversely proportional to the dipole measure, such that the
final real-emission and subtraction weights are calculated as:

.. code-block:: perl

   w_r -> w_r + sum_i [1-x(alpha_i)] w_{s,i}
   foreach i: w_{s,i} -> x(alpha_i) w_{s,i}

with the function :math:`x(\alpha)=(\frac{\alpha}{|\alpha_0|})^n` for
:math:`\alpha<\alpha_0` and :math:`1` otherwise.

The threshold can be set by the parameter
``NLO_SMEAR_THRESHOLD: <alpha_0>`` and the functional form of
alpha and thus interpretation of the threshold can be chosen by its
sign (positive: relative dipole kT in GeV, negative: dipole alpha).
In addition, the exponent n can be set by ``NLO_SMEAR_POWER: <n>``.

.. _Enforcing the renormalization scheme:

Enforcing the renormalization scheme
====================================

.. index:: LOOP_ME_INIT

Sherpa takes information about the renormalization scheme from the
loop ME generator.  The default scheme is MSbar, and this is assumed
if no loop ME is provided, for example when integrated subtraction
terms are computed by themselves.  This can lead to inconsistencies
when combining event samples, which may be avoided by setting
``AMEGIC: { LOOP_ME_INIT: 1 }``.

.. _Checking the pole cancellation:

Checking the pole cancellation
==============================

.. index:: CHECK_BORN
.. index:: CHECK_FINITE
.. index:: CHECK_POLES
.. index:: CHECK_THRESHOLD

The following options are all sub-settings for :option:`AMEGIC` and
can be specified as follows:

.. code-block:: yaml

   AMEGIC:
     <option>: <value>
     ...

To check whether the poles of the dipole subtraction and the
interfaced one-loop matrix element cancel phase space point by phase
space point ``CHECK_POLES: 1`` can be specified.  In the same way, the
finite contributions of the infrared subtraction and the one-loop
matrix element can be checked by setting ``CHECK_FINITE: 1``, and the
Born matrix element via ``CHECK_BORN: 1``.  The accuracy to which the
poles, finite parts and Born matrix elements are checked is set via
``CHECK_THRESHOLD: <accu>``.
