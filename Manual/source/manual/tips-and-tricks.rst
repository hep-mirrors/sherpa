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

   $ cmake -DHepMC3_DIR=/path/to/hepmc3 -DRIVET_DIR=/path/to/rivet

(Note: Both paths are equal if you used the Rivet bootstrap script.)
In the case that the packages are installed in standard locations,
you can instead use ``-DSHERPA_ENABLE_HEPMC3=ON``
and ``-DSHERPA_ENABLE_RIVET=ON``, respectively.

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
the ``rivet-mkhtml`` script. If you want to employ custom Rivet analyses you might
need to set the corresponding Rivet path variable, for example via

.. code-block:: shell-session

   $ export RIVET_ANALYSIS_PATH=$RIVET_ANALYSIS_PATH:<path to custom analysis lib>


The `RIVET:` block can be used with further options especially suitable for detailed
studies. Adding `JETCONTS: 1` will create separate histograms split by jet multiplicity as
created by the hard process. `SPLITSH: 1` creates histograms split by soft and
hard events, and `SPLITPM: 1` creates histograms split by events with positive and
negative event weights. Finally, `SPLITCOREPROCS: 1` will split by different
processes if multiple ones are specified in the runcard.

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
following option:

.. code-block:: yaml

   $ cmake -DMCFM_DIR=/path/to/MCFM

Or, if MCFM is installed in a standard location:

.. code-block:: yaml

   $ cmake -DSHERPA_ENABLE_MCFM=ON

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

To check whether the poles of the dipole subtraction and the
interfaced one-loop matrix element cancel for each phase space point,
specify
``AMEGIC: { CHECK_POLES: true }`` and/or ``COMIX: { CHECK_POLES: true }``.

In the same way, the
finite contributions of the infrared subtraction and the one-loop
matrix element can be checked using ``CHECK_FINITE``, and the
Born matrix element via ``CHECK_BORN``.  The accuracy to which the
poles, finite parts and Born matrix elements are checked is set via
``CHECK_THRESHOLD``.
These three settings are only supported by Amegic
and are thus set using
``AMEGIC: { <PARAMETER>: <VALUE> }``,
where ``<VALUE>`` is ``false`` or ``true`` for ``CHECK_FINITE``/``CHECK_BORN``,
or a number specifying the desired accuracy for ``CHECK_THRESHOLD``.


.. _Scale variations:

*****************************
A posteriori scale variations
*****************************

There are several ways to compute the effects of changing the scales
and PDFs of any event produced by Sherpa. They can computed
explicitly, cf. :ref:`Explicit scale variations`, on-the-fly, cf.
:ref:`On-the-fly event weight variations` (restricted to multiplicative
factors), or reconstructed a posteriori. The latter method needs
plenty of additional information in the event record and is (depending
on the actual calculation) available in two formats:

.. contents::
   :local:

.. _A posteriori scale and PDF variations using the HepMC GenEvent Output:

A posteriori scale and PDF variations using the HepMC GenEvent Output
=====================================================================

Events generated in a LO, LOPS, NLO, NLOPS, MEPS\@LO, MEPS\@NLO or
MENLOPS calculation can be written out in the HepMC format including
all information to carry out arbitrary scale variations a
posteriori. For this feature HepMC of at least version 2.06 is
necessary and both ``HEPMC_USE_NAMED_WEIGHTS: true`` and
``HEPMC_EXTENDED_WEIGHTS: true`` have to enabled. Detailed
instructions on how to use this information to construct the new event
weight can be found here
`<https://sherpa.hepforge.org/doc/ScaleVariations-Sherpa-2.2.0.pdf>`_.

.. _A posteriori scale and PDF variations using the ROOT NTuple Output:

A posteriori scale and PDF variations using the ROOT NTuple Output
==================================================================

.. index:: USR_WGT_MODE

Events generated at fixed-order LO and NLO can be stored in ROOT
NTuples that allow arbitrary a posteriori scale and PDF variations,
see :ref:`Event output formats`. An example for writing and reading in
such ROOT NTuples can be found here: :ref:`NTuple production`.  The
internal ROOT Tree has the following Branches:

``id``
  Event ID to identify correlated real sub-events.

``nparticle``
  Number of outgoing partons.

``E/px/py/pz``
  Momentum components of the partons.

``kf``
  Parton PDG code.

``weight``
  Event weight, if sub-event is treated independently.

``weight2``
  Event weight, if correlated sub-events are treated as single event.

``me_wgt``
  ME weight (w/o PDF), corresponds to 'weight'.

``me_wgt2``
  ME weight (w/o PDF), corresponds to 'weight2'.

``id1``
  PDG code of incoming parton 1.

``id2``
  PDG code of incoming parton 2.

``fac_scale``
  Factorisation scale.

``ren_scale``
  Renormalisation scale.

``x1``
  Bjorken-x of incoming parton 1.

``x2``
  Bjorken-x of incoming parton 2.

``x1p``
  x' for I-piece of incoming parton 1.

``x2p``
  x' for I-piece of incoming parton 2.

``nuwgt``
  Number of additional ME weights for loops and integrated subtraction terms.

``usr_wgt[nuwgt]``
  Additional ME weights for loops and integrated subtraction terms.

Computing (differential) cross sections of real correction events with statistical errors
=========================================================================================

Real correction events and their counter-events from subtraction terms are
highly correlated and exhibit large cancellations. Although a treatment of
sub-events as independent events leads to the correct cross section the
statistical error would be greatly overestimated. In order to get a realistic
statistical error sub-events belonging to the same event must be combined
before added to the total cross section or a histogram bin of a differential
cross section. Since in general each sub-event comes with it's own set of four
momenta the following treatment becomes necessary:

#. An event here refers to a full real correction event that may
   contain several sub-events. All entries with the same id belong to
   the same event.  Step 2 has to be repeated for each event.

#. Each sub-event must be checked separately whether it passes
   possible phase space cuts. Then for each observable add up
   ``weight2`` of all sub-events that go into the same histogram
   bin. These sums :math:`x_{id}` are the quantities to enter the actual
   histogram.

#. To compute statistical errors each bin must store the sum over all
   :math:`x_{id}` and the sum over all :math:`x_{id}^2`. The cross section
   in the bin is given by :math:`\langle x\rangle = \frac{1}{N} \cdot
   \sum x_{id}`, where :math:`N` is the number of events (not
   sub-events). The :math:`1-\sigma` statistical error for the bin is
   :math:`\sqrt{ (\langle x^2\rangle-\langle x\rangle^2)/(N-1) }`

Note: The main difference between ``weight`` and ``weight2`` is that they
refer to a different counting of events. While ``weight`` corresponds to
each event entry (sub-event) counted separately, ``weight2`` counts events
as defined in step 1 of the above procedure. For NLO pieces other than the real
correction ``weight`` and ``weight2`` are identical.

Computation of cross sections with new PDF's
============================================

Born and real pieces
--------------------

Notation:

.. code-block:: text

   f_a(x_a) = PDF 1 applied on parton a, F_b(x_b) = PDF 2 applied on
   parton b.

The total cross section weight is given by:

.. code-block:: text

   weight = me_wgt f_a(x_a)F_b(x_b)

Loop piece and integrated subtraction terms
-------------------------------------------

The weights here have an explicit dependence on the renormalization
and factorization scales.

To take care of the renormalization scale dependence (other than via
``alpha_S``) the weight ``w_0`` is defined as


.. code-block:: text

   w_0 = me_wgt + usr_wgts[0] log((\mu_R^new)^2/(\mu_R^old)^2) +
   usr_wgts[1] 1/2 [log((\mu_R^new)^2/(\mu_R^old)^2)]^2

To address the factorization scale dependence the weights ``w_1,...,w_8``
are given by

.. code-block:: text

   w_i = usr_wgts[i+1] + usr_wgts[i+9] log((\mu_F^new)^2/(\mu_F^old)^2)

The full cross section weight can be calculated as

.. code-block:: text

   weight = w_0 f_a(x_a)F_b(x_b)
             + (f_a^1 w_1 + f_a^2 w_2 + f_a^3 w_3 + f_a^4 w_4) F_b(x_b)
             + (F_b^1 w_5 + F_b^2 w_6 + F_b^3 w_7 + F_b^4 w_8) f_a(x_a)

where

.. code-block:: text

   f_a^1 = f_a(x_a) (a=quark), \sum_q f_q(x_a) (a=gluon),
   f_a^2 = f_a(x_a/x'_a)/x'_a (a=quark), \sum_q f_q(x_a/x'_a)x'_a (a=gluon),
   f_a^3 = f_g(x_a),
   f_a^4 = f_g(x_a/x'_a)/x'_a

The scale dependence coefficients ``usr_wgts[0]`` and ``usr_wgts[1]``
are normally obtained from the finite part of the virtual correction
by removing renormalization terms and universal terms from dipole
subtraction.  This may be undesirable, especially when the loop
provider splits up the calculation of the virtual correction into
several pieces, like leading and sub-leading color. In this case the
loop provider should control the scale dependence coefficients, which
can be enforced with option :option:`USR_WGT_MODE: false`.

.. warning::

   The loop provider must support this option or the scale dependence
   coefficients will be invalid!
