.. _NTuple production:

Production of NTuples
=====================

Root NTuples are a convenient way to store the result of cumbersome
fixed-order calculations in order to perform multiple analyses.  This
example shows how to generate such NTuples and reweighted them in
order to change factorisation and renormalisation scales.  Note that
in order to use this setup, Sherpa must be configured with option
:option:`-DSHERPA_ENABLE_ROOT`, see :ref:`Event output
formats`.  If Sherpa has not been configured with Rivet analysis
support, please disable the analysis using :option:`-a0` on the
command line, see :ref:`Command line`.

When using NTuples, one needs to bear in mind that every calculation
involving jets in the final state is exclusive in the sense that a
lower cut-off on the jet transverse momenta must be imposed.  It is
therefore necessary to check whether the event sample stored in the
NTuple is sufficiently inclusive before using it. Similar remarks
apply when photons are present in the NLO calculation or when cuts on
leptons have been applied at generation level to increase efficiency.
Every NTuple should therefore be accompanied by an appropriate
documentation.

NTuple compression can be customized using the parameter
:option:`ROOTNTUPLE_COMPRESSION`, which is used to call
``TFile::SetCompressionSettings``. For a detailed documentation of
available options, see `<http://root.cern.ch>`_

This example will generate NTuples for the process pp->lvj, where `l`
is an electron or positron, and `v` is an electron (anti-)neutrino.
We identify parton-level jets using the anti-k_T algorithm with
`R=0.4` :cite:`Cacciari2008gp`. We require the transverse momentum of
these jets to be larger than 20 GeV. No other cuts are applied at
generation level.

.. literalinclude:: /../../Examples/FixedOrder_NLO/NTuples/Sherpa.B-like.yaml
   :language: yaml

Things to notice:

* NTuple production is enabled by :option:`EVENT_OUTPUT:
  Root[NTuple_B-like]`, see :ref:`Event output formats`.

* The scale used is defined as in :cite:`Berger2009ep`.

* :option:`EW_SCHEME: 0` and :option:`WIDTH_SCHEME: Fixed` are used to
  set the value of the weak mixing angle to 0.23, consistent with EW
  precision measurements.

* :option:`DIPOLES:ALPHA: 0.03` is used to limit the active phase
  space of dipole subtractions.

* :option:`13:Massive: true` and :option:`15:Massive: 1` are used to
  limit the number of active lepton flavours to electron and positron.

* The option :option:`USE_HEPMC_SHORT: 1` is used in the Rivet
  analysis section as the events produced by Sherpa are not at
  particle level.

NTuple production
-----------------

Start Sherpa using the command line

.. code-block:: shell-session

   $ Sherpa Sherpa.B-like.yaml

Sherpa will first create source code for its matrix-element calculations.
This process will stop with a message instructing you to compile.
Do so by running

.. code-block:: shell-session

   $ ./makelibs -j4

Launch Sherpa again, using

.. code-block:: shell-session

   $ Sherpa Sherpa.B-like.yaml

Sherpa will then compute the Born, virtual and integrated subtraction
contribution to the NLO cross section and generate events. These
events are analysed using the Rivet library and stored in a Root
NTuple file called :kbd:`NTuple_B-like.root`.  We will use this NTuple
later to compute an NLO uncertainty band.

The real-emission contribution, including subtraction terms, to the
NLO cross section is computed using

.. code-block:: shell-session

   $ Sherpa Sherpa.R-like.yaml

Events are generated, analysed by Rivet and stored in the Root NTuple
file :kbd:`NTuple_R-like.root`.

The two analyses of events with Born-like and real-emission-like
kinematics need to be merged, which can be achieved using scripts like
:kbd:`yodamerge`.  The result can then be plotted and displayed.

Usage of NTuples in Sherpa
--------------------------

Next we will compute the NLO uncertainty band using Sherpa.  To this
end, we make use of the Root NTuples generated in the previous steps.
Note that the setup files for reweighting are almost identical to
those for generating the NTuples. We have simply replaced
:option:`EVENT_OUTPUT` by :option:`EVENT_INPUT`.

We re-evaluate the events with the scale variation as defined in the
:kbd:`Reweight` configuration files:

.. code-block:: shell-session

   $ Sherpa Sherpa.Reweight.B-like.yaml
   $ Sherpa Sherpa.Reweight.R-like.yaml

The contributions can again be combined using :kbd:`yodamerge`.
