.. _Examples:

########
Examples
########

Some example set-ups are included in Sherpa, in the
``<prefix>/share/SHERPA-MC/Examples/`` directory. These may be useful
to new users to practice with, or as templates for creating your own
Sherpa run-cards. In this section, we will look at some of the main
features of these examples.

.. contents::
   :local:
   :depth: 1

.. _VJets:

******************************
Vector boson + jets production
******************************

.. contents::
   :local:

To change any of the following LHC examples to production at different
collider energies or beam types, e.g. proton anti-proton at the Tevatron,
simply change the beam settings accordingly:

.. code-block:: yaml

   BEAMS: [2212, -2212]
   BEAM_ENERGIES: 980

.. include:: ./examples/w+jets-production.rst
.. include:: ./examples/z+jets-production.rst
.. include:: ./examples/bb-production.rst

.. _Jets:

**************
Jet production
**************


.. contents::
   :local:

.. include:: ./examples/jet-production.rst
.. include:: ./examples/jets-at-lepton-colliders.rst

.. _HJets:

*****************************
Higgs boson + jets production
*****************************

.. contents::
   :local:

.. include:: ./examples/h-production-in-gluon-fusion-with-interference-effects.rst
.. include:: ./examples/h+jets-production-in-gluon-fusion.rst
.. include:: ./examples/h+jets-production-in-gluon-fusion-with-finite-top-mass-effects.rst
.. include:: ./examples/h+jets-production-in-associated-production.rst
.. include:: ./examples/associated-t-anti-t-h-production-at-the-lhc.rst

.. _TopsJets:

**********************************
Top quark (pair) + jets production
**********************************


.. contents::
   :local:

.. include:: ./examples/Top-quark-(pair)-jets-production.rst
.. include:: ./examples/Production-of-a-top-quark-pair-in-association-with-a-W-boson.rst

.. _SingleTopChannels:

************************************************
Single-top production in the s, t and tW channel
************************************************


In this section, examples for single-top production in three different
channels are described.  For the channel definitions and a validation
of these setups, see :cite:`Bothmann2017jfv`.

.. contents::
   :local:

.. include:: ./examples/t-channel-single-top-production.rst
.. include:: ./examples/s-channel-single-top-production.rst
.. include:: ./examples/tw-channel-single-top-production.rst

.. _VVJets:

************************************
Vector boson pairs + jets production
************************************

.. contents::
   :local:

.. include:: ./examples/vector-boson-pairs-jets-production.rst

.. _MSSM/UFO:

**************************************
Event generation in the MSSM using UFO
**************************************

.. include:: ./examples/event-generation-in-the-mssm-using-ufo.rst

.. _DIS:

*************************
Deep-inelastic scattering
*************************

.. contents::
   :local:

.. include:: ./examples/deep-inelastic-scattering.rst

.. _FONLO:

**********************************************
Fixed-order next-to-leading order calculations
**********************************************

.. contents::
   :local:

.. include:: ./examples/production-of-ntuples.rst
.. include:: ./examples/minlo.rst

.. _SoftQCD:

*****************************************
Soft QCD: Minimum Bias and Cross Sections
*****************************************

.. contents::
   :local:

.. include:: ./examples/calculation-of-inclusive-cross-sections.rst
.. include:: ./examples/simulation-of-minimum-bias-events.rst

.. _BFactories:

******************************************
Setups for event production at B-factories
******************************************

.. contents::
   :local:

.. include:: ./examples/qcd-continuum.rst

.. _MEvalues:

*********************************************************************
Calculating matrix element values for externally given configurations
*********************************************************************

.. contents::
   :local:

.. include:: ./examples/computing-matrix-elements-for-individual-phase-space-points-using-the-python-interface.rst

.. _APIexamples:

**************************
Using the Python interface
**************************

.. contents::
   :local:

.. include:: ./examples/generate-events-using-scripts.rst
.. include:: ./examples/generate-events-with-mpi-using-scripts.rst
