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

.. include:: ./examples/V_plus_Jets/LHC_WJets/README.rst
.. include:: ./examples/V_plus_Jets/LHC_ZJets/README.rst
.. include:: ./examples/V_plus_Bs/LHC_Fusing/README.rst
.. include:: ./examples/V_plus_Bs/LHC_Zbb/README.rst

.. _Jets:

**************
Jet production
**************


.. contents::
   :local:

.. include:: ./examples/Jets_at_HadronColliders/LHC_Jets_MCatNLO/README.rst
.. include:: ./examples/Jets_at_LeptonColliders/LEP_Jets/README.rst

.. _HJets:

*****************************
Higgs boson + jets production
*****************************

.. contents::
   :local:

.. include:: ./examples/H_in_GluonFusion/LHC_HJets/README.rst
.. include:: ./examples/H_in_GluonFusion/LHC_HJets_Finite_MTop/README.rst
.. include:: ./examples/H_in_AssociatedProduction/LHC_WHJets/README.rst
.. include:: ./examples/H_in_TTBar/LHC_TTH_MCatNLO/README.rst

.. _TopsJets:

**********************************
Top quark (pair) + jets production
**********************************


.. contents::
   :local:

.. include:: ./examples/Tops_plus_Jets/LHC_Tops/README.rst
.. include:: ./examples/Tops_plus_V/LHC_TTW/README.rst

.. _SingleTopChannels:

************************************************
Single-top production in the s, t and tW channel
************************************************


In this section, examples for single-top production in three different
channels are described.  For the channel definitions and a validation
of these setups, see :cite:`Bothmann2017jfv`.

.. contents::
   :local:

.. include:: ./examples/SingleTop_Channels/README.rst

.. _VVJets:

************************************
Vector boson pairs + jets production
************************************

.. contents::
   :local:

.. include:: ./examples/VV_plus_Jets/LHC_2l2nuJets/README.rst
.. include:: ./examples/VV_plus_Jets/LHC_4lJets/README.rst
.. include:: ./examples/VV_plus_Jets/LHC_2l2nu2jJets_SameSign/README.rst
.. include:: ./examples/VV_plus_Jets/LHC_WqqZnunuJets/README.rst
.. include:: ./examples/VV_plus_Jets/LHC_2l2nuJets_GluonInitiated/README.rst
.. include:: ./examples/VV_plus_Jets/LHC_4lJets_GluonInitiated/README.rst
.. include:: ./examples/VV_plus_Jets/LHC_WZ_polarization_at_nLO_PS/README.rst
.. include:: ./examples/VV_plus_Jets/LHC_ssWW_polarization/README.rst

.. _BSMexamples:

**********
BSM setups
**********

.. contents::
   :local:

.. include:: ./examples/BSM/UFO_SMEFT/README.rst
.. include:: ./examples/BSM/UFO_MSSM/README.rst

.. _DIS:

*************************
Deep-inelastic scattering
*************************

.. contents::
   :local:

.. include:: ./examples/Jets_in_DIS/HERA/README.rst

.. _FONLO:

**********************************************
Fixed-order next-to-leading order calculations
**********************************************

.. contents::
   :local:

.. include:: ./examples/FixedOrder/MINLO_ppenuj/README.rst

.. _SoftQCD:

*****************************************
Soft QCD: Minimum Bias and Cross Sections
*****************************************

.. contents::
   :local:

.. include:: ./examples/Soft_QCD/LHC_7TeV_MinBias/README.rst

.. _BFactories:

******************************************
Setups for event production at B-factories
******************************************

.. contents::
   :local:

.. include:: ./examples/BFactory/Belle_QCDBackground/README.rst

.. _MEvalues:

*********************************************************************
Calculating matrix element values for externally given configurations
*********************************************************************

.. contents::
   :local:

.. include:: ./examples/API/ME2-Python/README.rst

.. _APIexamples:

**************************
Using the Python interface
**************************

.. contents::
   :local:

.. include:: ./examples/API/Events/README.rst
.. include:: ./examples/API/MPIEvents/README.rst

.. _UserhookExamples:

***************************************
Custom event processing with user hooks
***************************************

.. contents::
   :local:

.. include:: ./examples/Userhook/README.rst
