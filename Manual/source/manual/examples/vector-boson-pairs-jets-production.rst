.. _Dilepton missing energy and jets production:

Dilepton, missing energy and jets production
============================================

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_2l2nuJets/Sherpa.yaml
   :language: yaml

.. _Dilepton missing energy and jets production (gluon initiated):

Dilepton, missing energy and jets production (gluon initiated)
==============================================================

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_2l2nuJets_GluonInitiated/Sherpa.yaml
   :language: yaml

.. _Four lepton and jets production:

Four lepton and jets production
===============================

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_4lJets/Sherpa.yaml
   :language: yaml

.. _Four lepton and jets production (gluon initiated):

Four lepton and jets production (gluon initiated)
=================================================

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_4lJets_GluonInitiated/Sherpa.yaml
   :language: yaml

.. _WZ production with jets production:

WZ production with jets production
==================================

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_WqqZnunuJets/Sherpa.yaml
   :language: yaml

.. _Same sign dilepton missing energy and jets production:

Same sign dilepton, missing energy and jets production
======================================================

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_2l2nu2jJets_SameSign/Sherpa.yaml
   :language: yaml

.. _Polarized same-sign W-pair production in association with two jets:

Polarized same-sign :math:`\mathrm{W}^+` boson pair production in association with two jets at LO+PS
====================================================================================================
This is an example for the simulation of polarized cross sections for pure electroweak same-sign :math:`\mathrm{W}^+` boson pair production in
association with two jets at LO+PS.

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_ssWW_polarization/Sherpa.yaml
   :language: yaml

Things to notice:

* The ``COMIX_DEFAULT_GAUGE`` is set to a special value to get the polarization vectors mentioned in section
  :ref:`Simulation of polarized cross sections for intermediate particles` .
* The width of the :math:`\mathrm{W}^\pm` boson is set to zero to retain gauge invariance since it is considered as stable during the hard
  scattering process (``PROCESSES``). To preserve SU(2) Ward Identities also the width of the Z boson need to be set to
  zero then. 
* The process is divided into the production of the vector bosons (``PROCESSES``) and their decays (``HARD_DECAYS``), all matrix elements are 
  calculated with on-shell vector bosons (narrow-width approximation). ``Mass_Smearing`` and ``SPIN_CORRELATIONS`` are enabled per default.

Polarized :math:`\mathrm{W}^+` Z boson pair production -pair production at nLO+PS
=================================================================================
This is an example for the simulation of polarized cross sections for pure electroweak :math:`\mathrm{W}^+` Z boson pair production 
at nLO+PS. The resulting unpolarized cross section contains all NLO QCD corrections while for polarized
cross sections the effect of virtual corrections on the polarization fractions is neglected. 

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_WZ_polarization_at_nLO_PS/Sherpa.yaml
   :language: yaml

