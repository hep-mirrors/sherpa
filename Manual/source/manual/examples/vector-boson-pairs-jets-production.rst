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

Polarized same-sign W-pair production in association with two jets
======================================================
This is an example for the simulation of polarized cross sections for pure electroweak same-sign W-pair production in
association with two jets at leading order.

.. literalinclude:: /../../Examples/VV_plus_Jets/LHC_ssWW_polarization/Sherpa.yaml
   :language: yaml

Things to notice:

* The ``COMIX_DEFAULT_GAUGE`` is set to a special value to get the polarization vectors mentioned in section
  :ref:`Simulation of polarized cross sections for intermediate particles` .
* The width of the W is set to zero to retain gauge invariance since it is considered as stable during the hard
  scattering process (``PROCESSES``).
* The process is divided into the production of the vector bosons (``PROCESSES``) and their decays (``HARD_DECAYS``) so
  that the narrow-width-approximation is applied. Therefore, all matrix elements are calculated with on-shell vector
  bosons.
  The ``Mass_Smearing`` setting ensures that the actual off-shell kinematic is recovered. It is enabled per default and
  is only given here to show the full approximation.
