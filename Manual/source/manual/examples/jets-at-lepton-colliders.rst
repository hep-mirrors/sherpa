.. _LEP_Jets:

Jets at lepton colliders
========================

.. contents::
   :local:

This section contains two setups to describe jet production at LEP I,
either through multijet merging at leading order accuracy or at
next-to-leading order accuracy.

.. _LEP_Jets_LOmerging:

MEPS setup for ee->jets
-----------------------

This example shows a LEP set-up, with electrons and positrons
colliding at a centre of mass energy of 91.2 GeV.

.. literalinclude:: /../../Examples/Jets_at_LeptonColliders/LEP_Jets/Sherpa.yaml
   :language: yaml

Things to notice:

* The running of alpha_s is set to leading order and the value of
  alpha_s at the Z-mass is set.

* Note that initial-state radiation is enabled by default. See
  :ref:`ISR Parameters` on how to disable it if you want to evaluate
  the (unphysical) case where the energy for the incoming leptons is
  fixed.

.. _LEP_Jets_NLOmerging:

MEPS\@NLO setup for ee->jets
----------------------------

This example expands upon the above setup, elevating its description
of hard jet production to next-to-leading order.

.. literalinclude:: /../../Examples/Jets_at_LeptonColliders/LEP_Jets/Sherpa.NLO.yaml
   :language: yaml

Things to notice:

* the b-quark mass has been enabled for the matrix element calculation
  (the default is massless) because it is not negligible for LEP energies

* the b b-bar and b b b-bar b-bar processes are specified separately because the
  :option:`93` particle container contains only partons set massless in the matrix
  element calculation, see :ref:`Particle containers`.

* model parameters can be modified in the config file; in this example, the
  value of alpha_s at the Z mass is set.
