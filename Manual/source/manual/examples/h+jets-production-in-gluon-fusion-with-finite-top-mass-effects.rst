.. _LHC_HJets_FiniteMt:

H+jets production in gluon fusion with finite top mass effects
==============================================================

This is example is similar to :ref:`LHC_HJets` but with finite top quark
mass taken into account as described in :cite:`Buschmann2014sia` for
all merged jet multiplicities. Mass effects in the virtual corrections are
treated in an approximate way. In case of the tree-level contributions,
including real emission corrections, no approximations are made
concerning the mass effects.

.. literalinclude:: /../../Examples/H_in_GluonFusion/LHC_HJets_Finite_MTop/Sherpa.yaml
   :language: yaml

Things to notice:

* One-loop matrix elements from OpenLoops :cite:`Cascioli2011va`
  are used in order to correct for top mass effects. Sherpa must therefore
  be compiled with OpenLoops support to run this example. Also, the
  OpenLoops process libraries listed in the run card must be installed.

* The maximum jet multiplicities that can be merged in this setup
  are limited by the availability of loop matrix elements used to correct
  for finite top mass effects.

* The comments in :ref:`LHC_HJets` apply here as well.
