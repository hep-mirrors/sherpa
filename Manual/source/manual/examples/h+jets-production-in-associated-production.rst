.. _LHC_VHJets:

H+jets production in associated production
==========================================

This section collects example setups for Higgs boson production in
association with vector bosons

.. contents::
   :local:

.. _LHC_WHJets:

Higgs production in association with W bosons and jets
------------------------------------------------------

This is an example setup for Higgs boson production in association with
a W boson and jets, as used in :cite:`Hoeche2014rya`. It uses the
MEPS\@NLO method to merge pp->WH and pp->WHj at next-to-leading order
accuracy and adds pp->WHjj at leading order. The Higgs boson is decayed
to W-pairs and all W decay channels resulting in electrons or muons are
accounted for, including those with intermediate taus.

.. literalinclude:: /../../Examples/H_in_AssociatedProduction/LHC_WHJets/Sherpa.yaml
   :language: yaml

Things to notice:

* Two custom particle container, cf. :ref:`Particle containers`, have
  been declared, facilitating the process declaration.

* As the bottom quark is treated as being massless by default, a five
  flavour calculation is performed. The particle container ensures
  that no external bottom quarks, however, are considered to resolve
  the overlap with single top and top pair processes.

* OpenLoops :cite:`Cascioli2011va` is used as the provider of the
  one-loop matrix elements.

* To enable the decays of the Higgs, W boson and tau lepton the hard
  decay handler is invoked. For details on the hard decay handling and
  how to enable specific decay modes see :ref:`Hard decays`.

.. _LHC_ZHJets:

Higgs production in association with Z bosons and jets
------------------------------------------------------

This is an example setup for Higgs boson production in association with
a Z boson and jets, as used in :cite:`Hoeche2014rya`. It uses the
MEPS\@NLO method to merge pp->ZH and pp->ZHj at next-to-leading order
accuracy and adds pp->ZHjj at leading order. The Higgs boson is decayed
to W-pairs. All W and Z bosons are allowed to decay into electrons, muons
or tau leptons. The tau leptons are then allowed to decay into all possible
partonic channels, leptonic and hadronic, to allow for all possible
trilepton signatures, unavoidably producing two and four lepton events
as well.

.. literalinclude:: /../../Examples/H_in_AssociatedProduction/LHC_ZHJets/Sherpa.yaml
   :language: yaml

Things to notice:

* A custom particle container, cf. :ref:`Particle containers`, has
  been declared, facilitating the process declaration.

* As the bottom quark is treated as being massless by default, a five flavour
  calculation is performed. The particle container ensures that no external
  bottom quarks, however, are considered to resolve the overlap with single
  top and top pair processes.

* OpenLoops :cite:`Cascioli2011va` is used as the provider of the one-loop
  matrix elements.

* To enable the decays of the Higgs, W and Z bosons and tau lepton the hard
  decay handler is invoked. For details on the hard decay handling and how
  to enable specific decay modes see :ref:`Hard decays`.

.. _LHC_Hll_MC@NLO:

Higgs production in association with lepton pairs
-------------------------------------------------

This is an example setup for Higgs boson production in association with
an electron-positron pair using the MC\@NLO technique. The Higgs boson
is decayed to b-quark pairs. Contrary to the previous examples this setup
does not use on-shell intermediate vector bosons in its matrix element
calculation.

.. literalinclude:: /../../Examples/H_in_AssociatedProduction/LHC_Hll_MCatNLO/Sherpa.yaml
   :language: yaml

Things to notice:

* The central scale is set to the
  invariant mass of the Higgs boson and the lepton pair.

* As the bottom quark is set to be treated massively, a four flavour
  calculation is performed.

* OpenLoops :cite:`Cascioli2011va` is used as the provider of the one-loop
  matrix elements.

* To enable the decays of the Higgs the hard
  decay handler is invoked. For details on the hard decay handling and how
  to enable specific decay modes see :ref:`Hard decays`.
