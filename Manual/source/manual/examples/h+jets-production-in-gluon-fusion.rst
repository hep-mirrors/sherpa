.. _LHC_HJets:

H+jets production in gluon fusion
=================================

This is an example setup for inclusive Higgs production through
gluon fusion at hadron colliders used in :cite:`Hoeche2014lxa`.
The inclusive process is calculated
at next-to-leading order accuracy matched to the parton shower using
the MC\@NLO prescription detailed in :cite:`Hoeche2011fd`. The next
few higher jet multiplicities, calculated at next-to-leading order as well,
are merged into the inclusive sample using the MEPS\@NLO method -- an
extension of the CKKW method to NLO -- as described in :cite:`Hoeche2012yf`
and :cite:`Gehrmann2012yg`. Finally, even higher multiplicities, calculated
at leading order, are merged on top of that. A few things to note are
detailed below the example.

.. literalinclude:: /../../Examples/H_in_GluonFusion/LHC_HJets/Sherpa.yaml
   :language: yaml

Things to notice:

* The example can be converted into a simple MENLOPS setup by
  replacing ``2->1-3`` with ``2->1``, or into an MEPS setup with
  ``2->0``, to study the effect of incorporating higher-order matrix
  elements.

* Providers of the one-loop matrix elements for the respective
  multiplicities are set using ``Loop_Generator``. For the two
  simplest cases Sherpa can provide it internally. Additionally, MCFM
  is interfaced for the H+2jet process, cf. :ref:`MCFM interface`.

* To enable the Higgs to decay to a pair of photons, for example, the
  hard decays are invoked.  For details on the hard decay handling and
  how to enable specific decay modes see :ref:`Hard decays`.
