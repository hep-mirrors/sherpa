.. _LHC_ZJets:

Z+jets production
=================

This is an example setup for inclusive Z production at hadron
colliders.  The inclusive process is calculated at next-to-leading
order accuracy matched to the parton shower using the MC\@NLO
prescription detailed in :cite:`Hoeche2011fd`. The next few higher jet
multiplicities, calculated at next-to-leading order as well, are
merged into the inclusive sample using the MEPS\@NLO method -- an
extension of the CKKW method to NLO -- as described in
:cite:`Hoeche2012yf` and :cite:`Gehrmann2012yg`. Finally, even higher
multiplicities, calculated at leading order, are merged on top of
that.  A few things to note are detailed below the previous
:ref:`LHC_WJets` example and apply also to this example.


.. literalinclude:: /../../Examples/V_plus_Jets/LHC_ZJets/Sherpa.yaml
   :language: yaml
