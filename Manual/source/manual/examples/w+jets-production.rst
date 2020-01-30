.. _LHC_WJets:

W+jets production
=================

This is an example setup for inclusive W production at hadron
colliders.  The inclusive process is calculated at next-to-leading
order accuracy matched to the parton shower using the MC\@NLO
prescription detailed in :cite:`Hoeche2011fd`. The next few higher jet
multiplicities, calculated at next-to-leading order as well, are
merged into the inclusive sample using the MEPS\@NLO method -- an
extension of the CKKW method to NLO -- as described in
:cite:`Hoeche2012yf` and :cite:`Gehrmann2012yg`. Finally, even higher
multiplicities, calculated at leading order, are merged on top of
that.  A few more things to note are detailed below the example.


.. literalinclude:: /../../Examples/V_plus_Jets/LHC_WJets/Sherpa.yaml
   :language: yaml

Things to notice:

* The ``Order`` in the process definition in a multi-jet merged setup
  defines the order of the core process (here ``93 93 -> 90 91`` with
  two electroweak couplings). The additional strong couplings for
  multi-jet production are implicitly understood.

* The settings necessary for NLO accuracy are restricted to the
  2->2,3,4 processes using the ``2->2-4`` key below the ``# set up
  MC@NLO ...`` comment.  The example can be converted into a simple
  MENLOPS setup by using ``2->2`` instead, or into an MEPS setup by
  removing these lines altogether.  Thus one can study the effect of
  incorporating higher-order matrix elements.

* The number of additional LO jets can be varied through changing the
  integer within the curly braces in the ``Process`` definition, which
  gives the maximum number of additional partons in the matrix
  elements.

* ``OpenLoops`` is used here as the provider of the one-loop matrix
  elements for the respective multiplicities.

* Tau leptons are set massive in order to exclude them from the
  massless lepton container (``90``).

* As both Comix and Amegic are specified as matrix element generators
  to be used, Amegic has to be specified to be used for all MC\@NLO
  multiplicities using ``ME_Generator: Amegic``. Additionally, we
  specify ``RS_ME_Generator: Comix`` such that the subtracted
  real-emission bit of the NLO matrix elements is calculated more
  efficiently with Comix instead of Amegic.  This combination is
  currently the only one supported for NLO-matched/merged setups.


The jet criterion used to define the matrix element multiplicity in
the context of multijet merging can be supplied by the user. As an
example the source code file
:file:`./Examples/V_plus_Jets/LHC_WJets/My_JetCriterion.C` provides such
an alternative jet criterion. It can be compiled using ``SCons`` via
executing ``scons`` in that directory (edit the ``SConstruct`` file
accordingly). The newly created library is linked at run time using
the ``SHERPA_LDADD`` flag.  The new jet criterion is then evoked by
``JET_CRITERION``.
