.. _LHC_Jets:

Jet production
==============

.. contents::
   :local:

To change any of the following LHC examples to production at the Tevatron
simply change the beam settings to

.. code-block:: yaml

   BEAMS: [2212, -2212]
   BEAM_ENERGIES: 980

.. _LHC_Jets_MCatNLO:

MC\@NLO setup for dijet and inclusive jet production
----------------------------------------------------

This is an example setup for dijet and inclusive jet production at hadron
colliders at next-to-leading order precision matched to the parton shower
using the MC\@NLO prescription detailed in :cite:`Hoeche2011fd` and
:cite:`Hoeche2012fm`. A few things to note are detailed below the example.

.. literalinclude:: /../../Examples/Jets_at_HadronColliders/LHC_Jets_MCatNLO/Sherpa.yaml
   :language: yaml

Things to notice:

* Asymmetric cuts are implemented (relevant to the RS-piece of an
  MC\@NLO calculation) by requiring at least two jets with pT > 10
  GeV, one of which has to have pT > 20 GeV.

* Both the factorisation and renormalisation scales are set to the
  above defined scale factors times a quarter of the scalar sum of the
  transverse momenta of all anti-kt jets (R = 0.4, pT > 20 GeV) found
  on the ME-level before any parton shower emission. See :ref:`SCALES`
  for details on scale setters.

* The resummation scale, which sets the maximum scale of the
  additional emission to be resummed by the parton shower, is set to
  the above defined resummation scale factor times half of the
  transverse momentum of the softer of the two jets present at Born
  level.

* The external generator OpenLoops provides the one-loop matrix
  elements.

* The ``NLO_QCD_Mode`` is set to ``MC@NLO``.

.. _LHC_Jets_MEPS:

MEPS setup for jet production
-----------------------------

.. literalinclude:: /../../Examples/Jets_at_HadronColliders/LHC_Jets_MEPS/Sherpa.yaml
   :language: yaml

Things to notice:

* ``Order`` is set to :option:`{QCD: Any, EW: 0`}. This
  ensures that all final state jets are produced via
  the strong interaction.

* An ``NJetFinder`` selector is used to
  set a resolution criterion for the two jets of the core process.
  This is necessary because the "CKKW" tag does not apply any cuts to the core
  process, but only to the extra-jet matrix elements, see
  :ref:`Multijet merged event generation with Sherpa`.
