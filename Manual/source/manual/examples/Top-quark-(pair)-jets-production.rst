.. _Top quark pair production:

Top quark pair production
=========================

.. literalinclude:: /../../Examples/Tops_plus_Jets/LHC_Tops/Sherpa.yaml
   :language: yaml

Things to notice:

* We use OpenLoops to compute the virtual corrections
  :cite:`Cascioli2011va`.

* We match matrix elements and parton showers using the MC\@NLO
  technique for massive particles, as described in
  :cite:`Hoeche2013mua`.

* A non-default METS core scale setter is used, cf. :ref:`METS scale
  setting with multiparton core processes`

* We enable top decays through the internal decay module using
  :option:`HARD_DECAYS:Enabled: true`

* We calculate on-the-fly a 7-point scale variation,
  cf. :ref:`On-the-fly event weight variations`.

.. _Top quark pair production including approximate EW corrections:

Top quark pair production incuding approximate EW corrections
=============================================================

.. literalinclude:: /../../Examples/Tops_plus_Jets/LHC_Tops/Sherpa.EWapprox.yaml
   :language: yaml

Things to notice:

* In addition to the setup in :ref:`Top quark pair production` we add
  approximate EW corrections, cf. :cite:`Gutschow2018tuk`.

* Please note: this setup only works with OpenLoops v.2 or later.

* The approximate EW corrections are added as additional variations on the
  event weight.
