.. _Top quark pair production:

Top quark pair production
=========================

.. literalinclude:: /../../Examples/Tops_plus_Jets/LHC_Tops/Sherpa.yaml
   :language: yaml

Things to notice:

* We use OpenLoops to compute the virtual corrections
  :cite:`Cascioli2011va`.

* We match matrix elements and parton showers using the MC@@NLO
  technique for massive particles, as described in
  :cite:`Hoeche2013mua`.

* A non-default METS core scale setter is used, cf. :ref:`METS scale
  setting with multiparton core processes`

* We enable top decays through the internal decay module using
  :option:`HARD_DECAYS:Enabled: true`
