.. _t-channel single-top production:

t-channel single-top production
===============================

.. literalinclude:: /../../Examples/SingleTop_Channels/Sherpa.tj-t_channel.yaml
   :language: yaml

Things to notice:

* We use OpenLoops to compute the virtual corrections
  :cite:`Cascioli2011va`.

* We match matrix elements and parton showers using the MC\@NLO
  technique for massive particles, as described in
  :cite:`Hoeche2013mua`.

* A non-default METS core scale setter is used, cf. :ref:`METS scale
  setting with multiparton core processes`

* We enable top and W decays through the internal decay module using
  :option:`HARD_DECAYS:Enabled: true`.  The W is restricted to its
  leptonic decay channels.

* By setting :option:`Min_N_TChannels: 1`, only t-channel diagrams are
  used for the calculation

.. _t-channel single-top production with N_f=4:

t-channel single-top production with N_f=4
==========================================

.. literalinclude:: /../../Examples/SingleTop_Channels/Sherpa.tj-t_channel-nf4.yaml
   :language: yaml

Things to notice:

* We use an Nf4 PDF and use its definition of the bottom mass

* See :ref:`t-channel single-top production` for more comments
