.. _s-channel single-top production:

s-channel single-top production
===============================

.. literalinclude:: /examples/SingleTop_Channels/Sherpa.tj-s_channel.yaml
   :language: yaml

Things to notice:

* By excluding the bottom quark from the initial-state at Born level
  using :option:`PARTICLE_CONTAINERS`, and by setting
  :option:`Max_N_TChannels: 0`, only s-channel diagrams are used for
  the calculation

* See :ref:`t-channel single-top production` for more comments
.. _t-channel single-top production:

t-channel single-top production
===============================

.. literalinclude:: /examples/SingleTop_Channels/Sherpa.tj-t_channel.yaml
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

.. literalinclude:: /examples/SingleTop_Channels/Sherpa.tj-t_channel-nf4.yaml
   :language: yaml

Things to notice:

* We use an Nf4 PDF and use its definition of the bottom mass

* See :ref:`t-channel single-top production` for more comments
.. _tW-channel single-top production:

tW-channel single-top production
================================

.. literalinclude:: /examples/SingleTop_Channels/Sherpa.tW.yaml
   :language: yaml

Things to notice:

* By setting :option:`No_Decay: -6`, the doubly-resonant TTbar
  diagrams are removed. Only the singly-resonant diagrams remain as
  required by the definition of the channel.

* See :ref:`t-channel single-top production` for more comments
