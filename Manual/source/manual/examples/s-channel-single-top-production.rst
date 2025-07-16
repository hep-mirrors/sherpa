.. _s-channel single-top production:

s-channel single-top production
===============================

.. literalinclude:: /../../Examples/SingleTop_Channels/Sherpa.tj-s_channel.yaml
   :language: yaml

Things to notice:

* By excluding the bottom quark from the initial-state at Born level
  using :option:`PARTICLE_CONTAINERS`, and by setting
  :option:`Max_N_TChannels: 0`, only s-channel diagrams are used for
  the calculation

* See :ref:`t-channel single-top production` for more comments
