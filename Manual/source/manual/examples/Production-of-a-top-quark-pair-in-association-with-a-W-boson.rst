.. _Tops plus W in MC@@NLO:

Production of a top quark pair in association with a W-boson
============================================================

.. literalinclude:: /../../Examples/Tops_plus_V/LHC_TTW/Sherpa.yaml
   :language: yaml

Things to notice:

* Hard decays are enabled through :option:`HARD_DECAYS:Enabled: true`.

* Top quarks and W bosons are final states in the hard matrix elements, so
  their widths are set to zero using :option:`Width: 0` in their
  :option:`PARTICLE_DATA` settings.

* Certain decay channels are disabled using :option:`Status: 0` in the
  :option:`Channels` sub-settings of the :option:`HARD_DECAYS` setting.


