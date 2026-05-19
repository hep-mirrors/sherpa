.. _LHC_Zbb:

Zbb production
==============
.. literalinclude:: ./examples/V_plus_Bs/LHC_Zbb/Sherpa.yaml
      :language: yaml

Things to notice:

* The Z-boson is stable in the hard matrix elements.
  It is decayed using the internal decay module, indicated by
  the settings :option:`HARD_DECAYS:Enabled: true` and
  :option:`PARTICLE_DATA:23:Stable: 0`.

* fjcore from `FastJet <http://www.fastjet.fr/>`_ is used to regularize the hard cross section.
  We require two b-jets, indicated by :option:`Nb: 2`
  at the end of the :option:`FastjetFinder` options.

* Four-flavour PDF are used to comply with the calculational setup.

