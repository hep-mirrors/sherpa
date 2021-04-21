.. _LHC_Wbb:

W+bb production
===============

This example is currently broken. Please contact the :ref:`Authors` for more
information.

..
  .. literalinclude:: /../../Examples/V_plus_Bs/LHC_Wbb/Sherpa.yaml
      :language: yaml
..
   Things to notice:

   * The matrix elements are interfaced from :cite:`FebresCordero2006sj`.
     The shared library necessary for running this setup is built using
     `scons <http://www.scons.org/>`_.

   * The W-boson is stable in the hard matrix elements.
     It is decayed using the internal decay module, indicated by
     the settings :option:`HARD_DECAYS:Enabled: true` and
     :option:`PARTICLE_DATA:24:Stable: 0`.

   * fjcore from `FastJet <http://www.fastjet.fr/>`_ is used to regularize the hard cross section.
     We require two b-jets, indicated by the :option:`2`
     at the end of the :option:`FastjetFinder` options.

   * Four-flavour PDFs are used to comply with the calculational setup.


.. _LHC_Zbb:

Zbb production
==============
.. literalinclude:: /../../Examples/V_plus_Bs/LHC_Zbb/Sherpa.yaml
      :language: yaml

Things to notice:

* The matrix elements are interfaced via the
  Binoth Les Houches interface proposal :cite:`Binoth2010xt`,
  :cite:`Alioli2013nda`, :ref:`External one-loop ME`.

* The Z-boson is stable in the hard matrix elements.
  It is decayed using the internal decay module, indicated by
  the settings :option:`HARD_DECAYS:Enabled: true` and
  :option:`PARTICLE_DATA:23:Stable: 0`.

* fjcore from `FastJet <http://www.fastjet.fr/>`_ is used to regularize the hard cross section.
  We require two b-jets, indicated by :option:`Nb: 2`
  at the end of the :option:`FastjetFinder` options.

* Four-flavour PDF are used to comply with the calculational setup.
