.. _LHC_HInt:

H production in gluon fusion with interference effects
======================================================

.. index:: HIGGS_INTERFERENCE_MODE
.. index:: HIGGS_INTERFERENCE_ONLY
.. index:: HIGGS_INTERFERENCE_SPIN

This is a setup for inclusive Higgs production through gluon fusion at
hadron colliders. The inclusive process is calculated at
next-to-leading order accuracy, including all interference effects
between Higgs-boson production and the SM gg->yy background.  The
corresponding matrix elements are taken from :cite:`Bern2002jx` and
:cite:`Dixon2013haa`.

.. literalinclude:: /../../Examples/H_in_GluonFusion/LHC_HInterference/Sherpa.yaml
   :language: yaml

Things to notice:

* This calculation is at fixed-order NLO.

* All scales, i.e. the factorisation, renormalisation and resummation
  scales are set to the invariant mass of the di-photon pair.

* Dedicated phase space generators are used by setting
  :option:`Integrator: PS2` and :option:`RS_Integrator: PS3`,
  cf. :ref:`Integrator`.

To compute the interference contribution only, as was done in
:cite:`Dixon2013haa`, one can set :option:`HIGGS_INTERFERENCE_ONLY:
1`.  By default, all partonic processes are included in this
simulation, however, it is sensible to disable quark initial states at
the leading order. This is achieved by setting
:option:`HIGGS_INTERFERENCE_MODE: 3`.

One can also simulate the production of a spin-2 massive graviton in
Sherpa using the same input card by setting
:option:`HIGGS_INTERFERENCE_SPIN: 2`. Only the massive graviton case
is implemented, specifically the scenario where k_q=k_g. NLO
corrections are approximated, as the gg->X->yy and qq->X->yy loop
amplitudes have not been computed so far.
