.. _HEFT:

Effective Higgs Couplings
-------------------------

.. index:: FINITE_TOP_MASS
.. index:: FINITE_W_MASS
.. index:: DEACTIVATE_PPH
.. index:: DEACTIVATE_GGH

The HEFT describes the effective coupling of gluons and photons to
Higgs bosons via a top-quark loop, and a W-boson loop in case of
photons. This supplement to the Standard Model can be invoked by
configuring ``MODEL: HEFT``.

The effective coupling of gluons to the Higgs boson, g_ggH, can be
calculated either for a finite top-quark mass or in the limit of an
infinitely heavy top using the switch ``FINITE_TOP_MASS: true`` or
``FINITE_TOP_MASS: false``, respectively. Similarily, the
photon-photon-Higgs coupling, g_ppH, can be calculated both for finite
top and/or W masses or in the infinite mass limit using the switches
``FINITE_TOP_MASS`` and ``FINITE_W_MASS``. The default choice for both
is the infinite mass limit in either case. Note that these switches
affect only the calculation of the value of the effective coupling
constants. Please refer to the example setup :ref:`LHC_HJets_FiniteMt`
for information on how to include finite top quark mass effects on a
differential level.

Either one of these couplings can be switched off using the
``DEACTIVATE_GGH: true`` and ``DEACTIVATE_PPH: true`` switches.  Both
default to ``false``.
