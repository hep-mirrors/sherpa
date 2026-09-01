.. _Polarized same-sign W-pair production in association with two jets:

Polarized same-sign :math:`\mathrm{W}^+` boson pair production in association with two jets at LO+PS
====================================================================================================
This is an example for the simulation of polarized cross sections for pure electroweak same-sign :math:`\mathrm{W}^+` boson pair production in
association with two jets at LO+PS.

.. literalinclude:: ./examples/VV_plus_Jets/LHC_ssWW_polarization/Sherpa.yaml
   :language: yaml

Things to notice:

* The ``COMIX_DEFAULT_GAUGE`` is set to a special value to get the polarization vectors mentioned in section
  :ref:`Simulation of polarized cross sections for intermediate particles` .
* The width of the :math:`\mathrm{W}^\pm` boson is set to zero to retain gauge invariance since it is considered as stable during the hard
  scattering process (``PROCESSES``). To preserve SU(2) Ward Identities also the width of the Z boson need to be set to
  zero then. 
* The process is divided into the production of the vector bosons (``PROCESSES``) and their decays (``HARD_DECAYS``), all matrix elements are 
  calculated with on-shell vector bosons (narrow-width approximation). ``Mass_Smearing`` and ``SPIN_CORRELATIONS`` are enabled per default.

