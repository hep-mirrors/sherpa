.. _HERA:

DIS at HERA
===========

This is an example of a setup for hadronic final states in
deep-inelastic lepton-nucleon scattering at a centre-of-mass energy of
300 GeV.  Corresponding measurements were carried out by the H1 and
ZEUS collaborations at the HERA collider at DESY Hamburg.

.. literalinclude:: /../../Examples/Jets_in_DIS/HERA/Sherpa.yaml
   :language: yaml

Things to notice:

* the beams are asymmetric with the positrons at an energy of 27.5
  GeV, while the protons carry 820 GeV of energy.

* the multi-jet merging cut is set dynamically for each event,
  depending on the photon virtuality, see :cite:`Carli2009cg`.

* there is a selector cut on the photon virtuality. This cut
  implements the experimental requirements for identifying the
  deep-inelastic scattering process.
