.. _QCD continuum:

QCD continuum
=============

Example setup for QCD continuum production at the Belle/KEK collider.
Please note, it does not include any hadronic resonance.

.. literalinclude:: /../../Examples/BFactory/Belle_QCDBackground/Sherpa.yaml
   :language: yaml

Things to notice:

* Asymmetric beam energies, photon ISR is switched on per default.

* Full mass effects of c and b quarks computed.

* Strong coupling constant value set to 0.1188 and two loop (NLO)
  running.

.. _Signal process:

Signal process
==============

Example setup for B-hadron pair production on the Y(4S) pole.

.. literalinclude:: /../../Examples/BFactory/Belle_Signal/Sherpa.yaml
   :language: yaml

Things to notice:

* Same setup as :ref:`QCD continuum`, except for process specification.

* Production of both B0 and B+ pairs, in due proportion.

.. _Single hadron decay chains:

Single hadron decay chains
==========================

This setup is not a collider setup, but a simulation of a hadronic
decay chain.

.. literalinclude:: /../../Examples/BFactory/BDecay/Sherpa.yaml
   :language: yaml

Things to notice:

* ``EVENT_TYPE`` is set to ``HadronDecay``.

* ``DECAYER`` specifies the hadron flavour initialising the decay
  chain.

* A place holder process is declared such that the Sherpa framework
  can be initialised. That process will not be used.
