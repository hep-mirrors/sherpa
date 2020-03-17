.. _MINLO:

MINLO
=====

.. index:: MINLO

The following configuration file shows how to implement the MINLO
procedure from :cite:`Hamilton2012np`. A few things to note are
detailed below.  MINLO can also be applied when reading NTuples, see
:ref:`NTuple production`.  In this case, the scale and K factor must
be defined, see :ref:`SCALES` and :ref:`KFACTOR`.

.. literalinclude:: /../../Examples/FixedOrder_NLO/MINLO/Sherpa.yaml
   :language: yaml

Things to notice:

* The R parameter of the flavour-based kT clustering algorithm can be
  changed using ``MINLO:DELTA_R``.

* Setting ``MINLO: {SUDAKOV_MODE: 0``} defines whether to include
  power corrections stemming from the finite parts in the integral
  over branching probabilities.  It defaults to 1.

* The parameter ``MINLO:SUDAKOV_PRECISION`` defines the precision
  target for integration of the Sudakov exponent. It defaults to
  ``1e-4``.
