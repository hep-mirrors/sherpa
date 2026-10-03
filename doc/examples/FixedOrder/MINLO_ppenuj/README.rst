.. _MINLO:

MINLO
=====

.. index:: MINLO

The following configuration file shows how to implement the MINLO
procedure from :cite:`Hamilton2012np`. A few things to note are
detailed below.

.. literalinclude:: ./examples/FixedOrder/MINLO_ppenuj/Sherpa.yaml
   :language: yaml

Things to notice:

* The R parameter of the flavour-based kT clustering algorithm can be
  changed using ``MINLO:DELTA_R``.

* ``MINLO: {SUDAKOV_MODE: 3}`` is the default. This is a bitmask:
  bit ``1`` includes finite-bound power corrections in the integrated
  branching probabilities, and bit ``2`` enables the two-loop ``K`` term
  in the Sudakov integrand. The default enables both. Setting
  ``MINLO: {SUDAKOV_MODE: 0}`` disables these corrections, not the
  Sudakov reweighting itself.

* The parameter ``MINLO:SUDAKOV_PRECISION`` defines the precision
  target for integration of the Sudakov exponent. It defaults to
  ``1e-4``.
