.. _Integration:

***********
Integration
***********

The following parameters are used to steer the integration:

.. contents::
   :local:

.. _int_INTEGRATION_ERROR:

INTEGRATION_ERROR
=================

.. index:: INTEGRATION_ERROR

Specifies the relative integration error target.

.. _int_INTEGRATOR:

INTEGRATOR
==========

.. index:: INTEGRATOR

Specifies the integrator. The possible integrator types depend on the
matrix element generator. In general users should rely on the default
value and otherwise seek the help of the authors, see :ref:`Authors`.
Within AMEGIC++ the options ``AMEGIC: {INTEGRATOR: <type>,
RS_INTEGRATOR: <type>}`` can be used to steer the behaviour of the
default integrator.

* :option:`4`: building up the channels is achieved through respecting
  the peak structure given by the propagators. The algorithm works
  recursively starting from the initial state.

* :option:`5`: this is an extension of option 4. In the case of
  competing peaks (e.g. a Higgs boson decaying into W+W-, which
  further decay), additional channels are produced to account for all
  kinematical configurations where one of the propagating particles is
  produced on its mass shell.

* :option:`6`: in contrast to option 4 the algorithm now starts from
  the final state. The extra channels described in option 5 are
  produced as well.  This is the default integrator if both beams are
  hadronic.

* :option:`7`: Same as option :option:`4` but with tweaked
  exponents. Optimised for the integration of real-subtracted
  matrix-elements. This is the default integrator when at least one of
  the beams is not hadronic.

In addition, a few ME-generator independent integrators have been
implemented for specific processes:

* :option:`Rambo`: RAMBO :cite:`Kleiss1985gy`. Generates isotropic
  final states.

* :option:`VHAAG`: Vegas-improved HAAG integrator
  :cite:`vanHameren2002tc`.

* :option:`VHAAG_res`: is an integrator for a final state of a weak
  boson, decaying into two particles plus two or more jets based on
  HAAG :cite:`vanHameren2002tc`.  This integrator can be further
  configured using ``VHAAG`` sub-settings, i.e.  ``VHAAG:
  {<sub-setting>: <value>``}. The following sub-settings are
  available. ``RES_KF`` specifies the kf-code of the weak boson, the
  default is W (``24``).  ``RES_D1`` and ``RES_D2`` define the
  positions of the Boson decay products within the internal naming
  scheme, where ``2`` is the position of the first outgoing
  particle. The defaults are ``2`` and ``3``, respectively, which is
  the correct choice for all processes where the decay products are
  the only not strongly interacting final state particles.

.. _VEGAS_MODE:

VEGAS_MODE
==========

.. index:: VEGAS_MODE

Specifies the mode of the Vegas adaptive integration. :option:`0` disables
Vegas, :option:`2` enables it (default).

.. _FINISH_OPTIMIZATION:

FINISH_OPTIMIZATION
===================

.. index:: FINISH_OPTIMIZATION

Specifies whether the full Vegas optimization is to be carried out.
The two possible options are :option:`true` (default) and :option:`false`.

.. _PSI:

PSI
===

.. index:: PSI

The sub-settings for the phase space integrator can be customised as follows:

.. code-block:: yaml

   PSI:
     <sub-setting>: <value>
     # more PSI settings ...

The following sub-settings exist:

``ITMIN``
  The minimum number of points used for every optimisation cycle. Please note
  that it might be increased automatically for complicated processes.

``ITMAX``
  The maximum number of points used for every optimisation cycle. Please note
  that for complicated processes the number given might be insufficient for
  a meaningful optimisation.

``NPOWER``
  The power of two, by which the number of points increases with every step of the
  optimisation.

``NOPT``
  The number of optimization cycles.

``MAXOPT``
  The minimal number of integration cycles after the optimization is done.

``STOPOPT``
  The maximal number of additional cycles in the integration performed to reach the
  integration error goal.

``ITMIN_BY_NODE``
  Same as ``ITMIN``, but specified per node to allow tuning of
  integration performance in large-scale MPI runs.

``ITMAX_BY_NODE``
  Same as ``ITMAX``, but specified per node to allow tuning of
  integration performance in large-scale MPI runs.
