.. _Scale variations:

#############################
A posteriori scale variations
#############################

There are several ways to compute the effects of changing the scales
and PDFs of any event produced by Sherpa. They can computed
explicitly, cf. :ref:`Explicit scale variations`, on-the-fly, cf.
:ref:`On-the-fly event weight variations` (restricted to multiplicative
factors), or reconstructed a posteriori. The latter method needs
plenty of additional information in the event record and is (depending
on the actual calculation) available in two formats:

.. contents::
   :local:

.. _A posteriori scale and PDF variations using the HepMC GenEvent Output:

*********************************************************************
A posteriori scale and PDF variations using the HepMC GenEvent Output
*********************************************************************

Events generated in a LO, LOPS, NLO, NLOPS, MEPS\@LO, MEPS\@NLO or
MENLOPS calculation can be written out in the HepMC format including
all infomation to carry out arbitrary scale variations a
posteriori. For this feature HepMC of at least version 2.06 is
necessary and both ``HEPMC_USE_NAMED_WEIGHTS: true`` and
``HEPMC_EXTENDED_WEIGHTS: true`` have to enabled. Detailed
instructions on how to use this information to construct the new event
weight can be found here
`<https://sherpa.hepforge.org/doc/ScaleVariations-Sherpa-2.2.0.pdf>`_.

.. _A posteriori scale and PDF variations using the ROOT NTuple Output:

******************************************************************
A posteriori scale and PDF variations using the ROOT NTuple Output
******************************************************************

.. index:: USR_WGT_MODE

Events generated at fixed-order LO and NLO can be stored in ROOT
NTuples that allow arbitrary a posteriori scale and PDF variations,
see :ref:`Event output formats`. An example for writing and reading in
such ROOT NTuples can be found here: :ref:`NTuple production`.  The
internal ROOT Tree has the following Branches:

``id``
  Event ID to identify correlated real sub-events.

``nparticle``
  Number of outgoing partons.

``E/px/py/pz``
  Momentum components of the partons.

``kf``
  Parton PDG code.

``weight``
  Event weight, if sub-event is treated independently.

``weight2``
  Event weight, if correlated sub-events are treated as single event.

``me_wgt``
  ME weight (w/o PDF), corresponds to 'weight'.

``me_wgt2``
  ME weight (w/o PDF), corresponds to 'weight2'.

``id1``
  PDG code of incoming parton 1.

``id2``
  PDG code of incoming parton 2.

``fac_scale``
  Factorisation scale.

``ren_scale``
  Renormalisation scale.

``x1``
  Bjorken-x of incoming parton 1.

``x2``
  Bjorken-x of incoming parton 2.

``x1p``
  x' for I-piece of incoming parton 1.

``x2p``
  x' for I-piece of incoming parton 2.

``nuwgt``
  Number of additional ME weights for loops and integrated subtraction terms.

``usr_wgt[nuwgt]``
  Additional ME weights for loops and integrated subtraction terms.

Computing (differential) cross sections of real correction events with statistical errors
=========================================================================================

Real correction events and their counter-events from subtraction terms are
highly correlated and exhibit large cancellations. Although a treatment of
sub-events as independent events leads to the correct cross section the
statistical error would be greatly overestimated. In order to get a realistic
statistical error sub-events belonging to the same event must be combined
before added to the total cross section or a histogram bin of a differential
cross section. Since in general each sub-event comes with it's own set of four
momenta the following treatment becomes necessary:

#. An event here refers to a full real correction event that may
   contain several sub-events. All entries with the same id belong to
   the same event.  Step 2 has to be repeated for each event.

#. Each sub-event must be checked separately whether it passes
   possible phase space cuts. Then for each observable add up
   ``weight2`` of all sub-events that go into the same histogram
   bin. These sums :math:`x_{id}` are the quantities to enter the actual
   histogram.

#. To compute statistical errors each bin must store the sum over all
   :math:`x_{id}` and the sum over all :math:`x_{id}^2`. The cross section
   in the bin is given by :math:`\langle x\rangle = \frac{1}{N} \cdot
   \sum x_{id}`, where :math:`N` is the number of events (not
   sub-events). The :math:`1-\sigma` statistical error for the bin is
   :math:`\sqrt{ (\langle x^2\rangle-\langle x\rangle^2)/(N-1) }`

Note: The main difference between ``weight`` and ``weight2`` is that they
refer to a different counting of events. While ``weight`` corresponds to
each event entry (sub-event) counted separately, ``weight2`` counts events
as defined in step 1 of the above procedure. For NLO pieces other than the real
correction ``weight`` and ``weight2`` are identical.

Computation of cross sections with new PDF's
============================================

Born and real pieces
--------------------

Notation:

.. code-block:: text

   f_a(x_a) = PDF 1 applied on parton a, F_b(x_b) = PDF 2 applied on
   parton b.

The total cross section weight is given by:

.. code-block:: text

   weight = me_wgt f_a(x_a)F_b(x_b)

Loop piece and integrated subtraction terms
-------------------------------------------

The weights here have an explicit dependence on the renormalization
and factorization scales.

To take care of the renormalization scale dependence (other than via
``alpha_S``) the weight ``w_0`` is defined as


.. code-block:: text

   w_0 = me_wgt + usr_wgts[0] log((\mu_R^new)^2/(\mu_R^old)^2) +
   usr_wgts[1] 1/2 [log((\mu_R^new)^2/(\mu_R^old)^2)]^2

To address the factorization scale dependence the weights ``w_1,...,w_8``
are given by

.. code-block:: text

   w_i = usr_wgts[i+1] + usr_wgts[i+9] log((\mu_F^new)^2/(\mu_F^old)^2)

The full cross section weight can be calculated as

.. code-block:: text

   weight = w_0 f_a(x_a)F_b(x_b)
             + (f_a^1 w_1 + f_a^2 w_2 + f_a^3 w_3 + f_a^4 w_4) F_b(x_b)
             + (F_b^1 w_5 + F_b^2 w_6 + F_b^3 w_7 + F_b^4 w_8) f_a(x_a)

where

.. code-block:: text

   f_a^1 = f_a(x_a) (a=quark), \sum_q f_q(x_a) (a=gluon),
   f_a^2 = f_a(x_a/x'_a)/x'_a (a=quark), \sum_q f_q(x_a/x'_a)x'_a (a=gluon),
   f_a^3 = f_g(x_a),
   f_a^4 = f_g(x_a/x'_a)/x'_a

The scale dependence coefficients ``usr_wgts[0]`` and ``usr_wgts[1]``
are normally obtained from the finite part of the virtual correction
by removing renormalization terms and universal terms from dipole
subtraction.  This may be undesirable, especially when the loop
provider splits up the calculation of the virtual correction into
several pieces, like leading and sub-leading color. In this case the
loop provider should control the scale dependence coefficients, which
can be enforced with option :option:`USR_WGT_MODE: false`.

.. warning::

   The loop provider must support this option or the scale dependence
   coefficients will be invalid!
