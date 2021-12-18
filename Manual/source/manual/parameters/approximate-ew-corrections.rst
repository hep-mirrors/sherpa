.. _Approximate Electroweak Corrections:

***********************************
Approximate Electroweak Corrections
***********************************

As an alternative to the complete set of NLO EW corrections, methods restricted
to the leading effects due to EW loops are available in Sherpa. In particular
at energy scales :math:`Q` large compared to the masses of the EW gauge bosons,
contributions from virtual W- and Z-boson exchange and
corresponding collinear real emissions dominate. The leading contributions are
Sudakov-type logarithms of the form :cite:`Sudakov1954sw,Ciafaloni1998xg`.

.. math::

  \frac{\alpha}{4\pi \sin^2\theta_W}\log^2\left(\frac{Q^2}{M^2_W}\right)\quad\text{and}\quad
  \frac{\alpha}{4\pi \sin^2\theta_W}\log\left(\frac{Q^2}{M^2_W}\right)\,.

The one-loop EW Sudakov approximation, dubbed EWsud, has been developed for general processes
in :cite:`Denner2000jv,Denner2001gw`. A corresponding automated implementation in the
Sherpa framework, applicable to all common event generation modes of Sherpa,
including multijet-merged calculations, has been presented
in :cite:`Bothmann2020sxm` and :cite:`Bothmann2021led`.

Another available approximation, dubbed EWvirt, was devised in :cite:`Kallweit2015dum`.
It comprises exact renormalised NLO EW virtual corrections and integrated
approximate real-emission subtraction terms, thereby neglecting in particular hard
real-emission contributions. However, both methods qualify for a rather straightforward
inclusion of the dominant EW corrections in state-of-the-art matrix-element plus
parton-shower simulations.

In the following we will discuss how to enable the calculation of thew EWsud
and EWvirt corrections, and what options are available to steer their
evaluation, beginning with EWvirt.

.. contents::
   :local:

.. _EWVirt:

EWvirt
======

One option to enable EWvirt corrections is to use ``KFACTOR: EWvirt``.  Note
that this only works for LO calculations (both with and without the shower,
including MEPSatLO).  The EW virtual matrix element must be made available (for
all process multiplicities) using a suitable :ref:`Loop_Generator`.  The
EWvirt correction will then be directly applied to the nominal event weight.

The second option, which is only available for MEPSatNLO, applies the EWvirt
correction (and optionally subleading LO corrections) to all QCD NLO
multiplities. For this to work, one must use the the following syntax:

.. index:: ASSOCIATED_CONTRIBUTIONS_VARIATIONS

.. code-block:: yaml

   ASSOCIATED_CONTRIBUTIONS_VARIATIONS:
   - [EW]
   - [EW, LO1]
   - [EW, LO1, LO2]
   - [EW, LO1, LO2, LO3]

Each entry of ``ASSOCIATED_CONTRIBUTIONS_VARIATIONS`` defines a variation and
the different associated contributions that should be taken into account for
the corresponding alternative weight.
Note that the respective associated contribution must be listed
in the process setting :ref:`Associated_Contributions`.

The additional event weights can then be written into the event
output.  However, this is currently only supported for
``HepMC_GenEvent`` and ``HepMC_Short`` with versions >=2.06 and
``HEPMC_USE_NAMED_WEIGHTS: true``.  The alternative event weight
names are either ``ASS<contrib>`` or ``MULTIASS<contrib>``,
for additive and multiplicative combinations, correspondingly.


.. _EWSud:

EWsud
=====

The EWsud module must be enabled during configuration of Sherpa using the
``--enable-ewsud`` switch.

Similar to EWvirt, also with the EWsud corrections there is the option to use
it via ``KFACTOR: EWsud``, which will apply the corrections directly to the
nominal event weight, or as on-the-fly variations adding the following entry to
the list of variations (also cf. :ref:`On-the-fly event weight variations`):

.. code-block:: yaml

   VARIATIONS:
   - EWsud

Using the latter, corrections are provided as alternative event weights.
The most useful entries of the event weight list are accessed using the keys
`EWsud` and `EWsud_Exp`. The first is the nominal event weight corrected by the
NLL EWsud corrections, while the latter first exponentiates the corrections
prior to applying it to the nominal event weight, thus giving a resummed NLL
result.

The following configuration snippet shows the options steering the EWsud
calculation, along with their default values:

.. code-block:: yaml

   EWSud:
     THRESHOLD: 5.0
     INCLUDE_SUBLEADING: false

   ME_QED:
     CLUSTERING_THRESHOLD: 10.0

.. index:: EWSUD_THRESHOLD

* :option:`THRESHOLD` gives the minimal invariant mass (in units of the W mass)
  for each external pair of particles :math:`k` and :math:`l`, :math:`r_{kl}`,
  defining the high energy limit. If any of the invariant masses is below this
  value for a given event, then no EWsud correction is calculated.

.. index:: EWSUD_INCLUDE_SUBLEADING

* :option:`INCLUDE_SUBLEADING` determines whether a formally subleading term
  proportional to :math:`\log^2(r_{kl} / \hat s)` is included,
  where :math:`\hat s` is the Mandelstam variable for the partonic process,
  see :cite:`Bothmann2021led`.

.. index:: ME_QED_CLUSTERING_THRESHOLD

* :option:`CLUSTERING_THRESHOLD` determines the number of vector boson decay widths,
  for which a given lepton pair with the right quantum numbers is still allowed
  to be clustered prior to the calculation of the EWsud correction.
  For reasoning, see again :cite:`Bothmann2021led`.
