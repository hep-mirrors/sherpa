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

The one-loop EW Sudakov approximation, dubbed EWSud, has been developed for general processes
in :cite:`Denner2000jv,Denner2001gw`. A corresponding automated implementation in the
Sherpa framework, applicable to all common event generation modes of Sherpa,
including multijet-merged calculations, has been presented
in :cite:`Bothmann2020sxm` and :cite:`Bothmann2021led`.

Another available approximation, dubbed EWVirt, was devised in :cite:`Kallweit2015dum`.
It comprises exact renormalised NLO EW virtual corrections and integrated
approximate real-emission subtraction terms, thereby neglecting in particular hard
real-emission contributions. However, both methods qualify for a rather straightforward
inclusion of the dominant EW corrections in state-of-the-art matrix-element plus
parton-shower simulations.

In the following we will discuss how to enable the calculation of the EWSud
and EWVirt corrections, and what options are available to steer their
evaluation, beginning with EWVirt.

.. contents::
   :local:

.. _EWVirt:

EWVirt
======

One option to enable EWVirt corrections is to use ``KFACTOR: EWVirt``.  Note
that this only works for LO calculations (both with and without the shower,
including MEPSatLO).  The EW virtual matrix element must be made available (for
all process multiplicities) using a suitable :ref:`Loop_Generator`.  The
EWVirt correction will then be directly applied to the nominal event weight.

The second option, which is only available for MEPSatNLO, applies the EWVirt
correction (and optionally subleading LO corrections) to all QCD NLO
multiplicities. For this to work, one must use the the following syntax:

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
output. The alternative event weight names
are either ``ASSOCIATED_CONTRIBUTIONS.<contrib>``,
``ASSOCIATED_CONTRIBUTIONS.MULTI<contrib>``,
or ``ASSOCIATED_CONTRIBUTIONS.EXP<contrib>``
for additive, multiplicative, and exponentiated combinations, correspondingly.
See :ref:`On-the-fly event weight variations` for more information
on variation weights and the variation weight naming scheme.

.. _EWSud:

EWSud
=====

The EWSud module must be enabled during configuration of Sherpa using the
``-DSHERPA_ENABLE_EWSUD=ON`` switch.

Similar to EWVirt, also with the EWSud corrections there is the option to use
it via ``KFACTOR: EWSud``, which will apply the corrections directly to the
nominal event weight, or as on-the-fly variations adding the following entry to
the list of variations (also cf. :ref:`On-the-fly event weight variations`):

.. code-block:: yaml

   VARIATIONS:
   - EWSud

Using the latter, corrections are provided as alternative event weights.
The most useful entries of the event weight list are accessed using the keys
`EWSud.KFactor` and `EWSud.KFactorExp`.
The first is the nominal event weight corrected by the
NLL EWSud corrections, while the latter first exponentiates the corrections
prior to applying it to the nominal event weight, thus giving a resummed NLL
result.

In order for the ``EWSud`` corrections to make sense, Goldstone bosons need
to be made available. This is achieved by ensuring that the following is set

.. code-block:: yaml

   MODEL: SMGold

Additionally, a coupling order must be set to correctly initialize the
couplings for this model, see :ref:`Processes` for more details

.. code-block:: yaml

   PROCESSES:
     ...
     Order{QCD:xx, EW:yy, SMGold: 0}
     ...

The following configuration snippet shows the options steering the EWSud
calculation, along with their default values:

.. code-block:: yaml

   EWSUD:
     THRESHOLD: 1.0
     INCLUDE_SUBLEADING: true
     CLUSTERING_THRESHOLD: 10.0

.. index:: THRESHOLD

* :option:`THRESHOLD` . Strictly speaking the EWSudakov corrections are only
  valid in the high-energy limit, that is where all possible invariant masses,
  formed by pairing external particles, are much larger
  than the W mass. In practice, we need to define how much is much larger.
  The :option:`THRESHOLD` option, gives the minimal invariant mass (in units of
  :math:m_W) that each
  pairing of external particles can have to respect the high energy limit, and
  below which no EWSudakov correction is computed. To clarify, a large
  threshold, say for example 10 (10 times the W mass), would result in little to
  no corrections at all, except for regions of phase-space truly in the
  high-energy limit. This result is thus only expected to match exact EW
  corrections only when all invariants are larger than this
  threshold. Conversely a lower value, say 1, would apply the correction more
  uniformly at the price of violating the theoretically sound region where these
  corrections are derived, but is seen to better reproduce the effect of exact
  EW corrections across kinematical distributions.

.. index:: INCLUDE_SUBLEADING

* :option:`INCLUDE_SUBLEADING` determines whether a formally subleading term
  proportional to :math:`\log^2(r_{kl} / \hat s)` is included,
  where :math:`\hat s` is the Mandelstam variable for the partonic process,
  see :cite:`Bothmann2021led`. Note that depending on the value of
  :option:`THRESHOLD` these may become numerically significant. For lower threshold
  values, it is recommended to leave this option `true`, as default.

.. index:: CLUSTERING_THRESHOLD

* :option:`CLUSTERING_THRESHOLD` determines the number of vector boson decay widths,
  for which a given lepton pair with the right quantum numbers is still allowed
  to be clustered prior to the calculation of the EWSud correction.
  For reasoning, see again :cite:`Bothmann2021led`.

We next list all possible technical parameters under the scope of `EWSUD`. They
are mostly meant for internal or consistency checks and are advisable only to
expert users.

.. index:: RS

* :option:`RS` boolean flag to determine whether or not to apply the EWSudakov
  corrections to `RS` type events, defaults to `true`.

.. index:: CHECK

* :option:`CHECK` boolean flag to enable/disable internal checks on the
  logarithmic coefficients for various simple processes. Defaults to false and
  prevents normal running when set to true, in that it terminates the run after
  having checked the coefficients.

.. index:: CHECK_KFACTOR

* :option:`CHECK_KFACTOR` Same as `CHECK` but at the level of `KFACTOR`.

.. index:: CHECK_LOG_FILE

* :option:`CHECK_LOG_FILE` Specify a filename in which to store the result of
  `CHECK`, defaults to a null string.

.. index:: CHECKINVARIANTRATIOS

* :option:`CHECKINVARIANTRATIOS` boolean flag used to enforce a stricter
  definition of High Energy Limit, defaults to false.

.. index:: COEFF_REMODED_LIST

* :option:`COEFF_REMOVED_LIST` list of logarithmic coefficients that can be
  ignored, defaults to empty, meaning that all coefficients are included. The
  available options are: `LSC`, `Z`, `SSC`, `C`, `Yuk`, `PR` and `I`. See
  :cite:`Bothmann2020sxm` for further details.

.. index:: C_COEFF_IGNORES_VECTOR_BOSONS

* :option:`C_COEFF_IGNORES_VECTOR_BOSONS` boolean flag to control whether or not
  Vector Boson contributions should be included in the calculation of the `C`
  coefficient. Defaults to false, and can be used to check the `PR` logarithms,
  given that for some procs the contributions to `C` from vector bosons and the
  `PR` coefficients cancel.

.. index:: HIGH_ENERGY_SCHEME

* :option:`HIGH_ENERGY_SCHEME` different implementations of the High Energy
  limit conditions. At the moment only `Default` is fully implemented, all other
  available options imply that no check is enforced on the configurations, and a
  contribution is calculated independently on whether or we are in the high
  energy limit.

.. index:: PRINT_GRAPHS

* :option:`PRINT_GRAPHS` sets the name of the directory where to save graphs
  associated to processes generated by the `EWSudakov` calculation. Same as
  :ref:`Print_Graphs`.

**NOTE**
that at the moment EW Sudakov corrections do not work for processes that feature a four-vector boson vertex, such as a four-gluon vertex.
