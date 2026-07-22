.. _MPI Parameters:

*********************
Multiple interactions
*********************

The basic MPI model is described
in :cite:`Sjostrand1987su` while Sherpa's implementation details are
discussed in :cite:`Alekhin2005dx`.

The following parameters are used to steer the MPI setup:

.. contents::
   :local:


.. _MI_HANDLER:

MI_HANDLER
==========

.. index:: MI_HANDLER

Specifies the MPI handler. The two possible values
at the moment are :option:`None` and :option:`Amisic`.

.. _Amisic:

AMISIC
======

.. index:: Amisic
.. index:: Amisic:EVENT_TYPE
.. index:: Amisic:PT_0(ref)
.. index:: Amisic:PT_0(IR)
.. index:: Amisic:PT_Min(ref)
.. index:: Amisic:PT_Min(IR)
.. index:: Amisic:Eta
.. index:: Amisic:E(ref)
.. index:: Amisic:PT_Min
.. index:: Amisic:PT_0
.. index:: Amisic:MU_R_SCHEME
.. index:: Amisic:MU_F_SCHEME
.. index:: Amisic:MU_R_FACTOR
.. index:: Amisic:MU_F_FACTOR
.. index:: Amisic:E_Min
.. index:: Amisic:SIGMA_ND_NORM
.. index:: Amisic:nPT_bins
.. index:: Amisic:nS_bins
.. index:: Amisic:nB_bins
.. index:: Amisic:N_MaxScatters
.. index:: Amisic:TRIGGER
.. index:: Amisic:PomeronIntercept
.. index:: Amisic:PomeronSlope
.. index:: Amisic:TriplePomeronCoupling
.. index:: Amisic:ReggeonIntercept
.. index:: Amisic:Diffractive_cres
.. index:: Amisic:Diffractive_Mres
.. index:: Amisic:Diffractive_s1
.. index:: Amisic:ElasticSlope_c0
.. index:: Amisic:ElasticSlope_c1
.. index:: Amisic:TwoPionInterference
.. index:: Amisic:f_omega
.. index:: Amisic:phi_omega
.. index:: Amisic:f_nr
.. index:: Amisic:Lambda_nr
.. index:: Amisic:delta_nr

Amisic can simulate the interaction of three different combinations of incoming particles:
proton--proton, photon--proton and photon--photon collision. The parameters for the simulation of photonic multiple
interactions can be found in :cite:`Schuler:1993wr`. It has several parameters to control the simulation of the
multiple-parton interactions, they are listed below. Each of these parameters has to be set in the
subsetting ``AMISIC``, like so

.. code-block:: yaml

   AMISIC:
     PT_0: 2.5

The usual rules for yaml structure apply, c.f. :ref:`Input structure`.

:option:`EVENT_TYPE`
    Selects the class of event to generate. The default ``Perturbative`` runs the perturbative
    multiple-interaction model; ``Elastic``, ``DiffractiveA``, ``DiffractiveB``, ``DiffractiveAB``
    and ``QuasiElastic`` select the corresponding soft (non-perturbative) processes, and
    ``MinimumBias`` (``Non-Perturbative``) the inclusive soft sample. Defaults to ``Perturbative``.

Perturbative scales
-------------------

:option:`PT_0(ref)`
    Value :math:`p_\text{T,0}^\text{(ref)}` for the calculation of the IR regulator, see formula below. Defaults to ``2.05``.

:option:`PT_0(IR)`
    The absolute minimum (floor) of the IR regulator, see formula below. Defaults to ``0.5``.

:option:`PT_Min(ref)`
    Value :math:`p_\text{T,min}^\text{(ref)}` for the calculation of the IR cutoff, see formula below. Defaults to ``1.10``.

:option:`PT_Min(IR)`
    The absolute minimum (floor) of the IR cutoff, see formula below. Defaults to ``1.0``.

:option:`Eta`
    The energy-scaling exponent :math:`\eta` used to calculate the IR cutoff and regulator, :math:`p_\text{T,min}` and :math:`p_\text{T,0}`.
    Defaults to ``0.026``.

:option:`E(ref)`
    Reference energy to normalise the actual cms energy for the calculation of the IR cutoff and regulator.
    Defaults to ``7000``.

:option:`PT_Min`
    The IR cut-off for the 2->2 scatters. It is calculated as

    .. math::

        p_\text{T,min}^2 = \max\left[ \left(p_\text{T,min}^\text{(IR)}\right)^2,\;
        \left(p_\text{T,min}^\text{(ref)}\right)^2 \left( \frac{s}{(E_\text{cms}^\text{(ref)})^2} \right)^{2\eta} \right]

    but can also be set explicitly.

:option:`PT_0`
    IR regulator :math:`p_\text{T,0}` in the propagator and in the strong coupling. It is calculated as

    .. math::

        p_\text{T,0}^2 = \max\left[ \left(p_\text{T,0}^\text{(IR)}\right)^2,\;
        \left(p_\text{T,0}^\text{(ref)}\right)^2 \left( \frac{s}{(E_\text{cms}^\text{(ref)})^2} \right)^{2\eta} \right]

    but can also be set explicitly.

:option:`MU_R_SCHEME`
    Renormalisation-scale scheme. Defaults to ``PT``. More schemes have yet to be added.

:option:`MU_F_SCHEME`
    Factorisation-scale scheme. Defaults to ``PT``. More schemes have yet to be added.

:option:`MU_R_FACTOR`
    Factor to scale the renormalisation scale :math:`\mu_R`, defaults to ``1.0``.

:option:`MU_F_FACTOR`
    Factor to scale the factorisation scale :math:`\mu_F`, defaults to ``0.5``.

:option:`E_Min`
    Minimum residual beam energy (in GeV) required to continue producing further scatters.
    Defaults to ``0.25``.

Integration controls
--------------------

:option:`SIGMA_ND_NORM`
    Specifies the factor to scale the non-diffractive cross section calculated in the MPI initialisation.
    Defaults to ``0.44``.

:option:`nPT_bins`
    Controls the number of bins for the numerical integration of

    .. math:: \int_{p_T^2}^{s/4} dp_T^2 \frac{d \sigma}{dp_T^2}

    Defaults to ``200``.

:option:`nS_bins`
    Number of points to sample in the center-of-mass energy :math:`\sqrt{s}`. This is only used if the energy is not
    fixed, i.e. in the case of EPA photons. Defaults to ``40``.

:option:`nB_bins`
    Number of bins used to tabulate the impact-parameter dependence of the matter overlap. Defaults to ``20``.

:option:`N_MaxScatters`
    Upper limit on the number of additional scatters generated per event. Defaults to ``10000``.

:option:`TRIGGER`
    Optional list of (up to two) parton flavours (PDG codes) used to bias the minimum-bias generation
    towards events containing the requested final state. Empty by default (no trigger).

Non-diffractive and diffractive cross sections
----------------------------------------------

The total cross-section is calculated with

    .. math::

        \sigma_{tot} = X s^\epsilon + Y s^\eta

    where :math:`s` is the Mandelstam invariant.

:option:`PomeronIntercept`
    The parameter :math:`\epsilon` in the above equation, defaults to ``0.0808``.

:option:`ReggeonIntercept`
    The parameter :math:`\eta` in the above equation, defaults to ``-0.4525``.

The single- and double-diffractive cross-sections in the Regge picture have several free parameters:

:option:`PomeronSlope`
    The Pomeron slope :math:`\alpha^\prime`, default is ``0.25``.

:option:`TriplePomeronCoupling`
    The triple-Pomeron coupling :math:`g_{3\mathbb{P}}` at an input scale of 20 GeV, given in :math:`\text{mb}^{-0.5}`, with default ``0.318``.

:option:`Diffractive_cres`
    Strength of the low-mass resonance enhancement in the diffractive mass distribution, entering as
    :math:`1 + c_\text{res}\, M_\text{res}^2/(M_\text{res}^2+M_X^2)`. Defaults to ``2.0``.

:option:`Diffractive_Mres`
    Mass shift :math:`\Delta M_\text{res}` (in GeV) added to the diffracting hadron's mass to define the
    effective low-mass-resonance mass :math:`M_\text{res} = m_\text{had} + \Delta M_\text{res}` entering
    :math:`1 + c_\text{res}\, M_\text{res}^2/(M_\text{res}^2+M_X^2)`. Defaults to ``2.0``.

:option:`Diffractive_s1`
    Reference scale :math:`s_1` (in :math:`\text{GeV}^2`) at which the triple-Pomeron coupling :math:`g_{3\mathbb{P}}`
    is normalised, entering the single-/double-diffractive normalisation as :math:`s_1^{3\epsilon_\mathbb{P}/2}`
    and :math:`s_1^{\epsilon_\mathbb{P}}` respectively. Defaults to ``400`` (:math:`= (20\,\text{GeV})^2`).

:option:`ElasticSlope_c0`
    Coefficient :math:`c_0` of the energy-dependent (Regge-shrinkage) term :math:`c_0\,(s/s_0)^{\epsilon_\mathbb{P}}`
    in the elastic slope :math:`B(s)`. Defaults to ``2.24``.

:option:`ElasticSlope_c1`
    Constant term :math:`c_1` subtracted in the elastic slope
    :math:`B(s) = 2\,[\,b_A + b_B + c_0\,(s/s_0)^{\epsilon_\mathbb{P}} - c_1\,]`. Defaults to ``2.1``.

Photon dissociation (two-pion / VMD)
------------------------------------

These parameters control the :math:`\pi^+\pi^-` system produced when an incoming photon dissociates
through the vector-meson-dominance model.

:option:`TwoPionInterference`
    Treatment of the two-pion (:math:`\rho`/:math:`\omega`) system: ``0`` none, ``1`` :math:`\rho` only,
    ``2`` :math:`\rho+\omega`, ``3`` :math:`\rho+\omega+`\ continuum, ``4`` continuum only. Defaults to ``0``.

:option:`f_omega`
    Relative magnitude of the :math:`\omega` contribution in the :math:`\rho`--:math:`\omega` interference. Defaults to ``0.15``.

:option:`phi_omega`
    Relative phase of the :math:`\omega` contribution in the :math:`\rho`--:math:`\omega` interference. Defaults to ``-0.25``.

:option:`f_nr`
    Magnitude of the non-resonant (continuum) contribution to the two-pion mass distribution. Defaults to ``0.1675``.

:option:`Lambda_nr`
    Scale :math:`\Lambda_\text{nr}` (in GeV) of the non-resonant contribution. Defaults to ``0.15``.

:option:`delta_nr`
    Exponent :math:`\delta_\text{nr}` of the non-resonant contribution. Defaults to ``0.75``.

On-the-fly reweighting
----------------------

Sherpa can compute alternative event weights for variations of the MPI
model parameters on-the-fly, i.e. within a single event-generation run
that uses a nominal parameter set, without the need for dedicated runs at
each varied parameter point.

Reweighting is enabled by specifying a list of values for a parameter
instead of a single value.  The first entry of the list is the nominal
value, which is used for the actual event generation; the remaining
entries define the variations for which additional weights are computed.
For example,

.. code-block:: yaml

   AMISIC:
     SIGMA_ND_NORM: [0.5, 0.4, 0.6]

uses ``0.5`` as the nominal value of :option:`SIGMA_ND_NORM` and computes
two additional weights corresponding to ``0.4`` and ``0.6``.

The following MPI parameters support reweighting:

* :option:`SIGMA_ND_NORM`, the normalisation of the non-diffractive cross
  section, which through the constraint on the integrated interaction
  probability induces a rescaling of the effective matter-distribution
  widths;
* :option:`PT_0(ref)`, the IR regularisation scale;
* :option:`PT_Min(ref)`, the IR cutoff scale (only upward variations are possible);
* :option:`Eta`, the energy-scaling exponent, since varying it amounts to
  varying :option:`PT_0` and :option:`PT_Min` (only upward variations are possible);
* the matter-distribution parameters :option:`MATTER_RADIUS_1`,
  :option:`MATTER_RADIUS_2` and :option:`MATTER_FRACTION_1`, which are set
  in the :ref:`Remnants` block.

When several parameters are varied simultaneously, the variations are
matched by their position in the respective lists, i.e. the first variation
entry of each parameter forms one variation set, the second entries form
the next set, and so on.  If the number of variation entries differs
between parameters, missing values are automatically filled with the
corresponding nominal value.  For example,

.. code-block:: yaml

   AMISIC:
     SIGMA_ND_NORM: [0.5, 0.4, 0.6]
     PT_Min(ref):   [1.1, 1.2]

results in the two variation sets ``(0.5, 1.2)`` and ``(0.6, 1.1)``.

The same list syntax is used for the matter-distribution parameters in the
:ref:`Remnants` block and for the colour-reconnection
parameters in the :ref:`Colour_Reconnections` block.  Variations of all
three parameter groups are matched by list position and combined into a
single set of alternative weights.

The resulting variational event weights are labelled
``SoftPhysics.v1``, ``SoftPhysics.v2``, etc., where ``v1`` corresponds to
the first specified variation set.

.. _MI ISR parameters:

MI ISR parameters
=================

.. index:: MPI_PDF_SET
.. index:: MPI_PDF_SET_VERSIONS

The following two parameters can be used to overwrite the :ref:`ISR Parameters`
in the context of multiple interactions: ``MPI_PDF_SET``,
``MPI_PDF_SET_VERSIONS``.
