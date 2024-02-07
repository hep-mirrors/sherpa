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
.. index:: Amisic:PT_0(ref)
.. index:: Amisic:PT_0(IR)
.. index:: Amisic:PT_Min(ref)
.. index:: Amisic:Eta
.. index:: Amisic:E(ref)
.. index:: Amisic:PT_Min
.. index:: Amisic:PT_0
.. index:: Amisic:MU_R_SCHEME
.. index:: Amisic:MU_R_FACTOR
.. index:: Amisic:MU_F_FACTOR
.. index:: Amisic:SIGMA_ND_NORM
.. index:: Amisic:nPT_bins
.. index:: Amisic:nMC_points
.. index:: Amisic:nS_bins
.. index:: Amisic:PomeronIntercept
.. index:: Amisic:PomeronSlope
.. index:: Amisic:TriplePomeronCoupling
.. index:: Amisic:ReggeonIntercept

Amisic can simulate the interaction of three different combinations of incoming particles:
proton--proton, photon--proton and photon--photon collision. The parameters for the simulation of photonic multiple
interactions can be found in :cite:`Schuler:1993wr`. It has several parameters to control the simulation of the
multiple-parton interactions, they are listed below. Each of these parameters has to be set in the
subsetting ``AMISIC``, like so

.. code-block:: yaml

   AMISIC:
     PT_0: 2.5

The usual rules for yaml structure apply, c.f. :ref:`Input structure`.

:option:`PT_0(ref)`
    Value :math:`p_\text{T,0}^\text{(ref)}` for the calculation of the IR regulator, see formula below. Defaults to ``2.2``.

:option:`PT_0(IR)`
    The absolute minimum of the IR regulator, see formula below. Defaults to ``0.5``.

:option:`PT_Min(ref)`
    Value :math:`p_\text{T,min}^\text{(ref)}` for the calculation of the IR cutoff, see formula below. Defaults to ``3``.

:option:`Eta`
    The pseudorapidity :math:`\eta` used to calculate the IR cutoff and regulator, :math:`p_\text{T,min}` and :math:`p_\text{T,0}`.
    Defaults to ``0.16``.

:option:`E(ref)`
    Reference energy to normalise the actual cms energy for the calculation of the IR cutoff and regulator.
    Defaults to ``7000``.

:option:`PT_Min`
    The IR cut-off for the 2->2 scatters. It is calculated as

    .. math::

        p_\text{T,min} = p_\text{T,min}^\text{(ref)} \left( \frac{E_\text{cms}}{E_\text{cms}^\text{(ref)}} \right)^{2\eta}

    but can also be set explicitly.

:option:`PT_0`
    IR regulator :math:`p_\text{T,0}` in the propagator and in the strong coupling. It is calculated as

    .. math::

        p_\text{T,0} = p_\text{T,0}^\text{(ref)} \left( \frac{E_\text{cms}}{E_\text{cms}^\text{(ref)}} \right)^{2\eta}

    but can also be set explicitly.

:option:`MU_R_SCHEME`
    Defaults to ``PT`` scheme. More schemes have yet to be added.

:option:`MU_R_FACTOR`
    Factor to scale the renormalisation scale :math:`\mu_R`, defaults to ``0.5``.

:option:`MU_F_FACTOR`
    Factor to scale the factorisation scale :math:`\mu_F`, defaults to ``1.0``.

:option:`SIGMA_ND_NORM`
    Specifies the factor to scale the non-diffractive cross section calculated in the MPI initialisation.
    Defaults to ``0.4``.

:option:`nPT_bins`
    Controls the number of bins for the numerical integration of

    .. math:: \int_{p_T^2}^{s/4} dp_T^2 \frac{d \sigma}{dp_T^2}

    Defaults to ``200``.

:option:`nMC_points`
    Number of points to estimate the the cross-section during the integration. The error should behave as
    :math:`\frac{1}{\sqrt{n_\text{MC}}}`. Defaults to ``1000``.

:option:`nS_bins`
    Number of points to sample in the center-of-mass energy :math:`\sqrt{s}`. This is only used if the energy is not
    fixed, i.e. in the case of EPA photons. Defaults to ``100``.

The total cross-section is calculated with

    .. math::

        \sigma_{tot} = X s^\epsilon + Y s^\eta

    where :math:`s` is the Mandelstam invariant.

:option:`PomeronIntercept`
    The parameter :math:`\epsilon` in the above equation, defaults to ``0.0808``.

:option:`ReggeonIntercept`
    The parameter :math:`\eta` in the above equation, defaults to ``-0.4525``.

The single- and double-diffractive cross-sections in the Regge picture have two free parameters:

:option:`PomeronSlope`
    The parameter :math:`\alpha^\prime`, default is ``0.25``.

:option:`TriplePomeronCoupling`
    The parameter :math:`g_{3\mathbb{P}}` at an input scale of 20 GeV, given in :math:`\text{mb}^{-0.5}`, with default ``0.318``.

.. _MI ISR parameters:

MI ISR parameters
=================

.. index:: MPI_PDF_SET
.. index:: MPI_PDF_SET_VERSIONS

The following two parameters can be used to overwrite the :ref:`ISR Parameters`
in the context of multiple interactions: ``MPI_PDF_SET``,
``MPI_PDF_SET_VERSIONS``.
