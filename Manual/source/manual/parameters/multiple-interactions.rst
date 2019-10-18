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

.. _MI ISR parameters:

MI ISR parameters
=================

.. index:: MPI_PDF_SET
.. index:: MI_PDF_SET_VERSIONS

The following two parameters can be used to overwrite the :ref:`ISR Parameters`
in the context of multiple interactions: ``MPI_PDF_SET``,
``MPI_PDF_SET_VERSIONS``.

.. _TURNOFF:

TURNOFF
=======

.. index:: TURNOFF

Specifies the transverse momentum turnoff in GeV.

.. _SCALE_MIN:

SCALE_MIN
=========

.. index:: SCALE_MIN

Specifies the transverse momentum integration cutoff in GeV.


.. _PROFILE_FUNCTION:

PROFILE_FUNCTION
================

.. index:: PROFILE_FUNCTION

Specifies the hadron profile function. The possible values are
:option:`Exponential`, :option:`Gaussian` and :option:`Double_Gaussian`.
For the double gaussian profile, the relative core size and relative
matter fraction can be set using :ref:`PROFILE_PARAMETERS`.

.. _PROFILE_PARAMETERS:

PROFILE_PARAMETERS
==================

.. index:: PROFILE_PARAMETERS

The potential parameters for hadron profile functions, see
:ref:`PROFILE_FUNCTION`. For double gaussian profiles there are
two parameters, corresponding to the relative core size and relative
matter fraction.

.. _REFERENCE_SCALE:

REFERENCE_SCALE
===============

.. index:: REFERENCE_SCALE

Specifies the centre-of-mass energy at which the transverse momentum
integration cutoff is used as is, see :ref:`SCALE_MIN`.
This parameter should not be changed by the user. The default is
``1800``, corresponding to Tevatron Run I energies.

.. _RESCALE_EXPONENT:

RESCALE_EXPONENT
================

.. index:: RESCALE_EXPONENT

Specifies the rescaling exponent for fixing the transverse momentum
integration cutoff at centre-of-mass energies different from the
reference scale, see :ref:`SCALE_MIN`, :ref:`REFERENCE_SCALE`.

.. _TURNOFF_EXPONENT:

TURNOFF_EXPONENT
================

.. index:: TURNOFF_EXPONENT

Specifies the rescaling exponent for fixing the transverse momentum
turnoff at centre-of-mass energies different from the reference scale,
see :ref:`TURNOFF`, :ref:`REFERENCE_SCALE`.

.. _SIGMA_ND_FACTOR:

SIGMA_ND_FACTOR
===============

.. index:: SIGMA_ND_FACTOR

Specifies the factor to scale the non-diffractive cross section
calculated in the MPI initialisation.

.. _MI_RESULT_DIRECTORY:

MI_RESULT_DIRECTORY
===================

.. index:: MI_RESULT_DIRECTORY

Specifies the name of the directory where the MPI grid is stored. The
default comprises the beam particles, their energies and the PDF used.
In its default value, this information safeguards against using unsuitable
grids for the current calculation.

.. _MI_RESULT_DIRECTORY_SUFFIX:

MI_RESULT_DIRECTORY_SUFFIX
==========================

.. index:: MI_RESULT_DIRECTORY_SUFFIX

Supplements the default directory name for the MPI grid with a suffix.
