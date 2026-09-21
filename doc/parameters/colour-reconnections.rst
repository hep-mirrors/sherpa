.. _Colour_Reconnections:

********************
Colour_Reconnections
********************

The colour reconnections setup covers the non-perturbative reshuffling of
parton colours before the partons fragment into primordial hadrons.  In the
current implementation Sherpa collects all :math:`N` colourful partons 
entering colour reconnections and performs a sequence of :math:`N` iterations. 
In each iteration :math:`N` candidate reconnections between randomly chosen 
colour lines are proposed, and the one yielding the largest reduction of the 
total distance in momentum space is selected. The distance between two partons 
:math:`i` and :math:`j` is given by the Lund-inspired measure 
:math:`d_{ij} = p_ip_j-m_im_j`, where the gluon momenta are divided equally 
between both colours they carry.

The reshuffling of colour from old pairings :math:`\langle ij\rangle` and
:math:`\langle kl\rangle` to new pairings :math:`\langle il\rangle` and
:math:`\langle kj\rangle` is decided probabilistically, with

:math:`P(\langle ij\rangle\langle kl\rangle\to\langle il\rangle\langle kj\rangle) = R\cdot \left\{1-\exp[-\eta_Q^2\,(D_{ij}+D_{kl}-D_{il}-D_{kj})]\right\}`.

If the reconnection is accepted, the colour lines are swapped immediately, 
otherwise the original configuration is retained.  By default the logarithmic 
distance measure :math:`D_{ij}` is used (see :option:`PMODE` below).  For the 
power-law definition of the :math:`D_{ij}` we also calculate the average 
length of all colour connections :math:`ij` as 
:math:`\langle D_{ij}\rangle = \frac{1}{N}\sum_{ij}d_{ij}^\kappa`.
      
.. contents::
   :local:

.. index:: COLOUR_RECONNECTIONS:MODE
.. index:: COLOUR_RECONNECTIONS:PMODE
.. index:: COLOUR_RECONNECTIONS:Q_0
.. index:: COLOUR_RECONNECTIONS:ETA_Q
.. index:: COLOUR_RECONNECTIONS:RESHUFFLE
.. index:: COLOUR_RECONNECTIONS:KAPPA
	   
.. code-block:: yaml

   COLOUR_RECONNECTIONS:
     MODE:      On
     PMODE:     Log
     Q_0:       1.
     RESHUFFLE: 0.11

:option:`MODE` (default: On)
   Switches the colour reconnections on or off.

:option:`PMODE` (default: log)
   This switch defines how the distances of two partons in colour space are being
   calculated.  Available options define the distances :math:`D_{ij}` used in
   the decision whether colours are reshuffled as follows:
    
   * ``Log``: :math:`D_{ij} = \log(1+d_{ij}/Q_0^2)`
   * ``Power``: :math:`D_{ij} = \frac{d_{ij}^\kappa}{\langle D_{ij}\rangle}`

:option:`Q_0` (default: 1.)
   :math:`Q_0` in the logarithmic version of the momentum-space distance.

:option:`ETA_Q` (default: 0.63)
   :math:`\eta_Q` in the probability above.

:option:`RESHUFFLE` (default: 1/9)
   The colour suppression factor :math:`R` in the probability above.

:option:`KAPPA` (default: 1.)
   The exponent :math:`\kappa` in the equations above.
   
On-the-fly reweighting
----------------------

The colour-reconnection parameters :option:`ETA_Q` and :option:`RESHUFFLE` 
can be reweighted on-the-fly.

As for the MPI and beam-remnant parameters, reweighting is enabled by 
specifying a list of values instead of a single value.  The first entry is 
the nominal value used for event generation, the remaining entries define 
the variations.  For example,

.. code-block:: yaml

   COLOUR_RECONNECTIONS:
     MODE:  On
     ETA_Q: [0.6, 0.4, 0.8]

uses ``0.6`` as the nominal value of :option:`ETA_Q` and computes two 
additional weights corresponding to ``0.4`` and ``0.8``.

These colour-reconnection variations are matched by list position and 
combined with the variations of the MPI parameters (:ref:`Amisic`) and 
of the matter-distribution parameters (:ref:`Remnants`) (labelled as 
``SoftPhysics.v1``, ``SoftPhysics.v2``, etc.).

