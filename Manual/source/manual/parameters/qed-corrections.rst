.. _QED Corrections:

***************
QED corrections
***************

Higher order QED corrections are effected both on hard interaction
and, upon their formation, on each hadron's subsequent decay. The
Photons :cite:`Schonherr2008av` module is called in both cases for
this task. It employes a YFS-type resummation :cite:`Yennie1961ad` of
all infrared singular terms to all orders and is equipped with
complete first order corrections for the most relevant cases (all
other ones receive approximate real emission corrections built up by
Catani-Seymour splitting kernels).

.. contents::
   :local:

.. _General Switches:

General Switches
================

The relevant switches to steer the higher order QED corrections are
collected in the :option:`YFS` settings group and are modified like
this:

.. code-block:: yaml

   YFS:
     <option1>: <value1>
     <option2>: <value2>
     ...

The options are

.. contents::
   :local:

.. _MODE:

MODE
----

.. index:: MODE

The keyword ``MODE`` determines the mode of operation of Photons.
``MODE: None`` switches Photons off.  Consequently, neither the hard
interaction nor any hadron decay will be corrected for soft or hard
photon emission. ``MODE: Soft`` sets the mode to "soft only", meaning
soft emissions will be treated correctly to all orders but no hard
emission corrections will be included. With ``MODE: Full`` these hard
emission corrections will also be included up to first order in
alpha_QED. This is the default setting.

.. _USE_ME:

USE_ME
------

.. index:: USE_ME

The switch ``USE_ME`` tells Photons how to correct hard emissions to
first order in alpha_QED. If ``USE_ME: 0``, then Photons will use
collinearly approximated real emission matrix elements. Virtual
emission matrix elements of order alpha_QED are ignored. If, however,
YFS_USE_ME=1, then exact real and/or virtual emission matrix elements
are used wherever possible. These are presently available for V->FF,
V->SS, S->FF, S->SS, S->Slnu, S->Vlnu type decays, Z->FF decays and
leptonic tau and W decays. For all other decay types general
collinearly approximated matrix elements are used. In both approaches
all hadrons are treated as point-like objects. The default setting is
``USE_ME: 1``. This switch is only effective if ``MODE: Full``.

.. _IR_CUTOFF:

IR_CUTOFF
---------

.. index:: IR_CUTOFF

``IR_CUTOFF`` sets the infrared cut-off dividing the real emission in
two regions, one containing the infrared divergence, the other the
"hard" emissions.  This cut-off is currently applied in the rest frame
of the multipole of the respective decay. It also serves as a minimum
photon energy in this frame for explicit photon generation for the
event record. In contrast, all photons below with energy less than
this cut-off will be assumed to have negligible impact on the
final-state momentum distributions. The default is ``IR_CUTOFF: 1E-3``
(GeV). Of course, this switch is only effective if Photons is switched
on, i.e. ``MODE`` is not set to ``None``.

.. _QED Corrections to the Hard Interaction:

QED Corrections to the Hard Interaction
=======================================

The switches to steer QED corrections to the hard scattering are
collected in the :option:`ME_QED` settings group and are modified like
this:

.. code-block:: yaml

   ME_QED:
     <option1>: <value1>
     <option2>: <value2>
     ...

The following options can be customised:

.. contents::
   :local:

.. _ENABLED:

ENABLED
-------

.. index:: ENABLED

``ENABLED: false`` turns the higher order QED corrections to the
matrix element off. The default is :option:`true`. Switching QED
corrections to the matrix element off has no effect on :ref:`QED
Corrections to Hadron Decays`.  The QED corrections to the matrix
element will only be effected on final state not strongly interacting
particles. If a resonant production subprocess for an unambiguous
subset of all such particles is specified via the process declaration
(cf. :ref:`Processes`) this can be taken into account and dedicated
higher order matrix elements can be used (if ``YFS: { MODE: Full,
USE_ME: 1 }``).

.. _CLUSTERING_ENABLED:

CLUSTERING_ENABLED
------------------

.. index:: CLUSTERING_ENABLED

``CLUSTERING_ENABLED: false`` switches the phase space point dependent
identification of possible resonances within the hard matrix element
on or off, respectively. The default is :option:`true`.  Resonances
are identified by recombining the electroweak final state of the
matrix element into resonances that are allowed by the model.
Competing resonances are identified by their on-shell-ness, i.e.  the
distance of the decay product's invariant mass from the nominal
resonance mass in units of the resonance width.

.. _CLUSTERING_THRESHOLD:

CLUSTERING_THRESHOLD
--------------------

.. index:: CLUSTERING_THRESHOLD

Sets the maximal distance of the decay product invariant mass from the
nominal resonance mass in units of the resonance width in order for the
resonance to be identified. The default is
:option:`CLUSTERING_THRESHOLD: 10.0`.

.. _QED Corrections to Hadron Decays:

QED Corrections to Hadron Decays
================================

If the Photons module is switched on, all hadron decays are corrected for higher
order QED effects.
