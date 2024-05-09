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
Catani-Seymour splitting kernels). The module is also equipped with 
an algorithm to allow any photons produced to split into charged 
particle pairs.

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

.. _PHOTON_SPLITTER_MODE:

PHOTON_SPLITTER_MODE
--------------------

.. index:: PHOTON_SPLITTER_MODE

The parameter :OPTION:`PHOTON_SPLITTER_MODE` determines which particles, if any, may be produced in 
photon splittings:

  :option:`0`
    All photon splitting functions are turned off.
  :option:`1`
    Photons may split into electron-positron pairs;
  :option:`2`
    muons;
  :option:`4`
    tau leptons;
  :option:`8`
    light hadrons up to ``PHOTON_SPLITTER_MAX_HADMASS``.

The settings are additive, e.g. ``PHOTON_SPLITTER_MODE: 3``
allows splittings into electron-positron and muon-antimuon pairs.
The default is ``PHOTON_SPLITTER_MODE: 15`` (all splittings turned on).
This parameter is of course only effective if the Photons module is 
switched on using the ``MODE`` keyword.

.. _PHOTON_SPLITTER_MAX_HADMASS:

PHOTON_SPLITTER_MAX_HADMASS
---------------------------

.. index:: PHOTON_SPLITTER_MAX_HADMASS

``PHOTON_SPLITTER_MAX_HADMASS`` sets the mass (in GeV) of the heaviest 
hadron which may be produced in photon splittings. Note that vector 
splitting functions are currently not implemented: only fermions, 
scalars and pseudoscalars up to this cutoff will be considered. 
The default is 0.5 GeV.

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

QED Corrections for Lepton-Lepton Collisions
================================

The YFS resummation can be enabled for lepton-lepton scattering
by setting :option:`MODE` to either :option:`ISR` or :option:`ISRFSR`.


The options are

.. contents::
   :local:


.. _Beta

BETA
------------------
Higher order matrix element corrections can be
included by setting :option:`BETA` to either
:option:`1/2/3` to the desired order of accuray.
For example :option:`BETA: 0` disables all higer-order
corrections, while :option:`BETA: 3` includes QED
corrections up to :math:`\alpha^3L^3` in the exclusive
exponentiation (EEX) :cite:`Jadach:2000ir`


.. _Coulomb:

COULOMB
------------------
The Coulomb threshold corrections :cite:`Bardin:1993mc` :cite:`Fadin:1993kg`
to the :math:`W^+W^-` threshold can be included with :option:`COULOMB: True`.
Double counting of the virtual corrections with the YFS form-factor
is avoided by using analytical subtraction in the threshold limit :cite:`Krauss:2022ajk`.