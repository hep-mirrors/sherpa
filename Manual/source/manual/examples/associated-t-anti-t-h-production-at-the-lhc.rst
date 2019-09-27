.. _LHC_HTTBar:

Associated t anti-t H production at the LHC
===========================================



This set-up illustrates the interface to an external loop matrix element
generator as well as the possibility of specifying hard decays for
particles emerging from the hard interaction. The process generated
is the production of a Higgs boson in association with a top quark pair
from two light partons in the initial state. Each top quark decays into
an (anti-)bottom quark and a W boson. The W bosons in turn decay to either
quarks or leptons.

.. literalinclude:: /../../Examples/H_in_TTBar/LHC_TTH_MCatNLO/Sherpa.yaml
   :language: yaml

Things to notice:

* The virtual matrix elements can be interfaced either from OpenLoops
  or from :cite:`Reina2001sf`, :cite:`Reina2001bc`,
  :cite:`Dawson2002tg`, :cite:`Dawson2003zu`.  In the latter case, the
  shared library necessary for running this setup is built using
  `scons -f SConstruct-TTH <http://www.scons.org/>`_.

* The top quarks are stable in the hard matrix elements.
  They are decayed using the internal decay module, indicated by
  the settings in the :option:`HARD_DECAYS` and :option:`PARTICLE_DATA` blocks.

* Widths of top and Higgs are set to 0 for the matrix element calculation.
  A kinematical Breit-Wigner distribution will be imposed a-posteriori in the
  decay module.

* The Yukawa coupling of the b-quark has been set to a non-zero value to
  allow the H->bb decay channel even despite keeping the b-quark massless for a
  five-flavour-scheme calculation.

* Higgs BRs are not included in the cross section
  (:option:`Apply_Branching_Ratios: false`) as they will be LO only and not
  include loop-induced decays
