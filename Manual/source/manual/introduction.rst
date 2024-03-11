.. _Introduction:

############
Introduction
############


Sherpa is a Monte Carlo event generator for the Simulation of
High-Energy Reactions of Particles in lepton-lepton, lepton-photon,
photon-photon, lepton-hadron and hadron-hadron collisions.  This
document provides information to help users understand and apply
Sherpa for their physics studies. The event generator is introduced,
in broad terms, and the installation and running of the program are
outlined. The various options and parameters specifying the program
are compiled, and their meanings are explained. This document does not
aim at giving a complete description of the physics content of Sherpa.
To this end, the authors refer the reader to the original publications,
:cite:`Gleisberg2008ta` and :cite:`Bothmann:2019yzt`.


.. contents::
   :local:


.. _Introduction to Sherpa:

**********************
Introduction to Sherpa
**********************



`Sherpa <http://projects.hepforge.org/sherpa/>`_ :cite:`Bothmann:2019yzt` is a Monte Carlo event generator that provides
complete hadronic final states in simulations of high-energy particle
collisions. The produced events may be passed into detector
simulations used by the various experiments.  The entire code has been
written in C++, like its competitors `Herwig++
<http://projects.hepforge.org/herwig/>`_ :cite:`Bahr2008pv` and
`Pythia 8 <https://pythia.org/>`_
:cite:`Sjostrand2007gs`.

Sherpa simulations can be achieved for the following types of
collisions:

* for lepton--lepton collisions, as explored by the CERN LEP
  experiments,
* for lepton--photon collisions,
* for photon--photon collisions with both photons either resolved or
  unresolved,
* for deep-inelastic lepton-hadron scattering, as investigated
  by the HERA experiments at DESY, and,
* in particular, for hadronic interactions as studied at the
  Fermilab Tevatron or the CERN LHC.



The list of physics processes that can be simulated with Sherpa covers
all reactions in the Standard Model. Other models can be implemented
either using Sherpa's own model syntax, or by using the generic
interface :cite:`Hoeche2014kca` to the UFO output
:cite:`Degrande2011ua` :cite:`Darme:2023jdn` of `FeynRules`_
:cite:`Christensen2008py`, :cite:`Christensen2009jx`.  The Sherpa
program owes this versatility to the two inbuilt matrix-element
generators, AMEGIC++ and `Comix <http://comix.freacafe.de>`_, and to
it's phase-space generator Phasic :cite:`Krauss2001iv`, which
automatically calculate and integrate tree-level amplitudes for the
implemented models.  This feature enables Sherpa to be used as a
cross-section integrator and parton-level event generator as well.
This aspect has been extensively tested, see e.g.
:cite:`Gleisberg2003bi`, :cite:`Hagiwara2005wg`.


As a second key feature of Sherpa the program provides an
implementation of the merging approaches of :cite:`Hoeche2009rj` and
:cite:`Gehrmann2012yg`, :cite:`Hoeche2012yf`.  These algorithms yield
improved descriptions of multijet production processes, which
copiously appear at lepton-hadron colliders like HERA
:cite:`Carli2009cg`, or hadron-hadron colliders like the Tevatron and
the LHC, :cite:`Krauss2004bs`, :cite:`Krauss2005nu`,
:cite:`Gleisberg2005qq`, :cite:`Hoeche2009xc`.  An older approach,
implemented in previous versions of Sherpa and known as the CKKW
technique :cite:`Catani2001cc`, :cite:`Krauss2002up`, has been
compared in great detail in :cite:`Alwall2007fs` with other
approaches, such as the MLM merging prescription :cite:`Mangano2001xp`
as implemented in Alpgen :cite:`Mangano2002ea`, Madevent
:cite:`Stelzer1994ta`, :cite:`Maltoni2002qb`, or Helac
:cite:`Kanaki2000ey`, :cite:`Papadopoulos2005ky` and the CKKW-L
prescription :cite:`Lonnblad2001iq`, :cite:`Lavesson2005xu` of Ariadne
:cite:`Lonnblad1992tz`.

This manual contains all information necessary to get started with
Sherpa as quickly as possible. It lists options and switches of
interest for steering the simulation of various physics aspects of the
collision.  It does not describe the physics simulated by Sherpa or
the underlying structure of the program.  Many external codes can be
linked with Sherpa. This manual explains how to do this, but it does
not contain a description of the external programs.  You are
encouraged to read their corresponding documentations, which are
referenced in the text. If you use external programs with Sherpa, you
are encouraged to cite them accordingly.

The `MCnet Guidelines
<http://www.montecarlonet.org/index.php?p=Publications/Guidelines>`_
apply to Sherpa. You are kindly asked to cite :cite:`Bothmann:2019yzt`
if you have used the program in your work.

The Sherpa authors strongly recommend the study of the manuals and
many excellent publications on different aspects of event generation
and physics at collider experiments written by other event generator
authors.

This manual is organized as follows: in :ref:`Basic structure` the
modular structure intrinsic to Sherpa is introduced. :ref:`Getting
started` contains information about and instructions for the
installation of the package. There is also a description of the steps
that are needed to run Sherpa and generate events.  The :ref:`Input
structure` is then discussed, and the ways in which Sherpa can be
steered are explained.  All parameters and options are discussed in
:ref:`Parameters`.  Advanced :ref:`Tips and Tricks` are detailed, and
some options for :ref:`Customization` are outlined for those more
familiar with Sherpa.  There is also a short description of the
different :ref:`Examples` provided with Sherpa.

The construction of Monte Carlo programs requires several assumptions,
approximations and simplifications of complicated physics aspects. The
results of event generators should therefore always be verified and
cross-checked with results obtained by other programs, and they should
be interpreted with care and common sense.

.. _Basic structure:

***************
Basic structure
***************


Sherpa is a modular program. This reflects the paradigm of Monte Carlo
event generation, with the full simulation is split into well defined
event phases, based on QCD factorization theorems. Accordingly, each
module encapsulates a different aspect of event generation for
high-energy particle reactions. It resides within its own namespace
and is located in its own subdirectory of the same name. The main
module called ``SHERPA`` steers the interplay of all modules -- or
phases -- and the actual generation of the events.
Altogether, the following modules are currently distributed with the
Sherpa framework:

ATOOLS
  This is the toolbox for all other modules. Since the Sherpa
  framework does not rely on CLHEP etc., the ATOOLS contain classes
  with mathematical tools like vectors and matrices, organization
  tools such as read-in or write-out devices, and physics tools like
  particle data or classes for the event record.

METOOLS
  In this module some general methods for the evaluation of helicity
  amplitudes have been accumulated.  They are used in AMEGIC++, the
  EXTRA_XS module, and the new matrix-element generator Comix.  This
  module also contains helicity amplitudes for some generic matrix
  elements, that are, e.g., used by HADRONS++. Further, METOOLS also
  contains a simple library of tensor integrals which are used in the
  PHOTONS++ matrix element corrections.

BEAM
  This module manages the treatment of the initial beam spectra
  for different colliders. The three options which are currently
  available include a monochromatic beam, which requires no extra
  treatment, photon emission in the Equivalent Photon Approximation
  (EPA) and - for the case of an electron collider - laser
  backscattering off the electrons, leading to photonic initial
  states.

PDF
  The PDF module provides access to various parton density
  functions (PDFs) for the proton and the photon. In addition, it
  hosts an interface to the `LHAPDF
  <http://projects.hepforge.org/lhapdf>`_ package, which makes a full
  wealth of PDFs available. An (analytical) electron structure
  function is supplied in the PDF module as well.

MODEL
  This module sets up the physics model for the simulation.  It
  initializes particle properties, basic physics parameters (coupling
  constants, mixing angles, etc.) and the set of available interaction
  vertices (Feynman rules). By now, there exist explicit
  implementations of the Standard Model (SM), its Minimal
  Supersymmetric extension (MSSM), the ADD model of large extra
  dimensions, and a comprehensive set of operators parametrizing
  anomalous triple and quartic electroweak gauge boson couplings. An
  Interface to `FeynRules`_ is also
  available.

EXTRA_XS
  In this module a (limited) collection of analytic expressions
  for simple 2->2 processes within the SM are provided together with
  classes embedding them into the Sherpa framework. This also includes
  methods used for the definition of the starting conditions for
  parton-shower evolution, such as colour connections and the hard
  scale of the process.

AMEGIC++
  AMEGIC++ :cite:`Krauss2001iv` is Sherpa's original
  matrix-element generator. It employs the method of helicity
  amplitudes :cite:`Kleiss1985yh`, :cite:`Ballestrero1992dv` and works
  as a generator, which generates generators: During the
  initialization run the matrix elements for a given set of processes,
  as well as their specific phase-space mappings are created by
  AMEGIC++.  Corresponding C++ sourcecode is written to disk and
  compiled by the user using the ``makelibs`` script.
  The produced libraries are linked to the
  main program automatically in the next run and used to calculate
  cross sections and to generate weighted or unweighted events.
  AMEGIC++ has been tested for multi-particle production in the
  Standard Model :cite:`Gleisberg2003bi`. Its MSSM implementation has
  been validated in :cite:`Hagiwara2005wg`.

COMIX
  `Comix <http://comix.freacafe.de>`_ is a multi-leg tree-level
  matrix element generator, based on the color dressed Berends-Giele
  recursive relations :cite:`Duhr2006iq`.  It employs a new algorithm
  to recursively compute phase-space weights.  The module is a useful
  supplement to older matrix element generators like AMEGIC++ in
  the high multiplicity regime. Due to the usage of colour sampling it
  is particularly suited for an interface with parton shower
  simulations and can hence be easily employed for the ME-PS merging
  within Sherpa. It is Sherpa's default large multiplicity matrix
  element generator for the Standard Model.

PHASIC++
  All base classes dealing with the Monte Carlo phase-space
  integration are located in this module. For the evaluation of the
  initial-state (laser backscattering, initial-state radiation) and
  final-state integrals, the adaptive multi-channel method of
  :cite:`Kleiss1994qy`, :cite:`Berends1994pv` is used by default
  together with a Vegas optimization :cite:`Lepage1980dq` of the
  single channels. In addition, final-state integration accomplished
  by Rambo :cite:`Kleiss1985gy`, Sarge :cite:`Draggiotis2000gm` and
  HAAG :cite:`vanHameren2002tc` is supported.

CSSHOWER++
  This is the module hosting Sherpa's default parton shower,
  which was published in :cite:`Schumann2007mg`.  The corresponding
  shower model was originally proposed in :cite:`Nagy2005aa`,
  :cite:`Nagy2006kb`.  It relies on the factorisation of real-emission
  matrix elements in the CS subtraction framework
  :cite:`Catani1996vz`, :cite:`Catani2002hc`.  There exist four
  general types of CS dipole terms that capture the complete infrared
  singularity structure of next-to-leading order QCD amplitudes.  In
  the large-N_C limit, the corresponding splitter and spectator
  partons are always adjacent in colour space. The dipole functions
  for the various cases, taken in four dimensions and averaged over
  spins, are used as shower splitting kernels.

DIRE
  This is the module hosting Sherpa's alternative parton shower
  :cite:`Hoche2015sya`.  In the Dire model, the ordering variable
  exhibits a symmetry in emitter and spectator momenta, such that the
  dipole-like picture of the evolution can be re-interpreted as a
  dipole picture in the soft limit. At the same time, the splitting
  functions are regularized in the soft anti-collinear region using
  partial fractioning of the soft eikonal in the Catani-Seymour
  approach :cite:`Catani1996vz`, :cite:`Catani2002hc`. They are then
  modified to satisfy the sum rules in the collinear limit. This leads
  to an invariant formulation of the parton-shower algorithm, which is
  in complete analogy to the standard DGLAP case, but generates the
  correct soft anomalous dimension at one-loop order.

AMISIC++
  AMISIC++ contains classes for the simulation of multiple
  parton interactions according to :cite:`Sjostrand1987su`. In Sherpa
  the treatment of multiple interactions has been extended by allowing
  for the simultaneous evolution of an independent parton shower in
  each of the subsequent (semi-)hard collisions. 

REMNANTS:
  REMNANTS contains classes for the simulation of the beam remnants,
  including in particular the spatial form of the matter distribution
  which is relevant for the underlying event, and the treatment of
  the intrinsic transverse momentum.

RECONNECTIONS:
  RECONNECTIONS handles the colour reconnections preceding the hadronization.
  This module will experience future refinements.
  
AHADIC++
  AHADIC++ is Sherpa's hadronization package, for translating
  the partons (quarks and gluons) into primordial hadrons, to be
  further decayed in HADRONS++.  The algorithm bases on the cluster
  fragmentation ideas presented in :cite:`Gottschalk1982yt`,
  :cite:`Gottschalk1983fm`, :cite:`Webber1983if`,
  :cite:`Gottschalk1986bv` and implemented in the Herwig family of
  event generators.  The actual Sherpa implementation is based on
  :cite:`Chahal2022rid`.

HADRONS++
  HADRONS++ is the module for simulating hadron and tau-lepton
  decays.  The resulting decay products respect full spin correlations
  (if desired).  Several matrix elements and form-factor models have
  been implemented, such as the Kühn-Santamaría model, form-factor
  parametrizations from Resonance Chiral Theory for the tau and form
  factors from heavy quark effective theory or light cone sum rules
  for hadron decays.

PHOTONS++
  The PHOTONS++ module holds routines to add QED
  radiation to hadron and tau-lepton decays. This has been achieved by
  an implementation of the YFS algorithm :cite:`Yennie1961ad`. The
  structure of PHOTONS++ is such that the formalism can be
  extended to scattering processes and to a systematic improvement to
  higher orders of perturbation theory :cite:`Schonherr2008av`.  The
  application of PHOTONS++ therefore accounts for corrections
  that usually are added by the application of PHOTOS
  :cite:`Barberio1993qi` to the final state.

SHERPA
  Finally, SHERPA is the steering module that initializes,
  controls and evaluates the different phases during the entire
  process of event generation.  All routines for the combination of
  truncated showers and matrix elements, which are independent of the
  specific matrix element and parton shower generator are found in
  this module.


The actual executable of the Sherpa generator can be found in the
subdirectory ``<prefix>/bin/`` and is called ``Sherpa``. To run the
program, input files have to be provided in the current working
directory or elsewhere by specifying the corresponding path, see
:ref:`Input structure`. All output files are then written to this
directory as well.

.. _Feynrules: http://feynrules.irmp.ucl.ac.be/
