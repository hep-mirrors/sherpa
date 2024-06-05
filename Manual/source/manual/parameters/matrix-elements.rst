
.. _Matrix Elements:

***************
Matrix elements
***************


The following parameters are used to steer the matrix element calculation setup. To learn how to specify the hard scattering process and further process-specific options in its calculation, please also refer to :ref:`Processes`.

.. contents::
   :local:

.. _ME_GENERATORS:

ME_GENERATORS
=============

.. index:: ME_GENERATORS

The list of matrix element generators to be employed during the run.
When setting up hard processes, Sherpa calls these generators in order
to check whether either one is capable of generating the corresponding
matrix element. This parameter can also be set on the command line
using option :option:`-m`, see :ref:`Command line`.

The built-in generators are

:option:`Internal`
  Simple matrix element library, implementing a variety of 2->2 processes.

:option:`Amegic`
  The AMEGIC++ generator published under :cite:`Krauss2001iv`

:option:`Comix`
  The `Comix <http://comix.freacafe.de>`_ generator published under
  :cite:`Gleisberg2008fv`

It is possible to employ an external matrix element generator within
Sherpa.  For advice on this topic please contact the authors,
:ref:`Authors`.

.. _RESULT_DIRECTORY:

RESULT_DIRECTORY
================

.. index:: RESULT_DIRECTORY

This parameter specifies the name of the directory which is used by
Sherpa to store integration results and phase-space mappings. The
default is ``Results/``.  It can also be set using the command line
parameter :option:`-r`, see :ref:`Command line`. The directory will be
created automatically, unless the option
:option:`GENERATE_RESULT_DIRECTORY: false` is specified.  Its location
is relative to a potentially specified input path, see :ref:`Command
line`.

.. _EVENT_GENERATION_MODE:

EVENT_GENERATION_MODE
=====================

.. index:: EVENT_GENERATION_MODE
.. index:: OVERWEIGHT_THRESHOLD

This parameter specifies the event generation mode.  It can also be
set on the command line using option :option:`-w`, see :ref:`Command
line`.  The three possible options are:

:option:`Weighted`
  (alias :option:`W`) Weighted events.

:option:`Unweighted`
  (alias :option:`U`)
  Events with constant weight, which have been unweighted against the
  maximum determined during phase space integration.  In case of rare
  events with ``w > max`` the parton level event is repeated
  ``floor(w/max)`` times and the remainder is unweighted.  While this
  leads to unity weights for all events it can be misleading since the
  statistical impact of a high-weight event is not accounted for. In
  the extreme case this can lead to a high-weight event looking like a
  significant bump in distributions (in particular after the effects
  of the parton shower).

:option:`PartiallyUnweighted`
  (alias :option:`P`)
  Identical to :option:`Unweighted` events, but if the weight exceeds
  the maximum determined during the phase space integration, the event
  will carry a weight of ``w/max`` to correct for that. This is the
  recommended option to generate unweighted events and the default
  setting in Sherpa.

For :option:`Unweighted` and :option:`PartiallyUnweighted` events the user may
set :option:`OVERWEIGHT_THRESHOLD: <maxweight>` to cap the maximal over-weight
``w/max`` taken into account.


.. _SCALES:

SCALES
======

.. index:: SCALES

This parameter specifies how to compute the renormalisation and
factorisation scale and potential additional scales.

.. note::

   In a setup with the parton shower enabled, it is strongly recommended to
   leave this at its default value, :option:`METS`, and to instead customise the
   :option:`CORE_SCALE` setting as described in :ref:`METS scale setting with multiparton core processes`.

.. contents::
   :local:

Sherpa provides several built-in scale setting schemes. For each
scheme the scales are then set using expressions understood by the
:ref:`Interpreter`.  Each scale setter's syntax is

.. code-block:: yaml

   SCALES: <scale-setter>{<scale-definition>}

to define a single scale for both the factorisation and renormalisation scale.
They can be set to different values using

.. code-block:: yaml

   SCALES: <scale-setter>{<fac-scale-definition>}{<ren-scale-definition>}

In parton shower matched/merged calculations a third perturbative
scale is present, the resummation or parton shower starting scale. It
can be set by the user in the third argument like

.. code-block:: yaml

   SCALES: <scale-setter>{<fac-scale-definition>}{<ren-scale-definition>}{<res-scale-definition>}

If the final state of your hard scattering process contains QCD
partons, their kinematics fix the resummation scale for subsequent
emissions (cf. the description of the :option:`METS` scale setter
below).  With the CS Shower, you can instead specify your own
resummation scale also in such a case: Set ``SHOWER:RESPECT_Q2: true``
and use the third argument to specify your resummation scale as above.

.. note::

   For all scales their squares have to be given. See
   :ref:`Predefined scale tags` for some predefined scale tags.

More than three scales can be set as well to be subsequently used,
e.g.  by different couplings, see :ref:`COUPLINGS`.

.. _Scale setters:

Scale setters
-------------


The scale setter options which are currently available are

:option:`METS`
  ``METS`` is the default scale scheme in Sherpa and employed
  for multi-leg merging, both at leading and next-to-leading order.
  Since it is important and complex at the same time, it will be described in
  detail in the next section.


:option:`VAR`
  The variable scale setter is the simplest scale setter available. Scales
  are simply specified by additional parameters in a form which is understood
  by the internal interpreter, see :ref:`Interpreter`. If, for example the invariant
  mass of the lepton pair in Drell-Yan production is the desired scale,
  the corresponding setup reads

  .. code-block:: yaml

     SCALES: VAR{Abs2(p[2]+p[3])}

  Renormalisation and factorisation scales can be chosen differently.
  For example in Drell-Yan + jet production one could set

  .. code-block:: yaml

     SCALES: VAR{Abs2(p[2]+p[3])}{MPerp2(p[2]+p[3])}

:option:`FASTJET`
  This scale setter can be used to set a scale based on jet-, rather
  than parton-momenta, using `FastJet <http://www.fastjet.fr>`_.

  The final state parton configuration is first clustered using
  FastJet and resulting jet momenta are then added back to the list of
  non strongly interacting particles. The numbering of momenta
  therefore stays effectively the same as in standard Sherpa, except
  that final state partons are replaced with jets, if applicable (a
  parton might not pass the jet criteria and get "lost"). In
  particular, the indices of the initial state partons and all EW
  particles are unaffected. Jet momenta can then be accessed as
  described in :ref:`Predefined scale tags` through the identifiers
  ``p[i]``, and the nodal values of the clustering sequence can be
  used through ``MU_n2``.  The syntax is


  .. code-block:: yaml

     SCALES: FASTJET[<jet-algo-parameter>]{<scale-definition>}

  Therein the parameters of the jet algorithm to be used to define the
  jets are given as a comma separated list of

  * the jet algorithm ``A:kt,antikt,cambridge,siscone`` (default
    ``antikt``)

  * phase space restrictions, i.e. ``PT:<min-pt>``, ``ET:<min-et>``,
    ``Eta:<max-eta>``, ``Y:<max-rap>`` (otherwise unrestricted)

  * radial parameter ``R:<rad-param>`` (default ``0.4``)

  * f-parameter for Siscone ``f:<f-param>`` (default ``0.75``)

  * recombination scheme ``C:E,pt,pt2,Et,Et2,BIpt,BIpt2``
    (default ``E``)

  * b-tagging mode ``B:0,1,2`` (default ``0``)
    This parameter, if specified different from its default 0, allows
    to use b-tagged jets only, based on the parton-level constituents of the jets.
    There are two options: With ``B:1`` both b and anti-b quarks are
    counted equally towards b-jets, while for ``B:2`` they are added with a
    relative sign as constituents, i.e. a jet containing b and anti-b is not tagged.

  * scale setting mode ``M:0,1`` (default ``1``) It is possible to
    specify multiple scale definition blocks, each enclosed in curly
    brackets. The scale setting mode parameter then determines, how
    those are interpreted: In the ``M:0`` case, they specify
    factorisation, renormalisation and resummation scale separately in
    that order.  In the ``M:1`` case, the ``n`` given scales are used
    to calculate a mean scale such that
    :math:`\alpha_s^n(\mu_\text{mean})=\alpha_s(\mu_1)\dots\alpha_s(\mu_n)`
    This scale is then used for factorisation, renormalisation and
    resummation scale.

  Consider the example of lepton pair production in association with jets. The
  following scale setter

  .. code-block:: yaml

     SCALES: FASTJET[A:kt,PT:10,R:0.4,M:0]{sqrt(PPerp2(p[4])*PPerp2(p[5]))}

  reconstructs jets using the kt-algorithm with R=0.4 and a minimum
  transverse momentum of 10 GeV. The scale of all strong couplings is
  then set to the geometric mean of the hardest and second hardest
  jet. Note ``M:0``.

  Similarly, in processes with multiple strong couplings, their
  renormalisation scales can be set to different values, e.g.


  .. code-block:: yaml

     SCALES: FASTJET[A:kt,PT:10,R:0.4,M:1]{PPerp2(p[4])}{PPerp2(p[5])}

  sets the scale of one strong coupling to the transverse momentum of
  the hardest jet, and the scale of the second strong coupling to the
  transverse momentum of second hardest jet. Note ``M:1`` in this
  case.

  The additional tags :samp:`{MU_22}` .. :samp:`{MU_n2}`
  (n=2..njet+1), hold the nodal values of the jet clustering in
  descending order.

  Please note that currently this type of scale setting can only be done within
  the process block (:ref:`Processes`) and not within the (me) section.

..
..   :option:`QCD`
     The matrix element is clustered onto a core 2->2 configuration using a
     k_T-type algorithm with recombination into on-shell partons.
     Scales are defined as the minimum of the largest transverse momentum
     during clustering and the lowest invariant mass in the core process.


:option:`VBF`                                                                                                                                                                                                                                                                                                                
  Very similar to the :option:`METS` scale setter and thus also applicable in multi-leg merged setups, but
  catering specifically to topologies with two colour-separated parton lines like in VBF/VBS
  processes for the incoming quarks. 



.. _METS scale setting with multiparton core processes:

Scale setting in multi-parton processes (METS)
----------------------------------------------

.. index:: METS
.. index:: CORE_SCALE

``METS`` is the default scale setting in Sherpa, since it is employed
for multi-leg merging, both at leading and next-to-leading order.
It dynamically defines the three tags ``MU_F2``, ``MU_R2`` and
``MU_Q2`` as will be explained below. Those can then be employed in the
actual ``<scale-definition>`` in the scale setter. The default is

.. code-block:: yaml

   SCALES: METS{MU_F2}{MU_R2}{MU_Q2}

The tags may be omitted, i.e.

.. code-block:: yaml

   SCALES: METS

leads to an identical scale definition.

``METS`` is a very dynamic scheme and depends on two ingredients to construct
a scale that preserves the logarithmic accuracy of the parton evolution defined
by the parton shower:

1. A sequential recombination algorithm to cluster the multi-leg matrix element onto a core
   configuration (typically 2->2) using an inversion of the current parton shower.
   The clustered flavours are determined using run-time information from the matrix element
   generator. The clustering stops, when no combination
   is found that corresponds to a parton shower branching, or if
   two subsequent branchings are unordered in terms of the parton shower
   evolution parameter. That defines the core process.

2. A freely selectable scale in the core process, ``CORE_SCALE``.

These are then defined to calculate ``MU_R2`` from the core scale and the
individual clustering scales such that:

.. math::

   \alpha_s(\mu_{R}^2)^{n+k} = \alpha_s(\mu_{R,\text{core-scale}}^2)^k \alpha_s(k_{t,1}^2) \dots \alpha_s(k_{t,n}^2)

where :math:`n` is the order in strong coupling of the core process and :math:`k` is
the number of clusterings, :math:`k_{t,i}` are the relative transverse momenta
at each clustering.

The definition of ``MU_F2`` and ``MU_Q2`` are passed directly on from the core scale setter.

The functional form of the core scale can be defined by the user in the ``MEPS``
settings block as follows:
  
.. code-block:: yaml

   MEPS:
     CORE_SCALE: <core-scale-setter>{<core-fac-scale-definition>}{<core-ren-scale-definition>}{<core-res-scale-definition>}

Again, for core scale setters which define ``MU_F2``, ``MU_R2`` and
``MU_Q2`` the actual scale listing can be dropped.

Possible choices for the core scale setter are:

:option:`Default`
  The core scales are defined depending on the core process and as the list
  and functional form is regularly extended to more core processes, its definition
  is most easily seen in the `source code <https://gitlab.com/sherpa-team/sherpa/-/blob/master/PHASIC%2B%2B/Scales/Default_Core_Scale.C>`_

:option:`VAR`
  Variable core scale setter for free functional definition by the user.
  Syntax is identical to variable scale setter.

:option:`QCD`
  QCD core scale setter. Scales are set to harmonic mean of s, t and u. Only
  useful for 2->2 cores as alternative to the :option:`Default` core scale.

:option:`TTBar`
  Core scale setter for processes involving top quarks. Implementation details
  are described in Appendix C of :cite:`Hoeche2013mua`.

:option:`SingleTop`
  Core scale setter for single-top production in association with one jet.
  If the W is in the t-channel (s-channel), the squared scales are set to the
  Mandelstam variables ``t=2*p[0]*p[2]`` (``t=2*p[0]*p[1]``).

:option:`Photons`
  Core scale setter for photon(s)+jets production.
  Sets the following scales for the possible core processes:

    - :math:`\gamma\gamma`: :math:`\mu_f=\mu_r=\mu_q=m_{\gamma\gamma}`
    - :math:`\gamma j`: :math:`\mu_f=\mu_r=\mu_q=p_{\perp,\gamma}`
    - :math:`jj`: same as QCD core scale (harmonic mean of s, t, u)



Unordered cluster histories are by default not allowed. Instead, if during
clustering a new smaller scale is encountered, the previous maximal scale
will be used, or alternatively a user-defined scale specified, e.g.
  
.. code-block:: yaml

   MEPS:
     UNORDERED_SCALE: VAR{H_Tp2/sqr(N_FS-2)}

If instead you want to allow unordered histories you can also enable them with
``ALLOW_SCALE_UNORDERING: 1``.

Clusterings onto 2->n (n>2) configurations is possible and for complicated
processes can warrant the implementation of a custom core scale, cf. :ref:`Customization`.

Occasionally, users might encounter the warning message

.. code-block:: console

   METS_Scale_Setter::CalculateScale(): No CSS history for '<process name>' in <percentage>% of calls. Set \hat{s}.

As long as the percentage quoted here is not too high, this does not pose
a serious problem. The warning occurs when - based on the current colour
configuration and matrix element information - no suitable clustering is
found by the algorithm. In such cases the scale is set to the invariant mass
of the partonic process.

One final word of caution: The ``METS`` scale scheme might be subject to changes
to enable further classes of processes for merging in the future and integration
results might thus change slightly between different Sherpa versions.


.. _Custom scale implementation:

Custom scale implementation
---------------------------


When the flexibility of the :option:`VAR` scale setter above is not sufficient,
it is also possible to implement a completely custom scale scheme within Sherpa
as C++ class plugin. For details please refer to the :ref:`Customization`
section.

.. _Predefined scale tags:

Predefined scale tags
---------------------


There exist a few predefined tags to facilitate commonly used scale
choices or easily implement a user defined scale.

:option:`p[n]`
  Access to the four momentum of the nth particle. The initial state
  particles carry n=0 and n=1, the final state momenta start from
  n=2. Their ordering is determined by Sherpa's internal particle
  ordering and can be read e.g.  from the process names displayed at
  run time. Please note, that when building jets out of the final
  state partons first, e.g. through the ``FASTJET`` scale setter,
  these parton momenta will be replaced by the jet momenta ordered in
  transverse momenta. For example the process u ub -> e- e+ G G will
  have the electron and the positron at positions ``p[2]`` and
  ``p[3]`` and the gluons on positions ``p[4]`` and ``p[5]``. However,
  when finding jets first, the electrons will still be at ``p[2]`` and
  ``p[3]`` while the harder jet will be at ``p[4]`` and the softer one
  at ``p[5]``.

:option:`H_T2`
  Square of the scalar sum of the transverse momenta of
  all final state particles.

:option:`H_TM2`
  Square of the scalar sum of the transverse energies of
  all final state particles, i.e. contrary to ``H_T2`` ``H_TM2`` takes
  particle masses into account.

:option:`H_TY2(<factor>,<exponent>)`
  Square of the scalar sum of the transverse momenta of all final state particles
  weighted by their rapidity distance from the final state boost vector. Thus,
  takes the form

  .. code-block:: latex

     H_T^{(Y)} = sum_i pT_i exp [ fac |y-yboost|^exp ]

  Typical values to use would by ``0.3`` and ``1``.

:option:`H_Tp2`
  Scale setter for lepton-pair production in association with jets only,
  implements

  .. code-block:: latex

     H_T' = sqrt(m_ll^2 + pT(ll)^2) + sum_i pT_i (i not l)

:option:`DH_Tp2(<recombination-method>,<dR>)`
  Implements a version of ``H_Tp2`` which dresses charged particles first.
  The parameter ``<recombination-method>`` can take the following values:
  ``Cone``, ``kt``, ``CA`` or ``antikt``, while ``<dR>`` is
  the respective algorithm's angular distance parameter.

:option:`TAU_B2`
  Square of the beam thrust.

:option:`MU_F2, MU_R2, MU_Q2`
  Tags holding the values of the factorisation, renormalisation scale and
  resummation scale determined through backwards clustering in the
  ``METS`` scale setter.

:option:`MU_22, MU_32, ..., MU_n2`
  Tags holding the nodal values of the jet clustering in the ``FASTJET``
  scale setter, cf. :ref:`Scale setters`.




All of those objects can be operated upon by any operator/function known
to the :ref:`Interpreter`.

.. _Scale schemes for NLO calculations:

Scale schemes for NLO calculations
----------------------------------


For next-to-leading order calculations it must be guaranteed that the scale is
calculated separately for the real correction and the subtraction terms,
such that within the subtraction procedure the same amount is subtracted
and added back. Starting from version 1.2.2 this is the case for all
scale setters in Sherpa. Also, the definition of the scale must be
infrared safe w.r.t. to the radiation of an extra parton. Infrared safe
(for QCD-NLO calculations) are:


* any function of momenta of NOT strongly interacting particles

* sum of transverse quantities of all partons (e.g. ``H_T2``)

* any quantity referring to jets, constructed by an IR safe
  jet algorithm, see below.


Not infrared safe are

* any function of momenta of specific partons
* for processes with hadrons in the initial state: any quantity that depends on parton momenta along the beam axis, including the initial state partons itself.

Since the total number of partons is different for different pieces of
the NLO calculation any explicit reference to a parton momentum will
lead to an inconsistent result.

.. _Explicit scale variations:

Explicit scale variations
-------------------------

The (nominal) factorisation and renormalisation scales
in the fixed-order matrix elements can be scaled explicitly
simply by introducing a prefactor into the scale definition, e.g.

.. code-block:: yaml

   SCALES: VAR{0.25*H_T2}{0.25*H_T2}

for setting both the renormalisation and factorisation scales to
H_T/2.

However, to calculate several variations in a single event generation run,
you need to use :ref:`On-the-fly event weight variations`.
See the instructions given there
to find out how to vary factorisation and
renormalisation scale factors on-the-fly,
both in the matrix element and in the parton shower.

The starting scale of the parton shower resummation
in a ME+PS merged sample, ``MU_Q2``,
can at the moment not be varied on-the-fly.
To change the (nominal) starting scale explicitly,
a scale factor can be introduced
in the third argument of the METS scale setter:

.. code-block:: yaml

   SCALES: METS{MU_F2}{MU_R2}{4.0*MU_Q2}


.. _COUPLINGS:

COUPLINGS
=========

.. index:: COUPLINGS

Within Sherpa, strong and electroweak couplings can be computed at any scale
specified by a scale setter (cf. :ref:`SCALES`). The :option:`COUPLINGS` tag
links the argument of a running coupling to one of the respective scales.
This is better seen in an example. Assuming the following input

.. code-block:: yaml

   SCALES: VAR{...}{PPerp2(p[2])}{Abs2(p[2]+p[3])}
   COUPLINGS:
     - "Alpha_QCD 1"
     - "Alpha_QED 2"

Sherpa will compute any strong couplings at scale one,
i.e. ``PPerp2(p[2])`` and electroweak couplings at scale two,
i.e. ``Abs2(p[2]+p[3])``.  Note that counting starts at zero.

.. _KFACTOR:

KFACTOR
=======

.. index:: KFACTOR

This parameter specifies how to evaluate potential K-factors in the hard
process. This is equivalent to the :option:`COUPLINGS` specification of Sherpa
versions prior to 1.2.2. To list all available
K-factors, the tag ``SHOW_KFACTOR_SYNTAX: 1`` can be specified
on the command line. Currently available options are

:option:`None`
  No reweighting

:option:`VAR`
  Couplings specified by an additional parameter in a form which is understood
  by the internal interpreter, see :ref:`Interpreter`. The tags :kbd:`Alpha_QCD`
  and :kbd:`Alpha_QED` serve as links to the built-in running coupling implementation.

  If for example the process ``g g -> h g`` in effective theory is computed,
  one could think of evaluating two powers of the strong coupling at the Higgs mass scale
  and one power at the transverse momentum squared of the gluon.
  Assuming the Higgs mass to be 125 GeV, the corresponding reweighting would read

  .. code-block:: yaml

     SCALES:    VAR{...}{PPerp2(p[3])}
     COUPLINGS: "Alpha_QCD 1"
     KFACTOR:   VAR{sqr(Alpha_QCD(sqr(125))/Alpha_QCD(MU_12))}

  As can be seen from this example, scales are referred to as :kbd:`MU_<i>2`,
  where :kbd:`<i>` is replaced with the appropriate number.
  Note that counting starts at zero.

It is possible to implement a dedicated K-factor scheme within Sherpa.
For advice on this topic please contact the authors, :ref:`Authors`.

.. _YUKAWA_MASSES:

YUKAWA_MASSES
=============

.. index:: YUKAWA_MASSES

This parameter specifies whether the Yukawa couplings are evaluated
using running or fixed quark masses: ``YUKAWA_MASSES: Running`` is the
default since version 1.2.2 while ``YUKAWA_MASSES: Fixed`` was the
default until 1.2.1.

.. _Dipole subtraction:

Dipole subtraction
==================

.. index:: NLO_SUBTRACTION_MODE
.. index:: ALPHA
.. index:: ALPHA_FF
.. index:: ALPHA_FI
.. index:: ALPHA_IF
.. index:: ALPHA_II
.. index:: KAPPA
.. index:: NF_GSPLIT
.. index:: AMIN
.. index:: LIST

There is one general switch that governs the behaviour of the
Catani-Seymour subtraction :cite:`Catani1996vz` as implemented
in Sherpa :cite:`Gleisberg2007md,Schonherr2017qcj`.
`NLO_SUBTRACTION_MODE` defines which type of divergences will
be subtracted (both off the virtual and real amplitudes).
Options are:

:option:`QCD`
  Only QCD infrared divergences will be subtracted.
  This is the default as most users will be familiar with this
  setting corresponding to the abilities of older Sherpa
  versions.

:option:`QED`
  Only QED infrared divergences will be subtracted.

:option:`QCD+QED`
  All Standard Model infrared divergences will be subtracted.

Further, the following list of parameters can be used to optimise
the performance of the dipole subtraction.  These dipole
parameters are specified as subsettings to the ``DIPOLES`` setting,
like this:

.. code-block:: yaml

   DIPOLES:
     ALPHA: <alpha>
     NF_GSPLIT: <nf>
     # other dipole settings ...

The following parameters can be customised:

:option:`SCHEME`
  Defines the finite parts of the splitting functions used.
  Options are:

  `CS`,
  this is the standard Catani-Seymour subtraction definition
  of the splitting functions. This is the default for fixed-order
  calculations.

  `Dire`,
  this selects the modified splitting functions used in the Dire
  dipole shower. It is the default for MC@NLO calculations matching
  to Dire.

  `CSS`,
  this selects the modified splitting functions used in the CSS
  parton shower. It is the default for MC@NLO calculations matching
  to the CSS.

:option:`ALPHA`
  Specifies a dipole cutoff in the nonsingular region :cite:`Nagy2003tz`.
  Changing this parameter shifts contributions from the subtracted real
  correction piece (RS) to the piece including integrated dipole terms (I),
  while their sum remains constant. This parameter can be used to optimize
  the integration performance of the individual pieces.
  Also the average calculation time for the subtracted real correction
  is reduced with smaller choices of "ALPHA" due to the (on average)
  reduced number of contributing dipole terms. For most processes
  a reasonable choice is between 0.01 and 1 (default). See also
  :ref:`Choosing DIPOLES ALPHA`

:option:`ALPHA_FF, ALPHA_FI, ALPHA_IF, ALPHA_II`
  Specifies the above dipole alpha only for FF, FI, IF, or II dipoles, respectively.

:option:`AMIN`
  Specifies the cutoff of real correction terms in the infrared region
  to avoid numerical problems with the subtraction. The default is 1.e-8.

:option:`NF_GSPLIT`
  Specifies the number of quark flavours that are produced from
  gluon splittings. This number must be at least the number of massless
  flavours (default). If this number is larger than the number of massless
  quarks the massive dipole subtraction :cite:`Catani2002hc` is employed.

:option:`KAPPA`
  Specifies the kappa-parameter in the massive dipole subtraction formalism
  :cite:`Catani2002hc`. The default is 2.0/3.0.

:option:`LIST`
  If set to 1 all generated dipoles will be listed for all generated processes.
