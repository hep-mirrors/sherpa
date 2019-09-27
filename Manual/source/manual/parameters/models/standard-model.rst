.. _SM:

Standard Model
--------------

.. index:: 1/ALPHAQED(0)
.. index:: SIN2THETAW
.. index:: VEV
.. index:: CKM_Lambda
.. index:: EW_SCHEME
.. index:: CKM_Order
.. index:: CKM_Cabibbo
.. index:: CKM_A
.. index:: CKM_Eta
.. index:: CKM_Rho
.. index:: CKM_Matrix_Elements
.. index:: CKM_Output
.. index:: ALPHAS(MZ)
.. index:: ALPHAQED_DEFAULT_SCALE
.. index:: ORDER_ALPHAS
.. index:: ALPHAS_USE_PDF
.. index:: WIDTH_SCHEME
.. index:: PARTICLE_DATA_Mass
.. index:: PARTICLE_DATA_Massive
.. index:: PARTICLE_DATA_Width
.. index:: PARTICLE_DATA_Active
.. index:: PARTICLE_DATA_Stable
.. index:: PARTICLE_DATA_Yukawa

The SM inputs for the electroweak sector can be given in four
different schemes, that correspond to different choices of which SM
physics parameters are considered fixed and which are derived from the
given quantities. The input schemes are selected through the
``EW_SCHEME`` parameter, whose default is :option:`1`. The following
options are provided:

* Case ``0``:

  all EW parameters are explicitly given.  Here the W, Z and Higgs
  masses are taken as inputs, and the parameters ``1/ALPHAQED(0)``,
  ``ALPHAQED_DEFAULT_SCALE`` (cf. below), ``SIN2THETAW`` (weak mixing
  angle), ``VEV`` (Higgs field vacuum expectation value) and
  ``LAMBDA`` (Higgs quartic coupling) have to be specified.

  Note that this mode allows to violate the tree-level relations
  between some of the parameters and might thus lead to gauge
  violations in some regions of phase space.

* Case ``1``:

  all EW parameters are calculated from the W, Z and Higgs masses and
  the fine structure constant (taken from ``1/ALPHAQED(0)`` +
  ``ALPHAQED_DEFAULT_SCALE``, cf. below) using tree-level relations.

* Case ``2``:

  all EW parameters are calculated from the W, Z and Higgs masses and
  the fine structure constant (taken from ``1/ALPHAQED(MZ)``) using
  tree-level relations.

* Case ``3``:

  this choice corresponds to the G_mu-scheme. The EW parameters are
  calculated out of the weak gauge boson masses M_W, M_Z, the Higgs
  boson mass M_H and the Fermi constant ``GF`` using tree-level
  relations.

* Case ``4``:

  this choice corresponds to the scheme employed in the FeynRules/UFO
  setup.  The EW parameters are calculated out of the Z boson mass
  M_Z, the Higgs boson mass M_H, the Fermi constant ``GF`` and the
  fine structure constant (taken from ``1/ALPHAQED(0)`` +
  ``ALPHAQED_DEFAULT_SCALE``, cf. below) using tree-level
  relations. Note, the W boson mass is not an input parameter in this
  scheme.


The electro-weak coupling is by default not running, unless its
running has been enabled (cf. :ref:`COUPLINGS`). In EW schemes 0 and
1, the squared scale at which the fixed EW coupling is to be evaluated
can be specified by ``ALPHAQED_DEFAULT_SCALE``, which defaults to the
Z mass squared (note that this scale has to be specified in GeV
squared).  By default ``1/ALPHAQED(0): 137.03599976`` and
``ALPHAQED_DEFAULT_SCALE: 8315.18`` (:math:`=mZ^2`), which means that the MEs
are evaluated with a fixed value of ``alphaQED=1/128.802``.

To account for quark mixing the CKM matrix elements have to be
assigned.  For this purpose the Wolfenstein parametrization
:cite:`Wolfenstein1983yz` is employed. The order of expansion in the
lambda parameter is defined through

.. code-block:: yaml

   CKM:
     Order: <order>
     # other CKM settings ...

The default for ``Order`` is :option:`0`, corresponding to a unit
matrix.  The parameter convention for higher expansion terms reads:


* ``Order: 1``, the ``Cabibbo`` subsetting has to be set, it
  parametrizes lambda and has the default value :option:`0.22537`.

* ``Order: 2``, in addition the value of ``CKM_A`` has to be set, its
  default is :option:`0.814`.

* ``Order: 3``, the order lambda^3 expansion, ``Eta`` and ``Rho`` have
  to be specified. Their default values are :option:`0.353` and
  :option:`0.117`, respectively.


The CKM matrix elements V_ij can also be read in using

.. code-block:: yaml

   CKM:
     Matrix_Elements:
       i,j: <V_ij>
       # other CKM matrix elements ...
     # other CKM settings ...

Complex values can be given by providing two values: ``<V_ij> -> [Re,
Im]``.  Values not explicitly given are taken from the afore computed
Wolfenstein parametrisation. Setting ``CKM: @{Output: true@``} enables
an output of the CKM matrix.

The remaining parameter to fully specify the Standard Model is the
strong coupling constant at the Z-pole, given through
``ALPHAS(MZ)``. Its default value is :option:`0.118`. If the setup at
hand involves hadron collisions and thus PDFs, the value of the strong
coupling constant is automatically set consistent with the PDF fit and
can not be changed by the user. If Sherpa is compiled with LHAPDF
support, it is also possible to use the alphaS evolution provided in
LHAPDF by specifying ``ALPHAS: @{USE_PDF: 1@``}. The perturbative
order of the running of the strong coupling can be set via
``ORDER_ALPHAS``, where the default :option:`0` corresponds to
one-loop running and ``1``, ``2``, ``3`` to ``2,3,4``-loops,
respectively. If the setup at hand involves PDFs, this parameter is
set consistent with the information provided by the PDF set.

If unstable particles (e.g. W/Z bosons) appear as intermediate
propagators in the process, Sherpa uses the complex mass scheme to
construct MEs in a gauge-invariant way. For full consistency with this
scheme, by default the dependent EW parameters are also calculated
from the complex masses (:option:`WIDTH_SCHEME: CMS`), yielding
complex values e.g. for the weak mixing angle.  To keep the parameters
real one can set :option:`WIDTH_SCHEME: Fixed`. This may spoil gauge
invariance though.

With the following switches it is possible to change the properties of
all fundamental particles:

.. code-block:: yaml

   PARTICLE_DATA:
     <id>:
       <Property>: <value>
       # other properties for this particle ...
     # data for other particles

Here, ``<id>`` is the PDG ID of the particle for which one more
properties are to be modified. ``<Property>`` can be one of the
following:

:option:`Mass`
  Sets the mass (in GeV) of the particle.

  Masses of particles and corresponding anti-particles are always set
  simultaneously.

  For particles with Yukawa couplings, those are enabled/disabled
  consistent with the mass (taking into account the :option:`Massive`
  parameter) by default, but that can be modified using the
  :option:`Yukawa` parameter. Note that by default the Yukawa
  couplings are treated as running, cf. :ref:`YUKAWA_MASSES`.

:option:`Massive`
  Specifies whether the finite mass of the particle is to be considered
  in matrix-element calculations or not. Can be :option:`true` or
  :option:`false`.

:option:`Width`
  Sets the width (in GeV) of the particle.

:option:`Active`
  Enables/disables the particle with PDG id :option:`<id>`. Can be
  :option:`true` or :option:`false`.

:option:`Stable`
  Sets the particle either stable or unstable according
  to the following options:

  :option:`0`
    Particle and anti-particle are unstable

  :option:`1`
    Particle and anti-particle are stable

  :option:`2`
    Particle is stable, anti-particle is unstable

  :option:`3`
    Particle is unstable, anti-particle is stable



  This option applies to decays of hadrons (cf. :ref:`Hadron decays`)
  as well as particles produced in the hard scattering (cf. :ref:`Hard
  decays`).  For the latter, alternatively the decays can be specified
  explicitly in the process setup (see :ref:`Processes`) to avoid the
  narrow-width approximation.

:option:`Priority`
  Allows to overwrite the default automatic flavour
  sorting in a process by specifying a priority for the given
  flavour. This way one can identify certain particles which are part
  of a container (e.g. massless b-quarks), such that their position
  can be used reliably in selectors and scale setters.


.. note::

   :OPTION:`PARTICLE_DATA` can also be used to the properties of hadrons,
   you can use the same switches (except for :option:`Massive`), see
   :ref:`Hadronization`.
