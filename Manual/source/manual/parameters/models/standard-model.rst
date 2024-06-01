.. _SM:

Standard Model
--------------

.. index:: 1/ALPHAQED(0)
.. index:: 1/ALPHAQED(MW)
.. index:: 1/ALPHAQED(MZ)
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
.. index:: GMU_CMS_AQED_CONVENTION
.. index:: ORDER_ALPHAS
.. index:: ALPHAS_USE_PDF
.. index:: WIDTH_SCHEME
.. index:: PARTICLE_DATA_Mass
.. index:: PARTICLE_DATA_Massive
.. index:: PARTICLE_DATA_Width
.. index:: PARTICLE_DATA_Active
.. index:: PARTICLE_DATA_Stable
.. index:: PARTICLE_DATA_Yukawa

The SM inputs for the electroweak sector can be given in nine
different schemes, that correspond to different choices of which SM
physics parameters are considered fixed and which are derived from the
given quantities. The electroweak coupling is by default fixed, unless its
running has been enabled (cf. :ref:`COUPLINGS`).
The input schemes are selected through the ``EW_SCHEME`` parameter, whose
default is :option:`Gmu`. The following options are provided:

:option:`UserDefined`
  All EW parameters are explicitly given:  Here the W, Z and Higgs
  masses and widths are taken as inputs, and the parameters
  ``1/ALPHAQED(0)``, ``ALPHAQED_DEFAULT_SCALE``, ``SIN2THETAW`` (weak mixing
  angle), ``VEV`` (Higgs field vacuum expectation value) and
  ``LAMBDA`` (Higgs quartic coupling) have to be specified.

  By default, :option:`ALPHAQED_DEFAULT_SCALE: 8315.18` (:math:`=m_Z^2`), which
  means that the MEs are evaluated with a value of :math:`\alpha=\frac{1}{128.802}`.

  Note that this mode allows to violate the tree-level relations
  between some of the parameters and might thus lead to gauge
  violations in some regions of phase space.

:option:`alpha0`
  All EW parameters are calculated from the W, Z and Higgs masses
  and widths and
  the fine structure constant (taken from ``1/ALPHAQED(0)`` +
  ``ALPHAQED_DEFAULT_SCALE``, cf. below) using tree-level relations.

  By default, :option:`ALPHAQED_DEFAULT_SCALE: 0.0`, which means that the MEs
  are evaluated with a value of :math:`\alpha=\frac{1}{137.03599976}`.

:option:`alphamZ`
  All EW parameters are calculated from the W, Z and Higgs masses
  and widths and the fine structure constant (taken from ``1/ALPHAQED(MZ)``,
  default :option:`128.802`) using tree-level relations.

:option:`Gmu`
  This choice corresponds to the G_mu-scheme. The EW parameters are
  calculated out of the weak gauge boson masses M_W, M_Z, the Higgs
  boson mass M_H, their respective widths,
  and the Fermi constant ``GF`` using tree-level
  relations.

:option:`alphamZsW`
  All EW parameters are calculated from the Z and Higgs masses and
  widths, the fine structure constant (taken from ``1/ALPHAQED(MZ)``,
  default :option:`128.802`),
  and the weak mixing angle (``SIN2THETAW``) using
  tree-level relations. In particular, the W boson mass (and in the
  complex mass scheme also its width) is a derived quantity.

:option:`alphamWsW`
  All EW parameters are calculated from the W and Higgs masses and
  widths, the fine structure constant (taken from ``1/ALPHAQED(MW)``,
  default :option:`132.17`),
  and the weak mixing angle (``SIN2THETAW``) using
  tree-level relations. In particular, the Z boson mass (and in the
  complex mass scheme also its width) is a derived quantity.

:option:`GmumZsW`
  All EW parameters are calculated from the Z and Higgs masses and
  widths, the Fermi constant (``GF``),
  and the weak mixing angle (``SIN2THETAW``) using
  tree-level relations. In particular, the W boson mass (and in the
  complex mass scheme also its width) is a derived quantity.

:option:`GmumWsW`
  All EW parameters are calculated from the W and Higgs masses and
  widths, the Fermi constant (``GF``),
  and the weak mixing angle (``SIN2THETAW``) using
  tree-level relations. In particular, the Z boson mass (and in the
  complex mass scheme also its width) is a derived quantity.

:option:`FeynRules`
  This choice corresponds to the scheme employed in the FeynRules/UFO
  setup.  The EW parameters are calculated out of the Z boson mass
  M_Z, the Higgs boson mass M_H, the Fermi constant ``GF`` and the
  fine structure constant (taken from ``1/ALPHAQED(0)`` +
  ``ALPHAQED_DEFAULT_SCALE``, cf. below) using tree-level
  relations. Note, the W boson mass is not an input parameter in this
  scheme.

All ``Gmu``-derived schemes, where the EW coupling is a derived quantity,
possess an ambiguity on how to construct a real EW coupling in the
complex mass scheme. Several conventions are implemented and can
be accessed through ``GMU_CMS_AQED_CONVENTION``.

To account for quark mixing the CKM matrix elements have to be
assigned.  For this purpose the Wolfenstein parameterisation
:cite:`Wolfenstein1983yz` is employed. The order of expansion in the
lambda parameter is defined through

.. code-block:: yaml

   CKM:
     Order: <order>
     # other CKM settings ...

The default for ``Order`` is :option:`0`, corresponding to a unit
matrix.  The parameter convention for higher expansion terms reads:


* ``Order: 1``, the ``Cabibbo`` subsetting has to be set, it
  parameterises lambda and has the default value :option:`0.22537`.

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
Wolfenstein parameterisation. Setting ``CKM: {Output: true}`` enables
an output of the CKM matrix.

The remaining parameter to fully specify the Standard Model is the
strong coupling constant at the Z-pole, given through
``ALPHAS(MZ)``. Its default value is :option:`0.118`. If the setup at
hand involves hadron collisions and thus PDFs, the value of the strong
coupling constant is automatically set consistent with the PDF fit and
can not be changed by the user. Since Sherpa is compiled with LHAPDF
support, it is also possible to use the alphaS evolution provided in
LHAPDF by specifying ``ALPHAS: {USE_PDF: 1}``. The perturbative
order of the running of the strong coupling can be set via
``ORDER_ALPHAS``, where the default :option:`0` corresponds to
one-loop running and :option:`1`, :option:`2`, :option:`3` to :option:`2,3,4`-loops,
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
