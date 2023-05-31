.. _ISR Parameters:

**************
ISR parameters
**************

.. index:: BUNCH_1
.. index:: BUNCH_2
.. index:: ISR_E_ORDER
.. index:: ISR_E_SCHEME
.. index:: PDF_LIBRARY
.. index:: PDF_LIBRARY_1
.. index:: PDF_LIBRARY_2
.. index:: PDF_SET
.. index:: PDF_SET_1
.. index:: PDF_SET_2
.. index:: PDF_SET_VERSION
.. index:: PDF_SET_VERSION_1
.. index:: PDF_SET_VERSION_2
.. index:: SHOW_PDF_SETS

The following parameters are used to steer the setup of beam substructure and
initial state radiation (ISR).

:OPTION:`BUNCHES`
  Specify the PDG ID of the first (left) and second (right) bunch
  particle (or both if only one value is provided), i.e. the particle
  after eventual Beamstrahlung specified through the beam parameters,
  see :ref:`Beam Parameters`.  Per default these are taken to be
  identical to the values set using :OPTION:`BEAMS`, assuming the default
  beam spectrum is Monochromatic. In case the Simple Compton or Laser
  Backscattering spectra are enabled the bunch particles would have to
  be set to 22, the PDG code of the photon.

Sherpa provides access to a variety of structure functions.
They can be configured with the following parameters.

:OPTION:`PDF_LIBRARY`
  This parameter takes the list of PDF interfaces to load.  The
  following options are distributed with Sherpa:

  :option:`LHAPDFSherpa`
    Use PDF's from LHAPDF :cite:`Buckley2011ms`. This is the default.

  :option:`CT14Sherpa`
    Built-in library for some PDF sets from the CTEQ collaboration,
    cf. :cite:`Dulat2015mca`.

  :OPTION:`NNPDFSherpa`
    Built-in library for PDF sets from the NNPDF group, cf. :cite:`Ball2014uwa`.

  :option:`GRVSherpa`
    Built-in library for the GRV photon PDF :cite:`Gluck1991jc`, :cite:`Gluck1991ee`.

  :option:`GRSSherpa`
    Built-in library for the GRS photon PDF :cite:`Gluck1999ub`.

  :option:`SALSherpa`
    Built-in library for the SAL photon PDF :cite:`Slominski2005bw`.

  :option:`CJKSherpa`
    Built-in library for the CJK photon PDF :cite:`Cornet:2002iy`,
    :cite:`Cornet:2003ry`, :cite:`Cornet:2004ng`, :cite:`Cornet:2004nb`.

  :option:`SASGSherpa`
    Built-in library for the SaSgam photon PDF :cite:`Schuler1995fk`, :cite:`Schuler1996fc`.

  :option:`PDFESherpa`
    Built-in library for the electron structure function.  The
    perturbative order of the fine structure constant can be set using
    the parameter :OPTION:`ISR_E_ORDER` (default: 1). The switch
    :OPTION:`ISR_E_SCHEME` allows to set the scheme of respecting non-leading
    terms. Possible options are 0 ("mixed choice"), 1 ("eta choice"), or
    2 ("beta choice", default).

  :option:`None`
    No PDF. Fixed beam energy.

  Furthermore it is simple to build an external interface to an
  arbitrary PDF and load that dynamically in the Sherpa run. See
  :ref:`External PDF` for instructions.

  By default, Sherpa will try to install with the LHAPDF interface enabled.
  If this is not desired, for example in lepton-lepton collisions
  where LHAPDF is not used, the user can disable the interface
  with the cmake option ``-DCMAKE_ENABLE_LHAPDF=OFF``.
  Sherpa will then use the internal :OPTION:`PDF_LIBRARY` for
  hadronic collisions, with the default set being
  ``NNPDF31_nnlo_as_0118_mc``. Note that PDF variations
  and the evolution of ALPHAS, ``ALPHAS: {USE_PDF: 1}``,
  can only be used with LHAPDF enabled.

:OPTION:`PDF_SET`
  Specifies the PDF set for hadronic bunch particles. All
  sets available in the chosen :OPTION:`PDF_LIBRARY` can be figured by
  running Sherpa with the parameter :OPTION:`SHOW_PDF_SETS: 1`, e.g.:

  .. code-block:: shell-session

     $ Sherpa 'PDF_LIBRARY: CTEQ6Sherpa' 'SHOW_PDF_SETS: 1'

  If the two colliding beams are of different type, e.g. protons and
  electrons or photons and electrons, it is possible to specify two
  different PDF sets by providing two values: :option:`PDF_SET: [pdf1,
  pdf2]`. The special value ``Default`` can be used as a placeholder
  for letting Sherpa choose the appropriate PDF set (or none).

:OPTION:`PDF_SET_VERSIONS`
  This parameter allows to select a specific
  version (member) within the chosen PDF set. It is possible to
  specify two different PDF sets using :option:`PDF_SET_VERSIONS:
  [version1, version2]`

See :ref:`On-the-fly event weight variations`
to find out how to vary PDF sets and version on-the-fly,
both in the matrix element and in the parton shower.
