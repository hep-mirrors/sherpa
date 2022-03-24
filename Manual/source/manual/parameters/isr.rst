.. _ISR Parameters:

**************
ISR parameters
**************

.. index:: BUNCH_1
.. index:: BUNCH_2
.. index:: ISR_SMIN
.. index:: ISR_SMAX
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

:OPTION:`ISR_SMIN/ISR_SMAX`
  This parameter specifies the minimum fraction of cms energy squared
  after ISR. The reference value is the total centre of mass energy
  squared of the collision, `not` the centre of mass energy after
  eventual Beamstrahlung.
  The parameter can be specified using the internal interpreter, see
  :ref:`Interpreter`, e.g. as ``ISR_SMIN: sqr(20/E_CMS)``.

Sherpa provides access to a variety of structure functions.
They can be configured with the following parameters.

:OPTION:`PDF_LIBRARY`
  This parameter takes the list of PDF interfaces to load.  The
  following options are distributed with Sherpa:

  :option:`LHAPDFSherpa`
    Use PDF's from LHAPDF :cite:`Whalley2005nh`. The interface is only
    available if Sherpa has been compiled with support for LHAPDF, see
    :ref:`Installation`.

  :option:`CT14Sherpa`
    Built-in library for some PDF sets from the CTEQ collaboration,
    cf. :cite:`Dulat2015mca`. This is the default.

  :option:`CT12Sherpa`
    Built-in library for some PDF sets from the CTEQ collaboration,
    cf. :cite:`Gao2013xoa`.

  :option:`CT10Sherpa`
    Built-in library for some PDF sets from the CTEQ collaboration,
    cf. :cite:`Lai2010vv`.

  :option:`CTEQ6Sherpa`
    Built-in library for some PDF sets from the CTEQ collaboration,
    cf. :cite:`Nadolsky2008zw`.

  :OPTION:`NNPDF30Sherpa`
    Built-in library for PDF sets from the NNPDF group, cf. :cite:`Ball2014uwa`.

  :option:`MSTW08Sherpa`
    Built-in library for PDF sets from the MSTW group, cf. :cite:`Martin2009iq`.

  :option:`MRST04QEDSherpa`
    Built-in library for photon PDF sets from the MRST group, cf. :cite:`Martin2004dh`.

  :option:`MRST01LOSherpa`
    Built-in library for the 2001 leading-order PDF set from the MRST group, cf. :cite:`Martin2001es`.

  :option:`MRST99Sherpa`
    Built-in library for the 1999 PDF sets from the MRST group, cf. :cite:`Martin1999ww`.

  :option:`GRVSherpa`
    Built-in library for the GRV photon PDF :cite:`Gluck1991jc`, :cite:`Gluck1991ee`.

  :option:`GRSSherpa`
    Built-in library for the GRS photon PDF :cite:`Gluck1999ub`.

  :option:`SALSherpa`
    Built-in library for the SAL photon PDF :cite:`Slominski2005bw`.

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
