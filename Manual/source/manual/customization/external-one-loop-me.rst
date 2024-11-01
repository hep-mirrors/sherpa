.. _External one-loop ME:

********************
External one-loop ME
********************

.. index:: LHOLE_ORDERFILE
.. index:: LHOLE_CONTRACTFILE
.. index:: LHOLE_IR_REGULARISATION
.. index:: LHOLE_BOOST_TO_CMS
.. index:: LHOLE_OLP

Sherpa includes only a very limited selection of one-loop matrix
elements. To make full use of the implemented automated dipole
subtraction it is possible to link external one-loop codes to Sherpa
in order to perform full calculations at QCD next-to-leading order.

In general Sherpa can take care of any piece of the calculation except
one-loop matrix elements, i.e. the born ME, the real correction, the
real and integrated subtraction terms as well as the phase space
integration and PDF weights for hadron collisions. Sherpa will provide
sets of four-momenta and request for a specific parton level process
the helicity and colour summed one-loop matrix element (more specific:
the coefficients of the Laurent series in the dimensional
regularisation parameter epsilon up to the order epsilon^0).

An example setup for interfacing such an external one-loop code,
following the Binoth Les Houches interface proposal
:cite:`Binoth2010xt` of the 2009 Les Houches workshop, is provided in
:ref:`LHC_Zbb`. To use the LH-OLE interface, Sherpa has to be
configured with :option:`-DSHERPA_ENABLE_LHOLE=ON`.

The interface:

* During an initialization run Sherpa stores setup information
  (schemes, model information etc.) and requests a list of parton-level
  one-loop processes that are needed for the NLO calculation. This information
  is stored in a file, by default called ``OLE_order.lh``.
  The external one-loop code (OLE) should confirm these settings/requests
  and write out a file ``OLE_contract.lh``. Both filenames can be customised
  using ``LHOLE_ORDERFILE: <order-file>`` and
  ``<LHOLE_CONTRACTFILE: <contract-file>``. For the syntax of these files and
  more details see :cite:`Binoth2010xt`.

  For Sherpa the output/input of the order/contract file is handled
  in ``LH_OLE_Communicator.[CH]``.
  The actual interface is contained in ``LH_OLE_Interface.C``.
  The parameters to be exchanged with the OLE are defined in the
  latter file via

  .. code-block:: c++

     lhfile.AddParameter(...);

  and might require an update for specific OLE or processes. Per
  default, in addition to the standard options
  ``MatrixElementSquareType``, ``CorrectionType``,
  ``IRregularisation``, ``AlphasPower``, ``AlphaPower`` and
  ``OperationMode`` the masses and width of the W, Z and Higgs bosons
  and the top and bottom quarks are written out in free format, such
  that the respective OLE parameters can be easily synchronised.

* At runtime the communication is performed via function calls.
  To allow Sherpa to call the external code the functions

  .. code-block:: c++

     void OLP_Start(const char * filename);
     void OLP_EvalSubProcess(int,double*,double,double,double*);

  which are defined and called in ``LH_OLE_Interface.C`` must be
  specified.  For keywords and possible data fields passed with this
  functions see :cite:`Binoth2010xt`.

  The function ``OLP_Start(...)`` is called once when Sherpa is
  starting.  The function ``OLP_EvalSubProcess(...)`` will be called
  many times for different subprocesses and momentum configurations.

The setup (cf. example :ref:`LHC_Zbb`):

* The line ``Loop_Generator: LHOLE`` tells the code to use
  the interface for computing one-loop matrix elements.

* The switch ``SHERPA_LDADD`` has to be set to the appropriate
  library name (and path) of the one-loop generator.

* The IR regularisation scheme can be set via
  ``LHOLE_IR_REGULARISATION``. Possible values are ``DRED`` (default)
  and ``CDR``.

* Per default, Sherpa generates phase space points in the lab frame.
  If ``LHOLE_BOOST_TO_CMS: true`` is set, these phase space points are
  boosted to the centre of mass system before they are passed to the
  OLE.

* The original BLHA interface does not allow for run-time parameter
  passing. While this is discussed for an updated of the accord a
  workable solution is implemented for the use of GoSam and enabled
  through ``LHOLE_OLP: GoSam``. The ``LHOLE_BOOST_TO_CMS`` is also
  automatically active with this setup. This, of course, can be
  adapted for other one-loop programs if need be.

* Sherpa's internal analysis package can be used to generate a few
  histograms. Thus, then when installing Sherpa the option
  :option:`-DSHERPA_ENABLE_ANALYSIS=ON` must be include on the command line when
  Sherpa is configured, see :ref:`ANALYSIS`. However, this module is deprecated
  and will not be supported in the future.
