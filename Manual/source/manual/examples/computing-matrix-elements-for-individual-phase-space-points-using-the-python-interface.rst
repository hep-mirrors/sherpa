.. _Matrix Element values through Python interface:

Computing matrix elements for individual phase space points using the Python Interface
======================================================================================

Sherpa's Python interface (see :ref:`Python Interface`) can be used to
compute matrix elements for individual phase space points.  Access to
a designated class "MEProcess" is provided by interface to compute
matrix elements as illustrated in the example script.

Please note that the process in the script must be compatible with the
one specified in the Sherpa configuration file in the working
directory. A random phase space point for the process of interest can
be generated as shown in the example.

If AMEGIC++ is used as the matrix element generator, executing the
script will result in AMEGIC++ writing out libraries and exiting.
After compiling the libraries using ``./makelibs``, the script must be
executed again in order to obtain the matrix element.

.. literalinclude:: /../../Examples/API/ME2-Python/test.py.in
   :language: text

.. _Matrix Element values through C++ interface:

Computing matrix elements for individual phase space points using the C++ Interface
===================================================================================

.. index:: MOMENTA_DATA_FILE

Matrix elements values for user defined phase space points can also be
quarried using a small C++ executable provided in
``Examples/API/ME2``.  It can be compiled using the provided
``Makefile``. The test program is then run typing (note: the
``LD_LIBRARY_PATH`` must be set to include
``<Sherpa-installation>/lib/SHERPA-MC``)

.. code-block:: shell-session

   $ ./test <options>

where the usual options for Sherpa are passed. An example
configuration file, giving both the process and the requested phase
space points looks like

.. literalinclude:: /../../Examples/API/ME2-CPP/Sherpa.yaml
   :language: yaml

Please note that both the process and the beam specifications need to
be present in order for Sherpa to initialise properly. The momenta
need to be given in the proper ordering employed in Sherpa, which can
be read from the process name printed on screen. For each entry the
sequence is the following


.. code-block:: text

   [<pdg-id>, <E>, <px>, <py>, <pz>, triplet-index, antitriplet-index]

with the colour indices ranging from 1..3 for both the triplet and the
antitriplet index in the colour-flow basis. The colour information is only
needed if Comix is used for the calculation as Comix then also gives the
squared matrix element value for this colour configuration. Otherwise, the last
two arguments can be omitted. In any case, the colour-summed value is printed
to screen.
