.. _MPIEvents:

Generate events with MPI using scripts
======================================

This example shows how to generate events with Sherpa using a Python
wrapper script and MPI. For each event the weight, the number of
trials and the particle information is send to the MPI master node and
written into a single gzip'ed output file.  Note that you need the
`mpi4py <http://mpi4py.scipy.org>`_ module to run this Example. Sherpa
must be configured and installed using :option:`-DSHERPA_ENABLE_MPI`, see
:ref:`MPI parallelization`.

.. literalinclude:: /../../Examples/API/MPIEvents/test.py.in
   :language: python
