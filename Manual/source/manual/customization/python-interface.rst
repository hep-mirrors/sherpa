.. _Python Interface:

****************
Python Interface
****************

Certain Sherpa classes and methods can be made available to the Python
interpreter in the form of an extension module. This module can be
loaded in Python and provides access to certain functionalities of the
Sherpa event generator in Python. In order to build the module, Sherpa
must be configured with the option :option:`-DSHERPA_ENABLE_PYTHON=ON`.
Running ``make`` then invokes the automated interface generator SWIG
:cite:`Beazley2003` to create the Sherpa module using the Python C/C++
API. SWIG version 1.3.x or later is required for a successful build.
Problems might occur if more than one version of Python is present on
the system since automake currently doesn't always handle multiple
Python installations properly. If you have multiple Python versions
installed on your system, please set the ``PYTHON`` environment
variable to the Python 3 executable via

.. code-block:: shell-session

   $ export PYTHON=<path-to-python3>

before executing ``cmake`` script (see.  Certain Sherpa
classes and methods can be made available to the Python interpreter in
the form of an extension module. This module can be loaded in Python
and provides access to certain functionalities of the Sherpa event
generator in Python. It was designed specifically for the computation
of matrix elements in python (:ref:`APIexamples`) and its features are
currently limited to this purpose. In order to build the module,
Sherpa must be configured with the option
:option:`-DSHERPA_ENABLE_PYTHON=ON`.
Running ``make`` then invokes the automated interface generator SWIG
:cite:`Beazley2003` to create the Sherpa module using the Python C/C++
API. SWIG version 1.3.x or later is required for a successful build.
Problems might occur if more than one version of Python is present on
the system since automake currently doesn't always handle multiple
Python installations properly. A possible workaround is to temporarily
uninstall one version of python, configure and build Sherpa, and then
reinstall the temporarily uninstalled version of Python.

The following script is a minimal example that shows how to use the
Sherpa module in Python. In order to load the Sherpa module, the
location where it is installed must be added to the PYTHONPATH. There
are several ways to do this, in this example the sys module is
used. The sys module also allows it to directly pass the command line
arguments used to run the script to the initialization routine of
Sherpa. The script can thus be executed using the normal command line
options of Sherpa (see :ref:`Command line`). Furthermore it
illustrates how exceptions that Sherpa might throw can be taken care
of. If a run card is present in the directory where the script is
executed, the initialization of the generator causes Sherpa to compute
the cross sections for the processes specified in the run card. See
:ref:`Matrix Element values through Python interface` for an example
that shows how to use the Python interface to compute matrix elements
or :ref:`Events` to see how the interface can be used to generate
events in Python.

Note that if you have compiled Sherpa with MPI support, you need to
source the `mpi4py <http://mpi4py.scipy.org>`_ module using ``from
mpi4py import MPI``.


.. code-block:: python

   #!/usr/bin/python
   import sys
   sys.path.append('<sherpa-prefix>/lib/<your-python-version>/site-packages/>')
   import Sherpa

   # set up the generator
   Generator=Sherpa.Sherpa(len(sys.argv),sys.argv)

   # initialize the generator, pass command line arguments to initialization routine
   try:
     Generator.InitializeTheRun()

    # catch exceptions
    except Sherpa.Exception as exc:
      print exc
