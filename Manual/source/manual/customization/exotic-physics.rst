.. _Exotic physics:

**************
Exotic physics
**************


It is possible to add your own models to Sherpa in a straightforward way.
To illustrate, a simple example has been included in the directory
``Examples/Models/SM_ZPrime``, showing how to add a Z-prime boson to
the Standard Model.

The important features of this example include:

* The ``SM_Zprime.C`` file.

  This file contains the initialisation of the Z-prime boson. The
  properties of the Z-prime are set here, such as mass, width,
  electromagnetic charge, spin etc.

* The ``Interaction_Model_SM_Zprime.C`` file.

  This file contains the definition of the Z-prime boson's interactions.
  The right- and left-handed couplings to each of the fermions are
  set here.

* An example ``Makefile``.

  This shows how to compile the sources above into a shared library.

* The line ``SHERPA_LDADD: SMZprime`` in the config file.

  This line tells Sherpa to load the extra libraries created from the
  `*.C` files above.

* The line ``MODEL: SMZprime`` in the config file.

  This line tells Sherpa which model to use for the run.

* The following lines in the config file:

  .. code-block:: yaml

     PARTICLE_DATA:
       32:
         Mass: 1000
         Width: 50

  These lines show how you can overrule the choices you made for the
  properties of the new particle in the :file:`SM_Zprime.C` file. For
  more information on changing parameters in Sherpa, see :ref:`Input
  structure` and :ref:`Parameters`.

* The lines ``Zp_cpl_L: 0.3`` and ``Zp_cpl_R: 0.6`` set the couplings
  to left and right handed fermions in the config file.

To use this model, create the libraries for Sherpa to use by running

.. code-block:: shell-session

   $ make

in this directory. Then run Sherpa as normal:

.. code-block:: shell-session

   $ ../../../bin/Sherpa

To implement your own model, copy these example files anywhere and
modify them according to your needs.

Note: You don't have to modify or recompile any part of Sherpa to use
your model. As long as the ``SHERPA_LDADD`` parameter is specified as
above, Sherpa will pick up your model automatically.

Furthermore note: New physics models with an existing implementation
in FeynRules, cf. :cite:`Christensen2008py` and
:cite:`Christensen2009jx`, can directly be invoked using Sherpa's
support for the UFO model format, see :ref:`UFO Model Interface`.
