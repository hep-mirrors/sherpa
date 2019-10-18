.. _Calculation of inclusive cross sections:

Calculation of inclusive cross sections
=======================================

Note that this example is not yet updated to the new YAML input
format.  Contact the :ref:`Authors` for more information.

.. literalinclude:: /../../Examples/Soft_QCD/Run.Xsecs.dat
   :language: console

Things to notice:

* Inclusive cross sections (total, inelastic, low-mass
  single-diffractive, low-mass double-diffractive, elastic) and the
  elastic slope are calculated for varying centre-of-mass energies in
  pp collisions

* The results are written to the file
  InclusiveQuantities/xsecs_total.dat and to the screen.  The
  directory will automatically be created in the path from where
  Sherpa is run.

* The parameters of the model are not very well tuned.
