.. _Calculation of inclusive cross sections:

Calculation of inclusive cross sections
=======================================

Note that this example is not yet updated to the new YAML input
format.  Contact the :ref:`Authors` for more information.

.. code-block:: console

    (run){
      OUTPUT              = 2
      EVENT_TYPE          = MinimumBias
      SOFT_COLLISIONS     = Shrimps
      Shrimps_Mode        = Xsecs

      deltaY    =  1.5;
      Lambda2   =  1.7;
      beta_0^2  =  20.0;
      kappa     =  0.6;
      xi        =  0.2;
      lambda    =  0.3;
      Delta     =  0.4;
    }(run)

    (beam){
      BEAM_1 =  2212; BEAM_ENERGY_1 = 450.;
      BEAM_2 =  2212; BEAM_ENERGY_2 = 450.;
    }(beam)

    (me){
      ME_SIGNAL_GENERATOR = None
    }(me)

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
