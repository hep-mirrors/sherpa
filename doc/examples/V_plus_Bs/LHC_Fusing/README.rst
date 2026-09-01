.. _LHC_Fusing:

Zbb with fusing
===============

The fusing algorithm as introduced in :cite:`Hoche2019ncc` aims for a consistent, combined prediction of a 4FS simulation and a merged 5FS prediction. The new simulation is hereby formed as sum of two components.
First, the direct component is given by the 4FS calculation. It is reformulated as if it would be part of a merged simulation and while doing so, the Sudakov rejection by the parton shower will actually remove events.
The fragmentation component then is given by the original, merged 5FS simulation. All configurations which are already described by the direct component are removed by an event filter.

Since the event filter operates on combined evolution histories, treating actual emissions from the parton shower and clustered emissions from matrix elements on the same footing, it must be ensured that the same cluster options are used everywhere.

The following parameters summarise all settings which where used for :cite:`Hoche2019ncc`.

The direct component
--------------------

.. literalinclude:: ./examples/V_plus_Bs/LHC_Fusing/Direct.yaml
      :language: yaml

Things to notice:

* The PDF and the strong coupling constant are evaluated in the 5FS.
  This is needed for consistency and to allow cluster configurations with initial-state bottom quarks.
  The resulting mismatch is corrected by counter-terms, implemented as user hook.
  If two strong couplings are involved at Born level, :option:`FUSING_DIRECT_FACTOR` must be set to one.
  If there are four such couplings, as in the case of ttbb, it must be two.

* From a physical point of view, no merging cut is needed for massive bottom quarks.
  For technical reasons, it nevertheless must be set to an unreachable value such as :option:`CKKW 100000.0`.


The fragmentation component
---------------------------

.. literalinclude:: ./examples/V_plus_Bs/LHC_Fusing/Fragmentation.yaml
      :language: yaml

Things to notice:

* The event filter can either be performed directly (default), or the veto information can be stored in HepMC.
  The latter requires :option:`HEPMC_USE_NAMED_WEIGHTS=1` and  :option:`FUSING_FRAGMENTATION_STORE_AS_WEIGHT=1`.
  The correspending weight is set to zero in case there is a fusing-induced event veto but remains unchanged if not.
