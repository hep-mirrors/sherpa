
This example implements a custom event afterburner using the :option:`USERHOOKS` infrastructure in Sherpa. Written in C++ by the user, such a hook has access to the full Sherpa event structure and physics modules after the generation of the full hadron-level event. The event can be modified, rejected, or weights added before the event is passed to the analysis phases.

Example code for a simple user hook just counting some event statistics:

.. literalinclude:: /../../Examples/Userhook/Userhook_Example.C
   :language: c++

The user hook code has to be compiled into a shared library ``libSherpaUserhookExample.so``, e.g. with the following lines assuming ``Sherpa-config`` is in your ``$PATH``:

.. literalinclude:: /../../Examples/Userhook/Makefile
   :language: make

This user hook can then be enabled in the run card as follows:

.. literalinclude:: /../../Examples/Userhook/Sherpa.yaml
   :language: yaml

Note that the shared library with a corresponding name will be loaded dynamically at run time from the current working directory or a path in your ``$LD_LIBRARY_PATH``.
