This is an example for event generation in the MSSM using Sherpa's UFO
support. In the corresponding Example directory
``<prefix>/share/SHERPA-MC/Examples/BSM/UFO_MSSM/``, you will find a
directory ``MSSM`` that contains the UFO output for the MSSM
(`<https://feynrules.irmp.ucl.ac.be/wiki/MSSM>`_). To run the example,
generate the model as described in :ref:`UFO Model Interface` by
executing

.. code-block:: shell-session

   $ cd <prefix>/share/SHERPA-MC/Examples/BSM/UFO_MSSM/
   $ <prefix>/bin/Sherpa-generate-model MSSM

An example run card will be written to the working directory. Use this
run card as a template to generate events.
