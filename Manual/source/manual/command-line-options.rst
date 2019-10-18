.. _Command line:

####################
Command Line Options
####################

The available command line options for Sherpa.

.. option:: --run-data, -f <file>

   Read settings from input file ``<file>``

.. option:: --path, -p <path>

   Read input file from path ``<path>``.

.. option:: --sherpa-lib-path, -L <path>

   Set Sherpa library path to ``<path>``, see
   :REF:`SHERPA_CPP_PATH`.

.. option:: --events, -e <N_events>

   Set number of events to generate ``<N_events>``, see
   :REF:`param_EVENTS`.

.. option:: --event-type, -t <event_type>

   Set the event type to ``<event_type>``, see :ref:`EVENT_TYPE`.

.. option:: --result-directory, -r <path>

   Set the result directory to ``<path>``, see
   :REF:`RESULT_DIRECTORY`.

.. option:: --random-seed, -R <seed>

   Set the seed of the random number generator to ``<seed>``, see
   :REF:`RANDOM_SEED`.

.. option:: --me-generators, -m <generators>

   Set the matrix element generator list to ``<generators>``, see
   :REF:`ME_GENERATORS`. If you specify more than on generator, use the
   YAML sequence syntax, e.g. :option:`-m '[Amegic, Comix]'`.

.. option:: --mi-handler, -M <handler>

   Set multiple interaction handler to ``<handler>``, see :ref:`MI_HANDLER`.

.. option:: --event-generation-mode, -w <mode>

   Set the event generation mode to ``<mode>``,
   see :REF:`EVENT_GENERATION_MODE`.

.. option:: --shower-generator, -s <generator>

   Set the parton shower generator to ``<generator>``, see
   :REF:`SHOWER_GENERATOR`.

.. option:: --fragmentation, -F <module>

   Set the fragmentation module to ``<module>``, see
   :ref:`Fragmentation`.

.. option:: --decay, -D <module>

   Set the hadron decay module to ``<module>``, see :ref:`Hadron
   decays`.

.. option:: --analysis, -a <analyses>

   Set the analysis handler list to ``<analyses>``, see
   :REF:`ANALYSIS`.  If you specify more than one analysis handler, use
   the YAML sequence syntax, e.g. :option:`-a '[Rivet, Internal]'`.

.. option:: --analysis-output, -A <path>

   Set the analysis output path to ``<path>``, see
   :REF:`ANALYSIS_OUTPUT`.

.. option:: --output, -O <level>

   Set general output level ``<level>``, see :ref:`OUTPUT`.

.. option:: --event-output, -o <level>

   Set output level for event generation ``<level>``, see
   :REF:`OUTPUT`.

.. option:: --log-file, -l <logfile>

   Set log file name ``<logfile>``, see :ref:`LOG_FILE`.

.. option:: --disable-result-directory-generation, -g

   Do not create result directory, see :REF:`RESULT_DIRECTORY`.


.. option:: --disable-batch-mode, -b

   Switch to non-batch mode, see :REF:`BATCH_MODE`.

.. option:: --print-version-info, -V

   Print extended version information at startup.

.. option:: --version, -v

   Print versioning information.

.. option:: --help, -h

   Print a help message.

.. option:: 'PARAMETER: Value'

   Set the value of a parameter, see :ref:`Parameters`.

.. option:: 'Tags: {TAG: Value}'

   Set the value of a tag, see :ref:`Tags`.
