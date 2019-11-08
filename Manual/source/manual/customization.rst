.. _Customization:

#############
Customization
#############

.. index:: SHERPA_LDADD

Customizing Sherpa according to your needs.

Sherpa can be easily extended with certain user defined tools.
To this extent, a corresponding C++ class must be written,
and compiled into an external library:

.. code-block:: shell-session

   $ g++ -shared \
         -I`$SHERPA_PREFIX/bin/Sherpa-config --incdir` \
         `$SHERPA_PREFIX/bin/Sherpa-config --ldflags` \
         -o libMyCustomClass.so My_Custom_Class.C

This library can then be loaded in Sherpa at runtime with the option
``SHERPA_LDADD``, e.g.:

.. code-block:: yaml

   SHERPA_LDADD:
   - MyCustomClass

Several specific examples of features which can be extended in this way are
listed in the following sections.

.. contents::
   :local:

.. toctree::

   customization/exotic-physics
   customization/custom-scale-setter
   customization/external-one-loop-me
   customization/external-rng
   customization/external-pdf
   customization/python-interface
