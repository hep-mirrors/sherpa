.. _Models:

******
Models
******

.. index:: MODEL

The main switch ``MODEL`` sets the model that Sherpa uses throughout
the simulation run. The default is :option:`SM`, the built-in Standard
Model implementation of Sherpa. For BSM simulations, Sherpa offers an
option to use the Universal FeynRules Output Format (UFO)
:cite:`Degrande2011ua`.

Please note: AMEGIC can only be used for the built-in models (SM and
HEFT). For anything else, please use Comix. For more details on the
Sherpa capabilities to simulate BSM physics see :cite:`Hoeche2014kca`.

.. contents::
   :local:
   :depth: 1

Built-in Models
===============

.. contents::
   :local:
   :depth: 1

.. include:: ./models/standard-model.rst
.. include:: ./models/effective-higgs-coupling.rst

.. _UFO Model Interface:

UFO Model Interface
===================

.. index:: UFO_PARAM_CARD

.. include:: ./models/ufo-interface.rst
