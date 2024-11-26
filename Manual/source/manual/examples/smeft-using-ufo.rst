.. _BSM/UFO_SMEFTsim:

SMEFT using UFO
===============

This example uses a UFO model called SMEFTsim, which extends the Standard Model Lagrangian by introducing additional dimension-six operators. The new operators are scaled with their respective Wilson coefficients. Detailed information about the model can be found here: https://arxiv.org/abs/2012.11343.

In the corresponding Example directory
``<prefix>/share/SHERPA-MC/Examples/BSM/UFO_SMEFT/``
you will find an example run card to use this model.
Before using it, you will have to download the model definition e.g.
from (`<https://github.com/SMEFTsim/SMEFTsim/releases/latest>`_)
and make it available to Sherpa as described in :ref:`UFO Model Interface`
by executing

.. code-block:: shell-session

   $ cd <prefix>/share/SHERPA-MC/Examples/BSM/UFO_SMEFTsim/
   $ wget https://github.com/SMEFTsim/SMEFTsim/archive/refs/tags/v3.0.2.tar.gz
   $ tar xzf v3.0.2.tar.gz
   $ <prefix>/bin/Sherpa-generate-model SMEFTsim-3.0.2/UFO_models/SMEFTsim_topU3l_MwScheme_UFO
   $    # if it asks you about converting to Python3, reply "y"


The events produced with the example run card below for electroweak
same-sign W pair production will include all EFT contributions up to the
quadratic order. Sherpa has the ability to produce the Standard Model, linear
and quadratic contributions separately by utilizing the :option:`Order` and
:option:`Max_Order` settings as indicated:

.. literalinclude:: /../../Examples/BSM/UFO_SMEFT/Sherpa.yaml
   :language: yaml
