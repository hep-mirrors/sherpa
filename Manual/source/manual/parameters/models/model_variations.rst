Some Models can implement the ability to vary their external parameters on-the-fly and calculate the effects for the event generation.
In particular, all UFO-Models generated with the previously described interface can use these variations.
Although implemented for usage with UFO-Models, this feature can technically be implemented by any Model.
The variations are realized using reweighting of the events by recalculating the partonic cross-section.
To use the on-the-fly variation of a model parameter, add this to your configuration file:

.. code-block:: yaml

   MODEL_VARIATIONS:
      param_name1: [<value1>, <value2>, ...]
      param_name2:
         From: <value>
         To: <value>
         Step: <value>

The values to be varied over for a parameter can be specified by a list or in a range statement as shown above.
Multiple parameters for variation can be added, which will be combined for variation.
If the model does not implement the variations, they are ignored TODO.
Note that the variations are intended to work with Comix TODO.
As an example, consider the UFO-Model examples and vary a parameter from the parameter card.

