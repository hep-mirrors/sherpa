.. _Userhook

**********
User hooks
**********

The :option:`USERHOOKS` infrastructure in Sherpa can be used to implement a custom event afterburner. Written in C++ by the user, such a hook has access to the full Sherpa event structure and physics modules after the generation of the full hadron-level event. The event can be modified, rejected, or weights added before the event is passed to the analysis phases.

Please see :ref:`UserhookExamples` for an example how to use this.
