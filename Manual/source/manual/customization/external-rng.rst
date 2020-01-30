.. _External RNG:

************
External RNG
************

To use an external Random Number Generator (RNG) in Sherpa, you need
to provide an interface to your RNG in an external dynamic
library. This library is then loaded at runtime and Sherpa replaces
the internal RNG with the one provided.

*In this case Sherpa will not attempt to set, save, read or restore the RNG*

The corresponding code for the RNG interface is

.. code-block:: c++

   #include "ATOOLS/Math/Random.H"

   using namespace ATOOLS;

   class Example_RNG: public External_RNG {
   public:
     double Get()
     {
       // your code goes here ...
     }
   };// end of class Example_RNG

   // this makes Example_RNG loadable in Sherpa
   DECLARE_GETTER(Example_RNG,"Example_RNG",External_RNG,RNG_Key);
   External_RNG *ATOOLS::Getter<External_RNG,RNG_Key,Example_RNG>::operator()(const RNG_Key &) const
   { return new Example_RNG(); }
   // this eventually prints a help message
   void ATOOLS::Getter<External_RNG,RNG_Key,Example_RNG>::PrintInfo(std::ostream &str,const size_t) const
   { str<<"example RNG interface"; }

If the code is compiled into a library called libExampleRNG.so, then
this library is loaded dynamically in Sherpa using the command
``SHERPA_LDADD: ExampleRNG`` either on the command line or in
:file:`Sherpa.yaml`. If the library is bound at compile time, like e.g.
in cmt, you may skip this step.

Finally Sherpa is instructed to retrieve the external RNG by
specifying :option:`EXTERNAL_RNG: Example_RNG` on the command line or
in :file:`Sherpa.yaml`.
