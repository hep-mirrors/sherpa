%{
#include <fstream>
#include <stddef.h>
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;
 %}

namespace ATOOLS {

  class Random : public Terminator_Object {

  public:

    Random(long nid);  // initialization for Ran2()
    Random(unsigned int i1,unsigned int i2,unsigned int i3,
	   unsigned int i4,unsigned int i5,unsigned int i6);

    ~Random();
    
    void SetSeed(long nid);  // seed for Rnd2()
    void SetSeed(unsigned int i1,unsigned int i2,
		 unsigned int i3,unsigned int i4);
    inline long int GetSeed() { return m_id; }

    double Get();   

  };// end of class Random

}

