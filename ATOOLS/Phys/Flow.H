#ifndef ATOOLS_Phys_Flow_H
#define ATOOLS_Phys_Flow_H

#include <iostream>
#include <map>
#include <set>

namespace ATOOLS {

  class Particle;

  class Flow {
  private:

    static unsigned int s_qcd_counter;

    unsigned int m_codes[2];

    friend std::ostream& operator<<(std::ostream &ostr,const Flow &flow);

  public:

    Flow();
    Flow(const Flow &flow);
    ~Flow();

    void SetCode(const Flow &flow);

    /// Sets the code for code_index, which is expected to be either 1 or 2.
    /// If the code is left at its default, -1, a previously unused value is
    /// used as a colour code.
    void SetCode(const unsigned int code_index=1,const int code=-1);

    /// Returns the code for index, which is expected to be either 1 or 2.
    unsigned int Code(const unsigned int index=1) const;

    /// Returns the index for code, if it exists, or -1 otherwise
    int          Index(const unsigned int code) const;

    void         SwapColourIndices();

    inline static void ResetCounter() { s_qcd_counter=600; }
    inline static unsigned int Counter() { return ++s_qcd_counter; }

  };// end of class Flow

}// end of namespace ATOOLS

#endif
