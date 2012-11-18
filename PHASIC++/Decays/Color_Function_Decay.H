#ifndef Color_Function_Decay_h
#define Color_Function_Decay_h

#include "ATOOLS/Math/MyComplex.H"
#include <string>
#include <vector>

namespace MODEL {
  class Color_Function;
}

namespace PHASIC {

  class Color_Function_Decay : 
    public std::vector<std::pair<std::string,std::vector<int> > >
  {
    std::vector<int> m_internal;
    int m_max;
    static size_t m_nid;
  public:
    Color_Function_Decay();
    
    Color_Function_Decay(const MODEL::Color_Function& c1);

    void Conjugate();
    std::string String() const;
    std::vector<int> Multiply(const Color_Function_Decay& c);
    void Contract(int c1, int c2);
    void ReplaceIndex(int c1, int c2);
    void BumpIndices(size_t i, int bump);
    int HasIndex(int i) const;
    Complex Contract(const Color_Function_Decay& c);
    inline const std::vector<int>& Internal() const { return m_internal; }

    friend std::ostream &operator<<(std::ostream &os,
                                    const Color_Function_Decay &);

  };  

}  
#endif
