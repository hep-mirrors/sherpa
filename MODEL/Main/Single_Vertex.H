#ifndef MODEL__Main__Single_Vertex_H
#define MODEL__Main__Single_Vertex_H

#include "MODEL/Main/Color_Function.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Phys/Flavour.H"

namespace MODEL {

  class Single_Vertex {
  public:

    std::vector<ATOOLS::Flavour> in;
    std::vector<ATOOLS::Kabbala> cpl;
    std::vector<Color_Function> Color;
    std::vector<std::string> Lorentz;
    std::vector<size_t> id;
    std::vector<int> order;

    int dec;

    // constructors and destructor
    Single_Vertex();
    ~Single_Vertex(); 
    
    std::string PID() const;

    int Compare(const Single_Vertex *v) const;

    bool PureQCD() const;

    inline size_t NLegs() const { return in.size(); }
    inline void AddParticle(const ATOOLS::Flavour &fl,const int n=-1)
    { in.push_back(fl); id.push_back(n>=0?n:in.size()-1); }

    inline const Complex &Coupling(const size_t &i) const
    { return cpl[i].Value(); }

    bool operator==(const Single_Vertex&);

  };   

  std::ostream &operator<<(std::ostream &,const Single_Vertex &);

} // end of namespace

#endif
