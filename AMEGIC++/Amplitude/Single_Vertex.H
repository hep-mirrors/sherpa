#ifndef AMEGIC__Amplitude__Single_Vertex_H
#define AMEGIC__Amplitude__Single_Vertex_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Color_Function.H"
#include "AMEGIC++/Amplitude/Lorentz_Function.H"
#include <vector>

namespace AMEGIC {
  class Single_Vertex {
  public:
    int               nleg;        // number of legs
    ATOOLS::Flavour   in[4];       // flavours at the legs 
    short int         on;          // flag if vertex is on
    int               t,dec;
    std::vector<int>  order;
 
    std::vector<ATOOLS::Kabbala> cpl;      // coupling constants
    std::vector<MODEL::Color_Function>   Color;       // color structure 
    std::vector<MODEL::Lorentz_Function*> Lorentz;     // lorentz structure

    // constructors and destructor
    Single_Vertex();
    Single_Vertex(const Single_Vertex& v);
    ~Single_Vertex(); 
    
    Complex Coupling(size_t i) const;
    std::string PID() const;

    int CheckCoupling() const;
    int Compare(const Single_Vertex *v) const;
    // operators
    Single_Vertex&       operator= (const Single_Vertex&);
    bool                 operator==(const Single_Vertex&);
    friend std::ostream &operator<<(std::ostream&, const Single_Vertex&);
  };   

  /*! 
    \file
    \brief this file contains the class AMEGIC::Single_Vertex
  */
  /*! 
    \class Single_Vertex 
    \brief Contains all the information about a \f$ 1\to 2\f$ or \f$ 1\to 3\f$ 
    Feynman rule. Especially the three or four flavours, the associated color and 
    Lorentz structure, the left and right handed coupling constants and a
    representation by a charakter string.  
  */ 
} // end of namespace
#endif
