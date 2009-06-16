#ifndef MODEL_Interaction_Models_Single_Vertex_H
#define MODEL_Interaction_Models_Single_Vertex_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "MODEL/Interaction_Models/Lorentz_Function.H"
#include <vector>


namespace MODEL {
  class Single_Vertex {
  public:
    int               nleg;        // number of legs
    ATOOLS::Flavour   in[4];       // flavours at the legs 
    short int         on;          // flag if vertex is on
    std::string       Str;         // name of the vertex
    int               t;
 
    std::vector<ATOOLS::Kabbala> cpl;      // coupling constants
    std::vector<Color_Function>   Color;       // color structure 
    std::vector<Lorentz_Function*> Lorentz;     // lorentz structure

    // constructors and destructor
    Single_Vertex();
    Single_Vertex(const Single_Vertex& v);
    ~Single_Vertex(); 
    
    Complex Coupling(size_t i) const;
    // operators
    Single_Vertex&       operator= (const Single_Vertex&);
    const bool           operator==(const Single_Vertex&);
    friend std::ostream& operator<<(std::ostream&, const Single_Vertex&);
};  


 //Format of a Single_Vertex to transport via MPI:
  struct MPI_Single_Vertex {    
    int                   m_fl[4];
    /*
    double                m_cpl[8];        
    MPI_Lorentz_Function  m_lf;
    MPI_Color_Function    m_cf;
    */
  };
  
  //for MPI usage
  void Single_Vertex2MPI(const Single_Vertex * , MPI_Single_Vertex &);
  
  Single_Vertex * MPI2Single_Vertex(const MPI_Single_Vertex &);
  
  std::ostream & operator<<(std::ostream &, const MPI_Single_Vertex &);

  /*! 
    \file
    \brief this file contains the class MODEL::Single_Vertex
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
