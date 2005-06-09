#ifndef Single_Vertex_h
#define Single_Vertex_h

#include "MyComplex.H"
#include "Flavour.H"
#include "Color_Function.H"
#include "Lorentz_Function.H"
#include <map>
//#include <vector>


namespace AMEGIC {
  class Single_Vertex {
  public:
    short int on;                 // flag if vertex is on
    std::string    Str;           // name of the verte
    int nleg, ncf, nlf, t;        // number of legs, color functions, lorentz
                                  //  functions and ?

    ATOOLS::Flavour   in[4];      // flavours at the legs 
    Complex   cpl[4];             // coupling constants
    Color_Function*   Color;      // color structure 
    Lorentz_Function* Lorentz;    // lorentz structure

    // constructors and destructor
    Single_Vertex();
    Single_Vertex(const Single_Vertex& v);
    ~Single_Vertex(); 
    
    // operators
    Single_Vertex&       operator=(const Single_Vertex&);
    int                  operator==(Single_Vertex&);
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
    \brief this file contains the class AMEGIC::Single_Vertex
  */
  /*! 
    \class Single_Vertex 
    \brief Contains all the information about a \f$ 1\to 2\f$ or \f$ 1\to 3\f$ 
    Feynman rule. Especially the three or four flavours, the associated color and 
    Lorentz structure, the left and right handed coupling constants and a
    representation by a charakter string.  
  */ 
}
#endif
