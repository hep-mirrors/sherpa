#ifndef Vertex_h
#define Vertex_h

#include "MyComplex.H"
#include "Flavour.H"
#include "Color_Function.H"
#include "Lorentz_Function.H"
#include <map>
#include <vector>


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
    Single_Vertex& operator=(const Single_Vertex&);
    int            operator==(Single_Vertex&);
};
  
  class Interaction_Model_Base;

  class Vertex {
    int m_nvertex, m_n4vertex;       // number of single vertices included
    std::vector<Single_Vertex> m_v;  // list of 3-vertices 
    std::vector<Single_Vertex> m_v4; // list of 4-vertices 
   
    
    void        CheckEqual(ATOOLS::Flavour**,short int&);
    void        GenerateVertex();
    int         CheckExistence(Single_Vertex&);
    int         FermionRule(Single_Vertex& probe);
    int         SetVertex(Single_Vertex&,Single_Vertex&,int,int,int,int a=99);

    void        ColorExchange(Color_Function*,int,int,int,int);
    void        LorentzExchange(Lorentz_Function*,int,int,int,int);

  public:
    int on;
    Vertex(Interaction_Model_Base *);
    ~Vertex();

    void Print();
    

    inline Single_Vertex* operator[] (int n) {return &m_v[n];}  
    inline const Single_Vertex* operator[] (int n) const {return &m_v[n];}

    inline Single_Vertex* operator() (int n) {return &m_v4[n];}  
    inline const Single_Vertex* operator() (int n) const {return &m_v4[n];}

    inline const int MaxNumber() const {return m_nvertex;}

    inline const int MaxNumber4() const {return m_n4vertex;}

    void TexOutput();

    int  FindVertex(Single_Vertex*);
  };
  
  std::ostream& operator<<(std::ostream&, const Single_Vertex&);
  std::ostream& operator<<(std::ostream&, const Vertex&);

  typedef  std::vector<Single_Vertex *>  Vertex_List;
  typedef  std::map<ATOOLS::Flavour, Vertex_List> Vertex_Table;

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
    \brief this file contains the classes AMEGIC::Single_Vertex and
    AMEGIC::Vertex   
  */
  /*! 
    \class Single_Vertex 
    \brief Contains all the information about a \f$ 1\to 2\f$ or \f$ 1\to 3\f$ 
    Feynman rule. Especially the three or four flavours, the associated color and 
    Lorentz structure, the left and right handed coupling constants and a
    representation by a charakter string.  
  */

  /*! 
    \class Vertex 
    \brief This class generates all the interaction vertices present in a
    specific Model. 
  */
  
  /*!
    \fn Vertex::Vertex()
    
    \brief Builds the list of vertices and can produce a LaTeX output file.
    
    First of all , the lists m_v and m_v4 are initially filled via the Model routines: 
    
    - c_FFV, 
    - c_FFS, 
    - c_FFS, 
    - c_VVV, 
    - c_VVS, 
    - c_SSS, 
    - c_VVVV, 
    - c_SSVV, 
    - c_SSSS, 
    - c_FFT, 
    - c_VVT, 
    - c_SST, 
    - c_VVVT, 
    - c_FFVT, 
    - c_FFVT.
    

    Since during the matching of the vertices onto the raw topologies (Amplitude_Generator) 
    not only the standard form of a vertex is needed, but also all its permutations of external
    flavours, taking into account possible changes in the coupling constants. This task is 
    fulfilled by a call of the method GenerateVertex(). After all we are left with the complete
    lists of vertices contributing, m_v and m_v4.
    
  */
   
  /*!

  \fn void Vertex::GenerateVertex()
  
  hallo

  */
 /*!

  \fn int Vertex::SetVertex()
  
  
  */
 /*!

  \fn int Vertex::CheckExistence(Single_Vertex& probe)
  
  \brief Checks if the Single_Vertex probe already exists in the lists of
  vertices m_v or m_v4. 

  If the Single_Vertex probe possesses four legs its four flavours are compared to the four 
  flavours of all vertices in the list m_v4. If a complete agreement can 
  be found a zero is returned. 

  For a Single_Vertex with three external legs the three flavours of probe 
  are tested against all vertices of the list m_v. In addition the permutation 
  of the two outgoing flavours is tested. If an agreement is observed a zero
  is returned.

  If the vertex probe passes all tests a one is returned.

  */
  /*!
    \fn int Vertex::FermionRule(Single_Vertex& probe)

    \brief For three leg vertices involving fermions a specific
    ordering of the external flavous is used. 

    Vertices with three external legs involving fermions are chosen
    such that fermionic anti-particles have to occupy position in[1] and
    fermionic particles position in[2]. Note that Majorana fermions
    are explicitly left out by this rule.
    
    If one of the above conditions is not fulfilled the method returns a zero.
  */
 
}
#endif
