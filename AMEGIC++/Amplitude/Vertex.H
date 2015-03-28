#ifndef MODEL_Interaction_Models_Vertex_H
#define MODEL_Interaction_Models_Vertex_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Color_Function.H"
#include "MODEL/Main/Model_Base.H"
#include "AMEGIC++/Amplitude/Lorentz_Function.H"
#include "AMEGIC++/Amplitude/Single_Vertex.H"
#include <map>
#include <vector>


namespace AMEGIC {

  typedef  std::map<std::string,Complex> tscmap;

  typedef std::vector<Single_Vertex *> Vertex_List;
  typedef std::map<ATOOLS::Flavour, Vertex_List> Vertex_Table;

  class Vertex {
    int m_nvertex, m_n4vertex;       // number of single vertices included
    std::vector<Single_Vertex> m_v;  // list of 3-vertices 
    std::vector<Single_Vertex> m_v4; // list of 4-vertices 
    tscmap      m_cplmap;   
    
    void        CheckEqual(ATOOLS::Flavour**,short int&);
    void        GenerateVertex();
    int         CheckExistence(Single_Vertex&);  
    int         FermionRule(Single_Vertex& probe);
    int         SetVertex(Single_Vertex&,Single_Vertex&,int,int,int,int a=99,int mode=0);

    void        ColorExchange(MODEL::Color_Function*,int,int,int,int);
    void        LorentzExchange(MODEL::Lorentz_Function*,int,int,int,int);

  public:
    int on;
    Vertex(MODEL::Model_Base *);
    ~Vertex();

    void Print();
    

    inline Single_Vertex* operator[] (int n) {return &m_v[n];}  
    inline const Single_Vertex* operator[] (int n) const {return &m_v[n];}

    inline Single_Vertex* operator() (int n) {return &m_v4[n];}  
    inline const Single_Vertex* operator() (int n) const {return &m_v4[n];}

    inline int MaxNumber() const {return m_nvertex;}

    inline int MaxNumber4() const {return m_n4vertex;}
    
    inline tscmap* GetCouplings() { return &m_cplmap; }

    void TexOutput();

    int  FindVertex(Single_Vertex*);

    inline std::vector<Single_Vertex> &Vertices(int set)
    { return set?m_v4:m_v; }
  };
  
  std::ostream& operator<<(std::ostream&, const Vertex&);

  struct Amegic_Model {
  public:

    MODEL::Model_Base *p_model;
    Vertex *p_vertex;

    Amegic_Model(MODEL::Model_Base *model):
      p_model(model), p_vertex(new Vertex(model)) {}

    ~Amegic_Model() { delete p_vertex; }

  };


  /*! 
    \file
    \brief this file contains the class MODEL::Vertex   
  */

  /*! 
    \class Vertex 
    \brief This class generates all the interaction vertices present in a
    specific Model. 
  */
  
  /*!
    \fn Vertex::Vertex(Interaction_Model_Base*)
    
    \brief Builds the list of vertices and can produce a LaTeX output file.
    
    First of all , the lists m_v and m_v4 are initially filled via the Model routines of the model that's passed as the parameter: 
    
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
  
  \brief Generates the physically meaningfull permutations of the vertices created by the models.

  */
 /*!

  \fn int Vertex::SetVertex()
  
  
  */
 /*!

  \fn int Vertex::CheckExistence(Single_Vertex& probe)
  
  \brief Checks if the Single_Vertex probe already exists in the lists of
  vertices m_v or m_v4. 

  As a note here: In constrast to intuition this class returns true, if the vertex is not in the list.
  */
 
}
#endif
