#ifndef Vertex_h
#define Vertex_h

#include "MyComplex.H"
#include "Flavour.H"
#include "Color_Function.H"
#include "Lorentz_Function.H"


namespace AMEGIC {

  class Single_Vertex {
  public:
    /*!
      \brief Array of Flavour objects according to the particles 
      participating in the considered interaction
    */
    ATOOLS::Flavour   in[4];
    //!number of external legs the Single_Vertex possesses
    int       nleg;
    //!On/Off switch of a Single_Vertex 
    short int on;
    //!Array of complex couplings 
    Complex   cpl[4];
    /*!
      \brief 
      String expression associated to the coupling of the Single_Vertex object,
      here latex commands should be used
    */
    std::string    Str;
    //!the \f$SU(3)_C\f$ color structure of the interaction
    Color_Function*   Color;
    //!number of Color_Functions the Single_Vertex possesses
    int               ncf;
    //!the Lorentz structure of the interaction
    Lorentz_Function* Lorentz;
    //!number of Lorentz_Functions the SIngle_Vertex possesses
    int               nlf;
    
    /*!
      \brief The constructor initializes the number of external legs as 3 and 
      alocates a Color_Function and Lorentz_Function object
    */

    Single_Vertex() {
      nleg    = 3;
      Color   = 0;
      Lorentz = 0;
    }

    //! Copy constructor
    Single_Vertex(const Single_Vertex& v)
    {
      *this=v;
      
    }
    //! The destructor 
    
    ~Single_Vertex() {
      if (Color)   if (ncf==1) delete Color;
      else delete [] Color;
      if (Lorentz) if (nlf==1) delete Lorentz;
      else delete [] Lorentz;
    }
    
    

    Single_Vertex& operator=(const Single_Vertex& v) {
      if (Color)   if (ncf==1) delete Color;
      else delete [] Color;
      if (Lorentz) if (nlf==1) delete Lorentz;
      else delete [] Lorentz;
      if (this!=&v) {
	for (short int i=0;i<4;i++) in[i]  = v.in[i];
	for (short int j=0;j<4;j++) cpl[j] = v.cpl[j];
	
	nleg = v.nleg;
	Str  = v.Str;
	on   = v.on;
	ncf  = v.ncf;
	nlf  = v.nlf;
	if (ncf==1) Color   = new Color_Function(*v.Color);
	else {
	  Color = new Color_Function[ncf];
	  for (int i=0;i<ncf;i++)
	    Color[i] = v.Color[i];
	}
	if (nlf==1) Lorentz = new Lorentz_Function(*v.Lorentz);
	else {
	  Lorentz = new Lorentz_Function[nlf];
	  for (int i=0;i<nlf;i++)
	    Lorentz[i] = v.Lorentz[i];
	}
      }
      return *this;
    }
  };
  
  class Interaction_Model_Base;

  class Vertex {
    //! list of all vertices with three external legs
    Single_Vertex* m_v;
    //
    //!The total number of \f$ 1\to 2\f$ interactions present 
    int            m_nvertex;
    //! list of all vertices with four external legs
    Single_Vertex* m_v4;
    //!The total number of \f$ 1\to 3\f$ interactions present 
    int            m_n4vertex;
    
    void           CheckEqual(ATOOLS::Flavour**,short int&);
    void           GenerateVertex();
    int            CheckExistence(Single_Vertex&);
    int            FermionRule(Single_Vertex& probe);
    int            SetVertex(Single_Vertex&,Single_Vertex&,int,int,int,int a=99);

    void           ColorExchange(Color_Function*,int,int,int,int);
    void           LorentzExchange(Lorentz_Function*,int,int,int,int);
    //    void           Conjugate(Color_Function*);
  public:
    int on;
    Vertex(Interaction_Model_Base *);
    void Print();
    ~Vertex();
    
    //! The access operators for three leg vertices
    inline Single_Vertex* operator[] (int n) {return &m_v[n];}  
    inline const Single_Vertex* operator[] (int n) const {return &m_v[n];}
    //! The access operators for four leg vertices
    inline Single_Vertex* operator() (int n) {return &m_v4[n];}  
    inline const Single_Vertex* operator() (int n) const {return &m_v4[n];}
    /*! 
      \brief The total number of \f$ 1\to 2\f$ interactions present (nvertex) 
      is returned  
    */
    inline const int MaxNumber() const {return m_nvertex;}
    /*! 
      \brief The total number of \f$ 1\to 3\f$ interactions present (n4vertex) 
      is returned
    */
    inline const int MaxNumber4() const {return m_n4vertex;}
    /*!
      \brief LaTeX output is generated for all Feynman rules. The 
      string expression for each Single_Vertex is combine with a
      FeynMF picture. 
    */
    void TexOutput();
    //! Adds a Single_Vertex to the list of vertices
    void AddVertex(Single_Vertex*);
    int  FindVertex(Single_Vertex*);
  };
  
  std::ostream& operator<<(std::ostream&, const Single_Vertex&);
  std::ostream& operator<<(std::ostream&, const Vertex&);

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
