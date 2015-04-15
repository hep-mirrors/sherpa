#ifndef MODEL_Main_Model_Base_H
#define MODEL_Main_Model_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Matrix.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "MODEL/Main/Coupling_Data.H"
#include "PDF/Main/ISR_Handler.H"
#include <map>

namespace ATOOLS { 
  class Data_Reader;
} 

namespace MODEL {

  class Single_Vertex;

  typedef std::map<std::string,Complex> tscmap;

  typedef std::vector<Single_Vertex *> Vertex_List;
  typedef std::map<ATOOLS::Flavour, Vertex_List> Vertex_Table;
  typedef std::multimap<std::string,Single_Vertex*> Vertex_Map;
  typedef std::pair<std::string,Single_Vertex*> VMap_Key;
  typedef std::pair<Vertex_Map::const_iterator,
		    Vertex_Map::const_iterator> VMIterator_Pair;

  typedef std::map<std::string,int>                     ScalarNumbersMap; 
  typedef std::map<std::string,double>                  ScalarConstantsMap; 
  typedef std::map<std::string,Complex>                 ComplexConstantsMap; 
  typedef std::map<std::string,ATOOLS::Function_Base *> ScalarFunctionsMap; 
  typedef std::map<std::string,ATOOLS::CMatrix>         ComplexMatricesMap; 

  struct Model_Arguments {
    std::string m_path, m_file;
    bool m_elementary;
    Model_Arguments(const std::string &path,const std::string &file,const bool &elementary):
      m_path(path), m_file(file), m_elementary(elementary) {}
  };// end of struct Model_Arguments

  class Model_Base {
  public:

    typedef ATOOLS::Getter_Function<Model_Base,Model_Arguments> 
    Model_Getter_Function;

  protected:

    std::string               m_dir,m_file;
    std::string               m_name;
    bool                      m_elementary;

    ATOOLS::Data_Reader     * p_dataread;
    ScalarNumbersMap        * p_numbers;
    ScalarConstantsMap      * p_constants;
    ComplexConstantsMap     * p_complexconstants;
    ScalarFunctionsMap      * p_functions;

    Vertex_Table              m_vtable;
    Vertex_Map                m_vmap;
    ATOOLS::Flavour_Vector    m_fls;

    std::vector<Single_Vertex> m_v, m_ov;

    std::map<ATOOLS::Flavour,size_t> m_maxlegs;

    tscmap m_cplmap;   

    virtual void InitVertices() = 0;

    void RotateVertices();

    void CustomContainerInit();

    virtual void ParticleInit() = 0;
    void ReadParticleData();
    void AddStandardContainers();
    void InitMEInfo();

    // to be called in ModelInit 
    void SetAlphaQCD(const PDF::ISR_Handler_Map& isr);
    void SetRunningFermionMasses();
    
  public:
    Model_Base(std::string,std::string,bool);
    virtual ~Model_Base();

    virtual bool ModelInit(const PDF::ISR_Handler_Map& isr) = 0;

    static void ShowSyntax(const size_t mode);
    void InitializeInteractionModel();
    void GetCouplings(Coupling_Map &cpls);

    std::string Name() const { return m_name; } 

    inline tscmap *GetCouplings() { return &m_cplmap; }

    inline size_t MaxLegs(const ATOOLS::Flavour &fl) const
    { std::map<ATOOLS::Flavour,size_t>::const_iterator
	it(m_maxlegs.find(fl)); return it==m_maxlegs.end()?0:it->second; }

    virtual int                     ScalarNumber(const std::string);
    virtual double                  ScalarConstant(const std::string);
    virtual Complex                 ComplexConstant(const std::string);
    virtual ATOOLS::Function_Base * GetScalarFunction(const std::string);
    virtual double                  ScalarFunction(const std::string,double);

    int MaxNumber() const;

    const std::vector<Single_Vertex> &Vertices() const;
    const std::vector<Single_Vertex> &OriginalVertices() const;

    inline  const Vertex_Table *VertexTable() const { return &m_vtable; }

    virtual bool CheckFlavours(int nin, int nout, ATOOLS::Flavour* flavs);

    inline ATOOLS::Flavour_Vector IncludedFlavours() const { return m_fls; }

    inline VMIterator_Pair GetVertex(const std::string &key) const
    { return m_vmap.equal_range(key); }

  };

  extern Model_Base *s_model;

}

#endif


