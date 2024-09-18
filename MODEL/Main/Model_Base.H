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
    bool m_elementary;
    Model_Arguments(const bool &elementary):
      m_elementary(elementary)
    {}
  };// end of struct Model_Arguments

  struct EWParameters {
    double m_sw2_r, m_cw2_r, m_aew_r, m_mw_r, m_mz_r,
      m_mt_r, m_mh0_r, m_cvev_r;
  };

  struct ew_scheme {
    enum code {
      UserDefined = 0,
      alpha0      = 1,
      alphamZ     = 2,
      Gmu         = 3,
      alphamZsW   = 4,
      alphamWsW   = 5,
      GmumZsW     = 6,
      GmumWsW     = 7,
      FeynRules   = 10,
      Undefined   = 99
    };
  };// end of struct ew_scheme

  std::ostream &operator<<(std::ostream &str,const ew_scheme::code &c);
  std::istream &operator>>(std::istream &str,ew_scheme::code &c);

  class Model_Base {
  public:

    typedef ATOOLS::Getter_Function<Model_Base,Model_Arguments> 
    Model_Getter_Function;

  protected:

    std::string               m_name;
    bool                      m_elementary;
    bool                      m_hasnegativecouplingorders;

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

    const PDF::ISR_Handler_Map* p_isrhandlermap;

    virtual void InitVertices() = 0;

    void RotateVertices();

    void CustomContainerInit();

    virtual void ParticleInit() = 0;
    void ReadParticleData();
    // Needs to be called by subclasses after initialising particle data
    virtual void RegisterDefaults() const;
    void AddStandardContainers();
    void InitMEInfo();

    // Model initialization
    virtual bool ModelInit() = 0;
    void SetAlphaQED(const double& aqed_def);
    void SetAlphaQEDByScale(const double& scale2);
    void SetAlphaQEDByInput(const std::string& tag);
    void SetAlphaQCD(const PDF::ISR_Handler_Map& isr, const double& alphaS);
    void SetRunningFermionMasses();
    void SetRunningBosonMasses();

    void ReadExplicitCKM(ATOOLS::CMatrix& CKM);
    void OutputCKM();

  public:
    Model_Base(bool);
    virtual ~Model_Base();

    bool ModelInit(const PDF::ISR_Handler_Map& isr) {
      p_isrhandlermap = &isr;
      return ModelInit();
    }

    static void ShowSyntax(const size_t mode);
    void InitializeInteractionModel();
    void GetCouplings(Coupling_Map &cpls);

    std::string Name() const { return m_name; } 

    inline tscmap *GetCouplings() { return &m_cplmap; }

    inline size_t MaxLegs(const ATOOLS::Flavour &fl) const
    { std::map<ATOOLS::Flavour,size_t>::const_iterator
	it(m_maxlegs.find(fl)); return it==m_maxlegs.end()?0:it->second; }

    virtual int                     ScalarNumber(const std::string) const;
    virtual double                  ScalarConstant(const std::string) const;
    virtual Complex                 ComplexConstant(const std::string) const;
    virtual ATOOLS::Function_Base * GetScalarFunction(const std::string);
    virtual double                  ScalarFunction(const std::string,double);
    virtual std::string             MappedLorentzName(const std::string& label) const {return label;}

    int MaxNumber() const;

    const std::vector<Single_Vertex> &Vertices() const;
    const std::vector<Single_Vertex> &OriginalVertices() const;

    // Make this map public, such that OL interface can loop over it
    // in order to synchronize UFO model parameters
    const ScalarConstantsMap & ScalarConstants() const {return *p_constants; }

    inline  const Vertex_Table *VertexTable() const { return &m_vtable; }

    virtual bool CheckFlavours(int nin, int nout, ATOOLS::Flavour* flavs);

    inline ATOOLS::Flavour_Vector IncludedFlavours() const { return m_fls; }

    inline VMIterator_Pair GetVertex(const std::string &key) const
    { return m_vmap.equal_range(key); }

    /**
     * Map Order keys to the internal ordering of the couplings
     *
     * The default implementation returns 0 for "QCD" and 1 for "EW".
     * Overwrite this in subclasses if different/additional orders are used
     * therein.
     **/
    virtual size_t IndexOfOrderKey(const std::string&) const;

    virtual void ResetVerticesWithEWParameters(const EWParameters&) {};

    const PDF::ISR_Handler_Map& ISRHandlerMap() const { return *p_isrhandlermap; }

    void CheckForNegativeCouplingOrders();
    inline bool HasNegativeCouplingOrders() const 
    { return m_hasnegativecouplingorders; }
  };

  extern Model_Base *s_model;

}

#endif


