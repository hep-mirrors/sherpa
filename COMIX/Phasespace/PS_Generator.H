#ifndef COMIX_Phasespace_PS_Generator_H
#define COMIX_Phasespace_PS_Generator_H

#include "COMIX/Amplitude/Amplitude.H"

namespace PHASIC { class Cut_Data; }

namespace COMIX {

  class Process_Base;

  typedef std::vector<double> Double_Vector;

  class PS_Generator {
  public:

    typedef std::map<Current*,Current*>      CB_Map;
    typedef std::multimap<Current*,Current*> CB_MMap;
    typedef std::pair<Current*,Current*>     CB_Pair;

    typedef std::map<size_t,Current_Vector> TCC_Map;

  private:

    Process_Base *p_xs;

    Int_Matrix m_cl;

    Current_Matrix m_cur;
    Current_Vector m_ctt;

    size_t m_n, m_itmin, m_itmax;

    int m_zmode, m_ecmode, m_pmsinit, m_aniso;

    double m_thmass, m_chmass;

    Double_Vector m_smasses, m_swidths;

    CB_MMap m_cmap;
    CB_Map  m_cbmap;

    TCC_Map m_tccs;

    void CleanUp();
    void CalcJL();

    void AddSubChannel(NLO_subevtlist *const subs,const Vertex_Key &vkey);

    bool AddCurrent(Current *const ref,const ATOOLS::Flavour &fl,
		    const size_t &n,const int mode=0,
		    const double &m=0.0,const double &w=0.0,
		    Current *const scc=NULL);
    void AddExtraCurrent(Current *const cur,const size_t &n,
			 const double &m,const double &w,
			 Current *const scc=NULL);

    int DecayType(const Current *jc,
		  const Current *ja,const Current *jb) const;

    void AddSTCC();
    void AddSC();

  public:
    
    // constructor
    PS_Generator(Process_Base *const xs=NULL);

    // destructor
    ~PS_Generator();

    // member functions
    bool Construct(Amplitude *const ampl,
		   ATOOLS::NLO_subevtlist *const subs);
    void SetPrefMasses(PHASIC::Cut_Data *const cuts);

    bool Evaluate();

    void SetColors(const Int_Vector &rc,
		   const Int_Vector &ac);

    size_t NChannels() const;

    // inline functions
    inline const Current_Matrix &Graphs() const { return m_cur; }

    inline void SetZMode(const int &mode) { m_zmode=mode; }

    inline double ThresholdMass() const { return m_thmass; }

    inline const Double_Vector &ISRMasses() const { return m_smasses; }
    inline const Double_Vector &ISRWidths() const { return m_swidths; }

  };// end of class PS_Generator

}// end of namespace COMIX

#endif
