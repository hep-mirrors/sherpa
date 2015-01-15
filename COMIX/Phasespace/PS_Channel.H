#ifndef COMIX_Phasespace_PS_Channel_H
#define COMIX_Phasespace_PS_Channel_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/CXXFLAGS.H"

namespace PHASIC { class Color_Integrator; }

namespace COMIX {

  class Process_Base;
  class PS_Current;
  class PS_Vertex;

  typedef std::vector<PHASIC::Vegas*>          Vegas_Vector;
  typedef std::map<std::string,PHASIC::Vegas*> Vegas_Map;

  typedef std::map<const Current*,PHASIC::Vegas*> CVegas_Map;
  typedef std::map<size_t,PHASIC::Vegas*>         IVegas_Map;
  typedef std::map<size_t,CVegas_Map>             ICVegas_Map;

  typedef std::map<ATOOLS::NLO_subevt*,CVegas_Map>  SCVegas_Map;
  typedef std::map<ATOOLS::NLO_subevt*,ICVegas_Map> SICVegas_Map;

  typedef std::map<size_t,std::string> PSId_Map;

  typedef std::map<size_t,std::vector<int> > CId_Map;

  typedef std::map<size_t,const METOOLS::Current*> STCC_Map;

  class PS_Channel: public PHASIC::Single_Channel {
  protected:

    Process_Base   *p_xs;
    Current_Matrix *p_cur;

    Vertex_Vector  m_vtc;
    Vegas_Map      m_vmap;

    SCVegas_Map  m_pcmap;
    IVegas_Map   m_pimap;
    ICVegas_Map  m_sicmap;
    SICVegas_Map m_ticmap;

    PHASIC::Cut_Data *p_cuts;

    ATOOLS::Vec4D_Vector m_p;

    Double_Vector m_rns, m_wrns;
    Vegas_Vector  m_vgs, m_wvgs;

    size_t m_n, m_nr, m_num, m_lid, m_rid;
    double m_texp, m_stexp, m_sexp, m_thexp, m_mfac;
    double m_aexp, m_srbase;

    int m_bmode, m_omode, m_tmode, m_vmode, m_zmode, m_czmode;
    int m_nvints, m_vsopt, m_nopt;

    SP(PS_Generator) p_gen;

    PSId_Map *p_psid;
    CId_Map  *p_cid;
    
    const std::string      &GetPSId(const size_t &id);
    const std::vector<int> &GetCId(const size_t &id);

    inline size_t CIdCount(const size_t &id) { return GetCId(id).size(); }

    size_t SId(const size_t &id) const;
    double SCut(const size_t &id);

    void FillMoms(const size_t &aid,Int_Vector &cur,size_t n);

    double PropMomenta(const PS_Current *cur,const size_t &id,
		       const double &smin,const double &smax,const double *rn);
    double PropWeight(const PS_Current *cur,const size_t &id,
		      const double &smin,const double &smax,const double &s);

    void TChannelBounds(const size_t &aid,const size_t &lid,
			double &ctmin,double &ctmax,
			const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			const double &s1,const double &s2);
    void SingleTChannelBounds(const size_t &aid,const size_t &lid,
			      double &ctmin,double &ctmax,
			      const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			      const double &s1,const double &s2,const int mode);
    void TChannelMomenta(PS_Current *cur,ATOOLS::NLO_subevt *const dip,
			 const size_t &id,const size_t &aid,
			 const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			 ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2,
			 const double &s1,const double &s2,
			 const double *rns);
    double TChannelWeight(PS_Current *cur,ATOOLS::NLO_subevt *const dip,
			  const size_t &id,const size_t &aid,
			  const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			  ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2);

    void SChannelBounds(const size_t &id,const size_t &lid,
			double &ctmin,double &ctmax);
    void SChannelMomenta(PS_Current *cur,const int type,
			 const ATOOLS::Vec4D &pa,ATOOLS::Vec4D &p1,
			 ATOOLS::Vec4D &p2,const double &s1,
			 const double &s2,const double *rns);
    double SChannelWeight(PS_Current *cur,const int type,
			  ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2);

    bool GeneratePoint(PS_Current *const ja,PS_Current *const jb,
		       PS_Current *const jc,PS_Vertex *const v,size_t &nr);
    double GenerateWeight(PS_Current *const ja,PS_Current *const jb,
			  PS_Current *const jc,PS_Vertex *const v,size_t &nr);
    
    bool GenerateWeight(PS_Current *const cur);
    bool GenerateWeight();

    bool GeneratePoint(const size_t &id,size_t &nr,
		       Vertex_Vector &v);
    bool GeneratePoint(Vertex_Vector v);

    bool GenerateChannel(Current *const cur,
			 Vertex_Vector &v);
    bool GenerateChannel(Vertex_Vector &v);
    bool GenerateChannels();

    bool Zero(METOOLS::Vertex *const vtx) const;

    PHASIC::Vegas *GetVegas(const std::string &tag,int ni=-1);

    PHASIC::Vegas *GetPVegas(const PS_Current *cur,const size_t &id);
    PHASIC::Vegas *GetSVegas(const size_t &id,const PS_Current *cur);
    PHASIC::Vegas *GetTVegas(const size_t &id,const PS_Current *cur,
			     ATOOLS::NLO_subevt *const dip);

  public :

    // constructor
    PS_Channel(const size_t &nin,const size_t &nout,
	       ATOOLS::Flavour *fl,Process_Base *const ps); 

    // destructor
    ~PS_Channel();

    // member functions
    void GeneratePoint(ATOOLS::Vec4D *p,PHASIC::Cut_Data *cuts,double *rn);
    void GenerateWeight(ATOOLS::Vec4D *p,PHASIC::Cut_Data *cuts);

    void ISRInfo(int &type,double &m,double &w);
    void ISRInfo(std::vector<int> &ts,
		 std::vector<double> &ms,std::vector<double> &ws) const;

    int  ChNumber();
    void SetChNumber(int n);

    std::string ChID();

    size_t NChannels() const;

    void AddPoint(double value);

    void MPISync();
    void Optimize();
    void EndOptimize();
    bool OptimizationFinished();

    void WriteOut(std::string pid);
    void ReadIn(std::string pid);

  };// end of class PS_Channel

}// end of namespace COMIX

#endif
