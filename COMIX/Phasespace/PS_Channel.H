#ifndef COMIX_Phasespace_PS_Channel_H
#define COMIX_Phasespace_PS_Channel_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__Threading 
#include <pthread.h> 
#endif 

namespace PHASIC { class Color_Integrator; }

namespace COMIX {

  class Process_Base;
  class PS_Current;

  typedef std::vector<PHASIC::Vegas*>          Vegas_Vector;
  typedef std::map<std::string,PHASIC::Vegas*> Vegas_Map;

  std::string PSId(const size_t &id);

  size_t IdCount(const size_t &id);

#ifdef USING__Threading
  class PS_Channel;

  struct CDBG_PS_TID {
    pthread_t m_id;
    PS_Channel *p_psc;
    size_t m_s, m_n, m_b, m_e, m_i;
    pthread_mutex_t m_s_mtx, m_t_mtx;
    pthread_cond_t m_s_cnd, m_t_cnd;
    CDBG_PS_TID(PS_Channel *const psc): 
      p_psc(psc), m_s(2), m_b(0), m_e(0) {}
  };// end of struct CDBG_PS_TID

  typedef std::vector<CDBG_PS_TID*> CDBG_PS_TID_Vector; 
#endif

  class PS_Channel: public PHASIC::Single_Channel {
  protected:

    Process_Base      *p_xs;
    Current_Matrix *p_cur;

    Vertex_Vector  m_vtc;
    Vegas_Map      m_vmap;

    PHASIC::Cut_Data *p_cuts;

    Vec4D_Vector  m_p;

    Double_Vector m_rns, m_wrns;
    Vegas_Vector  m_vgs, m_wvgs;

    size_t m_n, m_nr, m_cpr, m_ccr, m_num, m_nopt, m_lid, m_rid;
    double m_texp, m_sexp, m_thexp, m_mfac;

    int    m_bmode, m_omode, m_vmode, m_zmode, m_czmode, m_ismode;
    size_t m_nvints, m_vsopt;
    
#ifdef USING__Threading 
    CDBG_PS_TID_Vector m_cts;
    pthread_mutex_t    m_vgs_mtx, m_wvgs_mtx;

    static void *TGenerateWeight(void *arg);
#endif 

    size_t SId(const size_t &id) const;
    double SCut(const size_t &id) const;

    void FillMoms(const size_t &aid,Int_Vector &cur,size_t n);

    double PropMomenta(const Current_Base *cur,const size_t &id,
		       const double &smin,const double &smax,
		       const double *rn);
    double PropWeight(const Current_Base *cur,const size_t &id,
		      const double &smin,const double &smax,const double &s);

    void TChannelBounds(const size_t &aid,const size_t &lid,
			double &ctmin,double &ctmax,
			const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			const double &s1,const double &s2);
    void TChannelMomenta(Current_Base *cur,const size_t &id,const size_t &aid,
			 const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			 ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2,
			 const double &s1,const double &s2,
			 const double *rns);
    double TChannelWeight(Current_Base *cur,const size_t &id,const size_t &aid,
			  const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			  ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2);

    void SChannelBounds(const size_t &id,const size_t &lid,
			double &ctmin,double &ctmax);
    void SChannelMomenta(Current_Base *cur,const size_t &id,
			 const ATOOLS::Vec4D &pa,ATOOLS::Vec4D &p1,
			 ATOOLS::Vec4D &p2,const double &s1,
			 const double &s2,const double *rns);
    double SChannelWeight(Current_Base *cur,const size_t &id,
			  ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2);

    bool GeneratePoint(Current_Base *const ja,Current_Base *const jb,
		       Current_Base *const jc,size_t &nr);
    double GenerateWeight(Current_Base *const ja,Current_Base *const jb,
			  Current_Base *const jc,size_t &nr);
    
    bool GenerateWeight(PS_Current *const cur);
    bool GenerateWeight();

    bool GeneratePoint(const size_t &id,size_t &nr,
		       Vertex_Vector &v);
    bool GeneratePoint(Vertex_Vector v);

    bool GenerateChannel(Current_Base *const cur,
			 Vertex_Vector &v);
    bool GenerateChannel(Vertex_Vector &v);
    bool GenerateChannels();

    bool Zero(Vertex_Base *const vtx) const;

    PHASIC::Vegas *GetVegas(const std::string &tag);

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

    void AddPoint(double value);

    void Optimize();
    void EndOptimize();

    void WriteOut(std::string pid);
    void ReadIn(std::string pid);

    // inline functions
    inline size_t ISMode() const { return m_ismode; }

  };// end of class PS_Channel

}// end of namespace COMIX

#endif
