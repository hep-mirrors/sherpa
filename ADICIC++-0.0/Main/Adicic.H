//bof
//Version: 4 ADICIC++-0.0/2006/05/24

//The main Adicic object.



#ifndef _Adicic_H_
#define _Adicic_H_ _Adicic_H_


#include <cassert>
#include <string>
#include "Run_Parameter.H"
#include "Blob_List.H"
#include "Particle_List.H"
#include "Model_Base.H"    //Will give running alpha_S etc.
#include "Cascade_Handler.H"


#define ADICIC_OUTPUT ADICIC_OUTPUT
#undef  ADICIC_OUTPUT





namespace ADICIC {



  class Adicic {

  private:

    int                 m_parcountcorr;    //Particle counter correction.
    int                 m_startnum;
    std::size_t         m_total;
    std::size_t         m_fail;
    std::size_t         m_noem;
    Cascade_Handler*    p_handler;
    Cascade             m_cascade;

    ATOOLS::Blob_List*               p_blist;
    std::list<ATOOLS::Particle_List> l_paic;
    std::vector<ATOOLS::Particle*>   v_conn;
    std::vector<ATOOLS::Particle*>   v_resi;

    bool FExtract(ATOOLS::Blob* pfs);
    bool IExtract(ATOOLS::Blob* pis, ATOOLS::Blob* pisme, ATOOLS::Blob* pfs);

  public:

    Adicic();
    Adicic(const std::string&, MODEL::Model_Base*);
    ~Adicic();

    inline void     PrepareCascade();
    inline Cascade& GetCascade();

    int  PerformShowers();
    bool ExtractPartons(ATOOLS::Blob_List*);

  };



  //===========================================================================



  inline void     Adicic::PrepareCascade() {
    assert(m_cascade.Clear());
    p_blist=NULL;
    l_paic.clear();
  }
  inline Cascade& Adicic::GetCascade()     {
    return m_cascade;
  }



}    //eo namespace ADICIC





#endif    //eo _Adicic_H_



//eof
