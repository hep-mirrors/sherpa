#ifndef AMISIC_Perturbative_MI_Process_Group_H
#define AMISIC_Perturbative_MI_Process_Group_H

#include "AMISIC++/Perturbative/MI_Process.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include <list>

namespace AMISIC {
  class MI_Process_Group {
  protected:
    std::string        m_name;
    scale_scheme::code m_muR_scheme;
    double             m_muR_fac, m_muF_fac;
    double             m_lastxs, m_pt02, m_scale, m_pref;

    MODEL::Running_AlphaS   * p_alphaS;
    MODEL::Running_AlphaQED * p_alpha;
    PDF::PDF_Base           * p_pdf[2];

    std::list<XS_Base * >    m_me2s;
    std::list<MI_Process * > m_processes;

    virtual double Scale(const double & pt2) const;
    virtual double Coupling(const double & scale) const=0;
    virtual double SoftCorrection(const double & pt2) const;
    virtual double PreCalculate(const double & shat,const double & that,
				const double & uhat);
  public:
    MI_Process_Group(const std::string & name);
    virtual ~MI_Process_Group();

    virtual double operator()(const double & shat,const double & that,
			      const double & uhat,const bool & out=false);
    virtual MI_Process * SelectProcess();
   
    inline const std::string & Name() const{ return m_name; }
    inline const double & LastXS()   const { return m_lastxs; }
    virtual void Output() const;
    
    inline void SetPDFs(PDF::PDF_Base * pdf0,PDF::PDF_Base * pdf1) {
      p_pdf[0] = pdf0; p_pdf[1] = pdf1;
    }
    inline void SetAlphaS(MODEL::Running_AlphaS * alpS)  { p_alphaS = alpS; }
    inline void SetAlpha(MODEL::Running_AlphaQED  * alp) { p_alpha  = alp; }
    inline void SetScale(const double & scale)           { m_scale = scale; }
  };

  class MI_GG_Processes: public MI_Process_Group {
  private:
    size_t m_Nqq;
    double Coupling(const double & scale) const;
  public:
    MI_GG_Processes();

    //double operator()(const double & shat,const double & that,
    //		      const double & uhat);
    MI_Process * SelectProcess();
  };

  class MI_QQB_Processes: public MI_Process_Group {
    double Coupling(const double & scale) const;
  public:
    MI_QQB_Processes();

    //double operator()(const double & shat,const double & that,
    //		      const double & uhat);
    // no SelectProcess(), fall back to generic solution:
    // smart specific solutions do not speed up the selection here
  };

  class MI_QQ_Processes: public MI_Process_Group {
    double Coupling(const double & scale) const;
  public:
    MI_QQ_Processes();

    //double operator()(const double & shat,const double & that,
    //		      const double & uhat);
    MI_Process * SelectProcess();
  };


  class MI_QG_Processes: public MI_Process_Group {
    double Coupling(const double & scale) const;
  public:
    MI_QG_Processes();

    //double operator()(const double & shat,const double & that,
    //		      const double & uhat);
    MI_Process * SelectProcess();
  };


  class MI_Q1Q2_Processes: public MI_Process_Group {
    double Coupling(const double & scale) const;
  public:
    MI_Q1Q2_Processes();

    //double operator()(const double & shat,const double & that,
    //		      const double & uhat);
    MI_Process * SelectProcess();
  };

  class MI_QG_QGamma_Processes: public MI_Process_Group {
    double Coupling(const double & scale) const;
  public:
    MI_QG_QGamma_Processes();

    //double operator()(const double & shat,const double & that,
    //		      const double & uhat);
    MI_Process * SelectProcess();
  };

  class MI_QQ_GGamma_Processes: public MI_Process_Group {
    double Coupling(const double & scale) const;
  public:
    MI_QQ_GGamma_Processes();

    //double operator()(const double & shat,const double & that,
    //		      const double & uhat);
    MI_Process * SelectProcess();
  };
}

#endif
