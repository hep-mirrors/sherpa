#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"

#include "EXTAMP/Main/External_ME_Interface.H"
#include "EXTAMP/Main/BVI_Process.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Color_Correlated_ME2.H"
#include "PHASIC++/Process/KP_Terms.H"

#include <assert.h>

namespace EXTAMP {

  double BVI_Process::m_NF = 5.0;

  BVI_Process::BVI_Process(const PHASIC::Process_Info& pi) :
    Process(pi)
  {

    m_subtype = ATOOLS::ToType<int>(rpa->gen.Variable("NLO_SUBTRACTION_SCHEME"));
      
    /* Load loop ME and color-correlated ME */
    PHASIC::Process_Info loop_pi(m_pinfo);
    loop_pi.m_fi.m_nlotype=ATOOLS::nlo_type::loop;
    loop_pi.m_maxcpl[1] = m_pinfo.m_mincpl[1];
    loop_pi.m_maxcpl[0] = m_pinfo.m_mincpl[0] + 1;

    p_loop_me = PHASIC::Virtual_ME2_Base::GetME2(loop_pi);  
    p_corr_me = PHASIC::Color_Correlated_ME2::GetME2(pi);

    if (!p_loop_me)
      THROW(not_implemented, "Couldn't find virtual ME for this process.");
    if (!p_corr_me)
      THROW(not_implemented, "Couldn't find color-correlated ME for this process.");

    p_loop_me->SetCouplings(m_cpls);
    p_corr_me->SetCouplings(m_cpls);

    /* Calculate number of massless quark flavours */
    m_NF = ATOOLS::Flavour(kf_quark).Size()/2.0;

    /* Initialize KP terms */
    std::vector<size_t> pl;
    for (size_t i(m_nin);i<m_nin+m_nout;++i)
      if (m_flavs[i].Strong()) pl.push_back(i);
    p_kpterms = new PHASIC::KP_Terms(this,ATOOLS::sbt::qcd,pl);
    /* No sensible default is set for Alpha */
    p_kpterms->SetAlpha(1.0);
    p_kpterms->SetIType(ToType<cs_itype::type>(rpa->gen.Variable("NLO_IMODE")));

    /* Tell the mewgt about multiple components */
    m_mewgtinfo.m_type = mewgttype::B|mewgttype::VI|mewgttype::KP ;

    m_beta0 =  11./3.*m_CA - 4.0/3.0*m_TR*m_NF;
  }

  
  BVI_Process::~BVI_Process()
  {
    if (p_loop_me) delete p_loop_me;
    if (p_corr_me) delete p_corr_me;
    if (p_kpterms) delete p_kpterms;
  }
  

  double BVI_Process::Partonic(const ATOOLS::Vec4D_Vector &p,
				const int mode)
  {
    DEBUG_FUNC(this);
    
    double B(0.0),V(0.0),I(0.0),KP(0.0);
				   
    /* Maybe move to PHASIC::Single_Process */
    ScaleSetter()->CalculateScale(p);
    double mur = p_scale->Scale(stp::ren,1);

    /* Calculate full ME color-correlation matrix */
    p_corr_me->Calc(p);

    /* Get squared Born from correlator ME */
    B = p_corr_me->GetBorn2();

    /* Calculate V */
    p_loop_me->SetRenScale(mur);
    p_loop_me->Calc(p,B);

    /* mode 0: loop ME is missing Born factor 
       mode 1: loop ME includes Born factor   
       mode 2: loop ME includes Born factor and integrated subtraction terms */
    switch(p_loop_me->Mode())
      {
      case 0:
	V = p_loop_me->AlphaQCD()/(2.0*M_PI) *  p_loop_me->ME_Finite() * B ; break;
      case 1:
	V = p_loop_me->AlphaQCD()/(2.0*M_PI) *  p_loop_me->ME_Finite()     ; break;
      default:
	THROW(not_implemented, "Loop ME mode not implemented: "+ATOOLS::ToString(p_loop_me->Mode()));
      }

    /* Calculate integrated subtraction terms I */
    I = Calc_I(p, mur);

    /* Calculate KP terms */
    if(m_flavs[0].Strong() || m_flavs[1].Strong())
      KP = Calc_KP(p);

    std::pair<double,double> scaleterms = Calc_ScaleDependenceTerms(p, mur);

    /* Now divide all components by the symfac */
    B  /= NormFac(); V  /= NormFac(); I  /= NormFac(); KP /= NormFac();
    scaleterms.first /= NormFac(); scaleterms.second /= NormFac();
    
    /* Store all XS components in ME weight info */
    m_mewgtinfo.m_B       = B  ;
    m_mewgtinfo.m_VI      = V+I;
    m_mewgtinfo.m_KP      = KP ;
    m_mewgtinfo.m_wren[0] = scaleterms.first;
    m_mewgtinfo.m_wren[1] = scaleterms.second;

    /* Results to debugging output */
    msg_Debugging() << "Results of " << METHOD << "() {"
		    << "\n  B       = " << B
		    << "\n  V       = " << V
		    << "\n  I       = " << I
		    << "\n  KP      = " << KP
		    << "\n  wren[0] = " << scaleterms.first
		    << "\n  wren[1] = " << scaleterms.second
		    << "\n}" << std::endl;

    /* Store born in m_lastbxs (used in PHASIC::Single_Process) */
    m_lastbxs        = B;

    /* Store full XS in m_lastxs (used in PHASIC::Single_Process) */
    return m_lastxs = (B + V + I + KP);
  }


  double BVI_Process::Calc_I(const ATOOLS::Vec4D_Vector& p,
			     const double& mur) const
  {
    double I(0.);
    /* Loop over all pairs of partons i,j */
    for(std::vector<size_t>::const_iterator i=PartonIndices().begin(); i!=PartonIndices().end(); ++i)
      for(std::vector<size_t>::const_iterator j=i+1; j!=PartonIndices().end(); ++j)
	{
	  const ATOOLS::Flavour& fl_i(m_flavs[*i]);
	  const ATOOLS::Flavour& fl_j(m_flavs[*j]);
	  
	  /* This assumes that p_corr_me->Calc(p) has already been called! */
	  double M_ij = p_corr_me->GetValue(*i,*j)/Ti2(fl_i); // <m|TiTj/Ti^2|m>
	  double M_ji = p_corr_me->GetValue(*j,*i)/Ti2(fl_j); // <m|TjTi/Tj^2|m>
	  
	  /* Eps_Scheme_Factor usually defaults to 4pi */
	  double logf = log(4.0*M_PI*mur/(2.0*p[*i]*p[*j])/p_loop_me->Eps_Scheme_Factor(p));
	  
	  /* finite part of 
	     [ 4\pi\mu^2_r / (2*p_i*p_k*Eps_Scheme_Factor) ]^\epsilon * V_{i/j}(\epsilon), 
	     assume poles to be included in virtuals, rendering them finite */
	  double Vi_fin = (Vi_eps0(fl_i, m_subtype) + Vi_eps1(fl_i)*logf +  0.5*Vi_eps2(fl_i)*sqr(logf));
	  double Vj_fin = (Vi_eps0(fl_j, m_subtype) + Vi_eps1(fl_j)*logf +  0.5*Vi_eps2(fl_j)*sqr(logf));
	  I += Vi_fin*M_ij + Vj_fin*M_ji;
	  
	}
    
    /* Do not divide by symfac at this stage, this is done for all
       components simultaneously in Partonic */
    return -p_corr_me->AlphaQCD()/(2.0*M_PI) * I;
  }


  std::pair<double,double> BVI_Process::
  Calc_ScaleDependenceTerms(const ATOOLS::Vec4D_Vector& p,
			    const double& mur) const
  {
    /* First item: first derivative of all terms with respect to logf,
       second item: second derivative of all terms with respect to logf */
    std::pair<double,double> terms =  std::make_pair(0.0,0.0);

    /* Terms in virtual correction */
    double B = p_corr_me->GetBorn2();
    assert(MinOrder(0) == MaxOrder(0)-1);
    switch(p_loop_me->Mode())
      {
      case 0:
    	terms.first  += B*p_loop_me->ME_E1() + B*MinOrder(0)*m_beta0/2.0;
    	terms.second += B*p_loop_me->ME_E2();
    	break;
      case 1:
    	terms.first  +=   p_loop_me->ME_E1() + B*MinOrder(0)*m_beta0/2.0;
    	terms.second +=   p_loop_me->ME_E2();
    	break;
      default:
    	THROW(not_implemented, "Not implemented");
      }
   
    /* Terms in integrated dipoles: compare to Calc_I */
    for(std::vector<size_t>::const_iterator i=PartonIndices().begin(); i!=PartonIndices().end(); ++i)
      for(std::vector<size_t>::const_iterator j=i+1; j!=PartonIndices().end(); ++j)
	{
	  const ATOOLS::Flavour& fl_i(m_flavs[*i]);
	  const ATOOLS::Flavour& fl_j(m_flavs[*j]);
	  
	  /* This assumes that p_corr_me->Calc(p) has already been called! */
	  double M_ij = p_corr_me->GetValue(*i,*j)/Ti2(fl_i); // <m|TiTj/Ti^2|m>
	  double M_ji = p_corr_me->GetValue(*j,*i)/Ti2(fl_j); // <m|TjTi/Tj^2|m>
	  
	  double logf = log(4.0*M_PI*mur/(2.0*p[*i]*p[*j])/p_loop_me->Eps_Scheme_Factor(p));

	  terms.first  -= M_ij*(Vi_eps1(fl_i) +  Vi_eps2(fl_i)*logf);
	  terms.first  -= M_ji*(Vi_eps1(fl_j) +  Vi_eps2(fl_j)*logf);

	  terms.second -= M_ij*Vi_eps2(fl_i);
	  terms.second -= M_ji*Vi_eps2(fl_j);
	}
    
    /* Do not divide by symfac at this stage, this is done for all
       components simultaneously in Partonic */
    terms.first  *= p_corr_me->AlphaQCD()/(2.0*M_PI);
    terms.second *= p_corr_me->AlphaQCD()/(2.0*M_PI);

    return terms;
  }

  
  double BVI_Process::Vi_eps0(const ATOOLS::Flavour& flav, int subtype)
  {
    if(subtype==1)
      {
	if(flav.IsGluon())
	  return m_CA*(50./9.-M_PI*M_PI/2.0 + 1./36.) +  m_NF  * m_TR * (-16./9.-1./18.);
	if(flav.IsQuark())
	  return m_CF*(5.0-sqr(M_PI)/2.0 - 1./4.);
      }
    else if(subtype==0||subtype==2)
      {
	if(flav.IsGluon())
	  return m_CA*(50./9.-M_PI*M_PI/2.0) +  m_NF  * m_TR * (-16./9.);
	if(flav.IsQuark())
	  return m_CF*(5.0-sqr(M_PI)/2.0);
      }
    THROW(not_implemented, "Subtraction scheme not implemented");
  }

  
  double BVI_Process::Vi_eps1(const ATOOLS::Flavour& flav)
  {
    if(flav.IsGluon())
      return m_CA*11./6.     +  m_NF  * m_TR * (-2./3.);
    if(flav.IsQuark())
      return m_CF*3.0/2.0;
    THROW(fatal_error, "Internal error");
  }

  
  double BVI_Process::Vi_eps2(const ATOOLS::Flavour& flav)
  {
    if(flav.IsGluon())
      return m_CA;
    if(flav.IsQuark())
      return m_CF;
    THROW(fatal_error, "Internal error");
  }


  double BVI_Process::Ti2(const ATOOLS::Flavour& flav)
  {
    if(flav.IsGluon())
      return m_CA;
    if(flav.IsQuark())
      return m_CF;
    THROW(fatal_error, "Internal error");
  }

  
  double BVI_Process::Calc_KP(const ATOOLS::Vec4D_Vector& p)
  {
    /* Calculate partonic momentum fractions of incoming partons*/
    m_eta0 = (p[0][3]>0.0?p[0].PPlus()/rpa->gen.PBeam(0).PPlus():
	      p[0].PMinus()/rpa->gen.PBeam(1).PMinus());
    m_eta1 = (p[1][3]<0.0?p[1].PMinus()/rpa->gen.PBeam(1).PMinus():
	      p[1].PPlus()/rpa->gen.PBeam(0).PPlus());

    /* Randomly select x0 \in [eta0,1] 
                       x1 \in [eta1,1] 
       and calc weight corresponding to 
       the volume of those intervals */
    double w(1.0);
    if(m_flavs[0].Strong())
      {
	m_x0 = m_eta0+ATOOLS::ran->Get()*(1.-m_eta0);
	//m_x0 = m_eta0+0.2               *(1.-m_eta0); // set ran to 0.2 for unit tests
	w *= (1.-m_eta0);
      }
    if(m_flavs[1].Strong())
      {
	m_x1 = m_eta1+ATOOLS::ran->Get()*(1.-m_eta1);
	//m_x1 = m_eta1+               0.8*(1.-m_eta1); // set ran to 0.8 for unit tests
	w *= (1.-m_eta1);
      }

    /* Set all relevant members of ME_Weight_Info */
    bool swap = p_int->Momenta()[0][3]<p_int->Momenta()[1][3];
    m_mewgtinfo.m_swap = swap;
    m_mewgtinfo.m_y1   = swap?m_x1:m_x0;
    m_mewgtinfo.m_y2   = swap?m_x0:m_x1;
    p_kpterms->FillMEwgts(m_mewgtinfo);

    /* Populate a 2D array with color correlated MEs for KP_Terms.
       This assumes that p_corr_me->Calc(p) has already been called!
       Diagonal entries must remain zero exept for the [0][0] entry,
       which holds the squared born by convention. Also note 
       KP_Terms expect: dsij[i][j] = M(PartonIndices()[i], PartonIndices()[j])
       NOT:             dsij[i][j] = M(i,                  j                 ) */
    std::vector<std::vector<double> >
      dsij(PartonIndices().size(),std::vector<double>(PartonIndices().size(), 0.0));
    for(size_t i(0); i<PartonIndices().size(); i++)
      for(size_t j(0); j<PartonIndices().size(); j++)
	if(i!=j) dsij[i][j] = p_corr_me->GetValue(PartonIndices()[i],PartonIndices()[j]);
    dsij[0][0] = p_corr_me->GetBorn2();

    /* Let KP_Terms class calculate */
    double cpl = p_corr_me->AlphaQCD()/(2.0*M_PI);
    p_kpterms->Calculate(p,dsij,m_x0,m_x1,m_eta0,m_eta1,w*cpl);

    /* Do not divide by symfac at this stage, this is done for all
       components simultaneously in Partonic */
    double muf2(ScaleSetter()->Scale(stp::fac,1));
    return p_kpterms->Get(p_int->ISR()->PDF(0),p_int->ISR()->PDF(1),
			  m_x0,m_x1,m_eta0,m_eta1,
			  muf2,muf2,1.0,1.0,m_flavs[0],m_flavs[1]);
  }

  
  double BVI_Process::KPTerms(int mode, double scalefac2)
  {
    /* Used by PHASIC::Single_Process for reweighting, so have to
       include the normalization factor here */
    double muf2(ScaleSetter()->Scale(stp::fac,1));
    return p_kpterms->Get(p_int->ISR()->PDF(0),p_int->ISR()->PDF(1),
			  m_x0,m_x1,m_eta0,m_eta1,
			  muf2,muf2,scalefac2,scalefac2,m_flavs[0],m_flavs[1])/NormFac();
  }



  void BVI_Process::SetNLOMC(PDF::NLOMC_Base *const mc)
  {
    p_kpterms->SetNLOMC(mc);
    m_subtype = mc->SubtractionType();
  }

}
