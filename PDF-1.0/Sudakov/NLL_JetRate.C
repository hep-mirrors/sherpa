#include "NLL_JetRate.H"
#include "Message.H"

using namespace SHERPA;
using namespace ATOOLS;

// static variable
Func_List NLL_JetRate::integrals;

NLL_JetRate::NLL_JetRate(int njet ,double q0, double Q, int mass_flag) :
  m_massive(mass_flag),m_njet(njet),m_qmin(q0),m_qmax(Q),m_qmin_b(q0*0.5),m_qmax_b(Q),m_gauss(this),
  p_deltaq(NULL), p_deltag(NULL), p_gammaq(NULL),  p_gammag(NULL),  p_gammaf(NULL), 
  p_deltaqm(NULL), p_deltagm(NULL),  p_gammaqm(NULL),  p_gammafm(NULL)
{
  m_mode = nll::Rate;
  m_table =0;
}

NLL_JetRate::~NLL_JetRate() 
{
  for (Func_List::iterator i=integrals.begin(); i!=integrals.end();++i) 
    delete i->second;
  integrals.clear();
}


void NLL_JetRate::Init() 
{
  msg_Debugging()<<" NLL_JetRate::Init "<<std::endl
		 <<" qmin ="<<m_qmin<<std::endl
		 <<" qmax ="<<m_qmax<<std::endl;
  m_qmin_b=m_qmin;
  
  if (m_massive==0 || m_massive==1 || m_massive==2 || m_massive==3) {
    InitMassless();
  }
  if (m_massive==1 || m_massive==2 || m_massive==3 || m_massive==4) {
    InitMassive();
  }
}

void NLL_JetRate::InitMassless()
{

  if (m_njet<3) return;
  if (m_njet<4) return;

  InitIntegral(nll::GgDg,nll::IntGgDg,"intGgDg");
  InitIntegral(nll::GfDf,nll::IntGfDf,"intGfDf");

  //  InitIntegral(nll::GgDG,nll::IntGgDG,"intGgDG");
  //  InitIntegral(nll::GFDF,nll::IntGFDF,"intGFDF");

  if (m_njet<5) return;

  m_table=1;
  InitIntegral(nll::GgDgIntGgDg,nll::IntGgDgIntGgDg,"intGgDgIntGgDg");
  InitIntegral(nll::GgDgIntGfDf,nll::IntGgDgIntGfDf,"intGgDgIntGfDf");
  //  InitIntegral(nll::GfDfIntGstarDg,nll::IntGfDfIntGstarDg,"intGfDfIntGstarDg");
  InitIntegral(nll::GqqgDg,nll::IntGqqgDg,"intGqqgDg");
  InitIntegral(nll::GfDfIntGqqgDg,nll::IntGfDfIntGqqgDg,"intGfDfIntGqqgDg");

}

void NLL_JetRate::InitMassive()
{
  m_qmin_b=m_qmin;

  if (m_njet<3) return;
  if (m_njet<4) return;

  InitIntegral(nll::GgDG,nll::IntGgDG,"intGgDG");
  InitIntegral(nll::GFDF,nll::IntGFDF,"intGFDF");

  if (m_njet<5) return;
  std::cout<<" Massive 5 Jet not implemented "<<std::endl;

}

void NLL_JetRate::InitIntegral(nll::code id, nll::code int_id, std::string name) 
{
  if (id+10 != int_id) {
    msg.Out()<<" WARNING: something fishy in NLL_JetRate::InitIntegral("<<id<<","<<int_id<<","<<name<<")"<<std::endl;
  } 

  // check if already exists!

  // ---
  Fast_Function * fun = new Fast_Function;
  m_mode=int_id;
  //  fun->Init(*this,m_qmin,m_qmax,640);
  fun->Init(*this,m_qmin,m_qmax,180);
  integrals[name]=fun;

  msg_Debugging()<<"NLL_JetRate::InitIntegral : "<<std::endl
		 <<Integrate(id,m_qmax)<<" vs. "<<(*integrals[name])(m_qmax)
		 <<" and "<<Integrate(id,sqrt(0.5)*m_qmax)<<" vs. "
		 <<(*integrals[name])(sqrt(0.5)*m_qmax)<<std::endl;

}

void NLL_JetRate::SetSudakovs(NLL_Sudakov_Base * deltaq, NLL_Sudakov_Base * deltag,
			      NLL_Sudakov_Base * deltaqm, NLL_Sudakov_Base * deltagm)
{
  p_deltaq  = deltaq;
  p_deltag  = deltag;
  p_deltaqm = deltaqm;
  p_deltagm = deltagm;
}


void NLL_JetRate::SetGammas(NLL_Branching_Probability_Base * gammaq ,NLL_Branching_Probability_Base * gammag,
	       NLL_Branching_Probability_Base * gammaf,
	       NLL_Branching_Probability_Base * gammaqm ,NLL_Branching_Probability_Base * gammafm) 
{
  p_gammaq  = gammaq;
  p_gammag  = gammag;
  p_gammaf  = gammaf;  

  p_gammaqm = gammaqm;
  p_gammafm = gammafm;  
}


void NLL_JetRate::Rates(double & r2, double & r3, double & r4, double & r5)
{
  if (m_qmin!=m_qmin_b) {
    msg.Error()<<"Warning in NLL_JetRate::Rates : "<<std::endl
	       <<"  q0 changed from "<<m_qmin<<" to "<<m_qmin_b<<std::endl
	       <<"  deleting all integration tables "<<std::endl;
    for (Func_List::iterator i=integrals.begin(); i!=integrals.end();++i) 
      delete i->second;
    integrals.clear();
    Init();
  }

  r2=r3=r4=r5=0.;

  if (m_njet<2) return;

  double rb =0.,sudqm=0.,intGQDG=0.,intGQDGintggDG_GFDF=0.;
  double intGQDGintggDG_gfdf_GFDF=0.;
  double sudq = (*p_deltaq)(m_qmax,m_qmin);

  // ---- out gammaqm ----
  if (1) {
    double q=m_qmin * sqrt(m_qmax/m_qmin);
    //    double gq = (*p_gammaq)(q,m_qmax);
    //    double gg = (*p_gammag)(q,m_qmax);
    //    double dq = (*p_deltaq)(q,m_qmin);
    //    double dg = (*p_deltag)(q,m_qmin);
    double gqm=0, dqm=0, dgm =0;
    if (p_gammaqm) gqm=(*p_gammaqm)(q,m_qmax);
    if (p_deltaqm) dqm=(*p_deltaqm)(q,m_qmin);
    if (p_deltagm) dgm=(*p_deltagm)(q,m_qmin);
  }

  // ----------------------------------------
  // *** R_2^0 == tR_2^0
  if (m_massive==0 || m_massive==2 || m_massive==3) {
    r2   = sqr(sudq);
  }

  
  if (m_massive==1 || m_massive==3 || m_massive==4 ) {
    sudqm = (*p_deltaqm)(m_qmax,m_qmin);
    // *** R_2^M == tR_2^M
    double r2M = sqr(sudqm);
    if (m_massive==1 || m_massive==4) r2=r2M;
    else {
      double sig1234 =  9418.8 + 7325.37+ 9418.8 + 7325.37;
      double sig5     = 9317.61;

      rb = sig5/(sig1234+sig5);
      // **** R_2 ****
      r2 = (1.-rb)*r2 + rb*r2M;
      std::cout<<"  r2   = "<<r2<<std::endl;
    }
  }

  if (m_njet==2) return ; 

  // ----------------------------------------
  double intgqdg = 0.;
  double intgqDG = 0.;
  if (m_massive==0) { 
    // *** tR_3^0
    intgqdg = Integrate(nll::GqDg,m_qmax);
    r3   = 2. * sqr(sudq) * intgqdg;
    std::cout<<" tr3^0 = "<<r3<<std::endl;
  }
  else 
  if (m_massive==2 || m_massive==3) {
    // *** R_3^0
    intgqDG = Integrate(nll::GqDG,m_qmax); 
    r3   = 2. * sqr(sudq) * intgqDG;
    std::cout<<"  r3^0 = "<<r3<<std::endl;
  }
  if (m_massive==1 || m_massive==3 || m_massive==4) {
    // *** R_3^m ~ tR_3^m
    intGQDG = Integrate(nll::GQDG,m_qmax);
    double r3M = 2. * sqr(sudqm) * intGQDG;
    std::cout<<"  r3^M = "<<r3M<<std::endl;
    if (m_massive==1 || m_massive==4) r3=r3M;
    else {
      // **** R_3 ****
      r3 = (1.-rb)*r3 + rb*r3M;
      std::cout<<"  r3   = "<<r3<<std::endl;
    }
  }
  if (m_njet==3) return ;

  // ----------------------------------------
  double  intgqdgintggdg_gfdf = 0.;
  double  intgqDGintggDG_gfdf_GFDF = 0.;
  m_table=1;
  if (m_massive==0) { 
    // *** tR_4^0
    intgqdgintggdg_gfdf = Integrate(nll::GqDgIntGgDg_GfDf,m_qmax);
    r4 = 2. * sqr(sudq) * ( sqr(intgqdg) + intgqdgintggdg_gfdf);
    std::cout<<" tr4^0 = "<<r4<<std::endl; 
  }
  else if (m_massive==2 || m_massive==3) {
    // *** R_4^0
    intgqDGintggDG_gfdf_GFDF = Integrate(nll::GqDGIntGgDG_GfDf_GFDF,m_qmax);
    r4 = 2. * sqr(sudq) * ( sqr(intgqDG) + intgqDGintggDG_gfdf_GFDF);
    std::cout<<"  r4^0 = "<<r4<<std::endl; 
  }

  if (m_massive==1 || m_massive==3 || m_massive==4) {
    // *** tR_4^m
    double r4M=0.;
    if (m_massive==4) {
      intGQDGintggDG_GFDF = Integrate(nll::GQDGIntGgDG_GFDF,m_qmax);
      r4M = 2. * sqr(sudqm) * ( sqr(intGQDG) + intGQDGintggDG_GFDF);
      r4=r4M;
      std::cout<<" tr4^M = "<<r4<<std::endl; 
    }
    else {
      // *** R_4^m
      intGQDGintggDG_gfdf_GFDF = Integrate(nll::GQDGIntGgDG_GfDf_GFDF,m_qmax);
      r4M = 2. * sqr(sudqm) * ( sqr(intGQDG) + intGQDGintggDG_gfdf_GFDF);
      std::cout<<"  r4^M = "<<r4M<<std::endl;
      if (m_massive==1) r4=r4M;
      else {
	// **** R_4 ****
	r4 = (1.-rb)*r4 + rb*r4M;
	std::cout<<"  r4   = "<<r4<<std::endl;
      }
    }
  }

  if (m_njet==4) return ; 


  // ----------------------------------------
  // *** tR_5^0
  /*
  m_table=0;
  double r5b = Integrate(nll::GqDgIntgg2_gggg_ggff_ffsg,m_qmax);
  */
  m_table=1;
  /*
  double r5b = Integrate(nll::GqDgIntgg2_gggg_ggff_ffsg,m_qmax);
  r5 = sqr(sudq) * ( 4./3. * intgqdg * ( sqr(intgqdg) + 3. * intgqdgintggdg_gfdf) 
		     + r5b);
  std::cout<<" r5 = "<<r5<<std::endl;
  */
  double r5c = Integrate(nll::GqDgIntffgg_gg2_ggff_gggg_ggtg,m_qmax);
  r5 = sqr(sudq) * ( 4./3. * intgqdg * ( sqr(intgqdg) + 3. * intgqdgintggdg_gfdf) 
		     + r5c);
  std::cout<<" r5 = "<<r5<<std::endl;
  //  std::cout<<" r5b ="<<r5b<<std::endl;
  std::cout<<" r5c ="<<r5c<<std::endl;
  
}

double NLL_JetRate::Integrate(nll::code mode,double q) 
{
  nll::code smode=m_mode;
  double sqmax  = m_qmax;

  m_mode = mode;
  m_qmax = q;
  double sum=m_gauss.Integrate(m_qmin,m_qmax,1.e-4,1);
  m_qmax = sqmax;
  m_mode = smode;

  return sum;
}

double NLL_JetRate::Integrate(nll::code mode,double q,double q_b) 
{
  nll::code smode=m_mode;
  double sqmax   = m_qmax;
  double sqmax_b = m_qmax_b;

  m_mode   = mode;
  m_qmax   = q;
  m_qmax_b = q_b;
  double sum=m_gauss.Integrate(m_qmin,m_qmax,1.e-4,1);
  m_qmax   = sqmax;
  m_qmax_b = sqmax_b;
  m_mode = smode;

  return sum;
}

double NLL_JetRate::operator()(double q)
{
  switch (m_mode) {
  case nll::Rate:
    m_qmax=q;
    return 0.;
  case nll::GqDg:
    return (*p_gammaq)(q,m_qmax) * (*p_deltag)(q,m_qmin);
  case nll::GgDg:
    return (*p_gammag)(q,m_qmax) * (*p_deltag)(q,m_qmin);
  case nll::GqDG:
    return (*p_gammaq)(q,m_qmax) * (*p_deltagm)(q,m_qmin);
  case nll::GfDf:
    {
      double dq = (*p_deltaq)(q,m_qmin);
      double dg = 0;
      if (m_massive==0) dg= (*p_deltag)(q,m_qmin);
      else dg= (*p_deltagm)(q,m_qmin);
      return (*p_gammaf)(q,m_qmax) * dq*dq/dg;
    }
  case nll::GstarDg: // Note: this part depents on 2 variables! The integral cannot be put in an 1D-fast-table!
    return ( 2.*(*p_gammaq)(q,m_qmax) - (*p_gammag)(q,m_qmax) + (*p_gammag)(q,m_qmax_b)) * (*p_deltag)(q,m_qmin);
  case nll::GqqgDg:
    return ( 2.*(*p_gammaq)(q,m_qmax) - (*p_gammag)(q,m_qmax)) * (*p_deltag)(q,m_qmin);
  case nll::IntGqDg:
    return Integrate(nll::GqDg,q);
  case nll::IntGgDg:
    return Integrate(nll::GgDg,q);
  case nll::IntGfDf:
    return Integrate(nll::GfDf,q);
  case nll::IntGqqgDg:
    return Integrate(nll::GqqgDg,q);
  case nll::GqDgIntGgDg_GfDf:
    if (m_table)
      return (*p_gammaq)(q,m_qmax) * (*p_deltag)(q,m_qmin) * 
	( (*integrals["intGgDg"])(q) + (*integrals["intGfDf"])(q) );
    else
      return (*p_gammaq)(q,m_qmax) * (*p_deltag)(q,m_qmin) * 
	( Integrate(nll::GgDg,q) + Integrate(nll::GfDf,q) );
  case nll::GgDgIntGgDg:
    if (m_table)
      return (*p_gammag)(q,m_qmax) * (*p_deltag)(q,m_qmin) * (*integrals["intGgDg"])(q);
    else
      return (*p_gammag)(q,m_qmax) * (*p_deltag)(q,m_qmin) * Integrate(nll::GgDg,q) ;
  case nll::GgDgIntGfDf:
    if (m_table)
      return (*p_gammag)(q,m_qmax) * (*p_deltag)(q,m_qmin) * (*integrals["intGfDf"])(q);
    else
      return (*p_gammag)(q,m_qmax) * (*p_deltag)(q,m_qmin) * Integrate(nll::GfDf,q) ;
  case nll::GfDfIntGstarDg:
    {
      double dq = (*p_deltaq)(q,m_qmin);
      double dg = (*p_deltag)(q,m_qmin);
      return (*p_gammaf)(q,m_qmax) * dq*dq/dg * Integrate(nll::GstarDg,q);
    }
  case nll::GfDfIntGqqgDg:
    {
      double dq = (*p_deltaq)(q,m_qmin);
      double dg = (*p_deltag)(q,m_qmin);
      if (m_table) 
        return (*p_gammaf)(q,m_qmax) * dq*dq/dg * (*integrals["intGqqgDg"])(q);
      else
        return (*p_gammaf)(q,m_qmax) * dq*dq/dg * Integrate(nll::GqqgDg,q);
    }
  case nll::IntGgDgIntGgDg:
    return Integrate(nll::GgDgIntGgDg,q);
  case nll::IntGgDgIntGfDf:
    return Integrate(nll::GgDgIntGfDf,q);
  case nll::IntGfDfIntGstarDg:
    return Integrate(nll::GfDfIntGstarDg,q,q);
  case nll::IntGfDfIntGqqgDg:
    return Integrate(nll::GfDfIntGqqgDg,q);
  case nll::GqDgIntgg2_gggg_ggff_ffsg:
    if (m_table)
      return (*p_gammaq)(q,m_qmax) * (*p_deltag)(q,m_qmin) * 
	( sqr( (*integrals["intGgDg"])(q)) 
	  + 2. * (*integrals["intGgDgIntGgDg"])(q) 
	  + 4. * (*integrals["intGgDgIntGfDf"])(q)
	  + 2. * (*integrals["intGfDfIntGstarDg"])(q));
    else
      return (*p_gammaq)(q,m_qmax) * (*p_deltag)(q,m_qmin) * 
	( sqr( Integrate(nll::GgDg,q)) 
	  + 2. * Integrate(nll::GgDgIntGgDg,q) 
	  + 4. * Integrate(nll::GgDgIntGfDf,q)
	  + 2. * Integrate(nll::GfDfIntGstarDg,q,q));
  case nll::GqDgIntffgg_gg2_ggff_gggg_ggtg:
    if (m_table)
      return (*p_gammaq)(q,m_qmax) * (*p_deltag)(q,m_qmin) * 
	( sqr( (*integrals["intGgDg"])(q)) 
	  + 2. * (*integrals["intGgDgIntGgDg"])(q) 
	  + 2. * (*integrals["intGgDgIntGfDf"])(q)
	  + 2. * (*integrals["intGfDf"])(q) * (*integrals["intGgDg"])(q)
	  + 2. * (*integrals["intGfDfIntGqqgDg"])(q));
    else
      return (*p_gammaq)(q,m_qmax) * (*p_deltag)(q,m_qmin) * 
	( sqr( Integrate(nll::GgDg,q)) 
	  + 2. * Integrate(nll::GgDgIntGgDg,q) 
	  + 2. * Integrate(nll::GgDgIntGfDf,q)
	  + 2. * Integrate(nll::GfDf,q) * Integrate(nll::GgDg,q)
	  + 2. * Integrate(nll::GfDfIntGqqgDg,q));
  case nll::GQDG:
    return (*p_gammaqm)(q,m_qmax) * (*p_deltagm)(q,m_qmin);
  case nll::GgDG:
    return (*p_gammag)(q,m_qmax) * (*p_deltagm)(q,m_qmin);
  case nll::GFDF:
    {
      double dq = (*p_deltaqm)(q,m_qmin);
      double dg = 0;
      if (m_massive==0) dg= (*p_deltag)(q,m_qmin);
      else dg= (*p_deltagm)(q,m_qmin);
      return (*p_gammafm)(q,m_qmax) * dq*dq/dg;
    }
  case nll::IntGQDG:
    return Integrate(nll::GQDG,q);
  case nll::IntGgDG:
    return Integrate(nll::GgDG,q);
  case nll::IntGFDF:
    return Integrate(nll::GFDF,q);
  case nll::GQDGIntGgDG_GFDF:
    if (m_table)
      return (*p_gammaqm)(q,m_qmax) * (*p_deltagm)(q,m_qmin) * 
	( (*integrals["intGgDG"])(q) + 5.* (*integrals["intGFDF"])(q) );
    else
      return (*p_gammaqm)(q,m_qmax) * (*p_deltagm)(q,m_qmin) * 
	( Integrate(nll::GgDG,q) + 5.* Integrate(nll::GFDF,q) );
  case nll::GqDGIntGgDG_GfDf_GFDF:
    if (m_table)
      return (*p_gammaq)(q,m_qmax) * (*p_deltagm)(q,m_qmin) * 
	( (*integrals["intGgDG"])(q) + (*integrals["intGfDf"])(q) + (*integrals["intGFDF"])(q) );
    else
      return (*p_gammaq)(q,m_qmax) * (*p_deltagm)(q,m_qmin) * 
	( Integrate(nll::GgDG,q) + Integrate(nll::GfDf,q) + Integrate(nll::GFDF,q) );
  case nll::GQDGIntGgDG_GfDf_GFDF:
    if (m_table)
      return (*p_gammaqm)(q,m_qmax) * (*p_deltagm)(q,m_qmin) * 
	( (*integrals["intGgDG"])(q) + (*integrals["intGfDf"])(q) + (*integrals["intGFDF"])(q) );
    else
      return (*p_gammaqm)(q,m_qmax) * (*p_deltagm)(q,m_qmin) * 
	( Integrate(nll::GgDG,q) + Integrate(nll::GfDf,q) + Integrate(nll::GFDF,q) );
  default:
    msg.Out()<<"ERROR in NLL_JetRate::operator()("<<m_mode<<std::endl ;
    return 0.;
  }
  
}

