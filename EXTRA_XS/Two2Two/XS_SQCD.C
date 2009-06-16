#include "EXTRA_XS/Two2Two/XS_SQCD.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Phys/Flow.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1q2_sQ1sQ2>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsQuark() && flavours[1].IsQuark() &&
	flavours[2].IsSquark() && flavours[3].IsSquark() &&
	flavours[0]!=flavours[1] && 
	((!flavours[0].IsAnti() && !flavours[1].IsAnti()) || 
	 (flavours[0].IsAnti() && flavours[1].IsAnti())) &&
	((SuperPartner(flavours[2],flavours[0])  && 
	  SuperPartner(flavours[3],flavours[1])) ||
	 (SuperPartner(flavours[3],flavours[0])  && 
	  SuperPartner(flavours[2],flavours[1]))) &&
	((abs(flavours[2].Kfcode())>2000000 && abs(flavours[3].Kfcode())>2000000) ||
	 (abs(flavours[2].Kfcode())<2000000 && abs(flavours[3].Kfcode())<2000000))) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1q2_sQ1sQ2(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1q2_sQ1sQ2::XS_q1q2_sQ1sQ2(const size_t nin,const size_t nout, 
			       const ATOOLS::Flavour *fl,
			       XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=fl[1].IsAnti();
  m_g3=std::abs(p_model->Constant("g_3"));
  m_mgluino2=sqr(Flavour(kf_Gluino).Mass());
  m_nstrong=4;
}

double XS_q1q2_sQ1sQ2::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  return sqr(m_g3*m_g3)*2.0/9.0*m_mgluino2*s/sqr(t-m_mgluino2);
}

bool XS_q1q2_sQ1sQ2::SetColours(double s,double t,double u) 
{ 
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  int r = !SuperPartner(p_flavours[2],p_flavours[0]);
  if (m_a==m_p) {
    /*
      0-----\   /-----2, if fl[0]==~fl[2]
             \ /
	      X  t
             / \
      1-----/   \-----3, if fl[1]==~fl[3]

      shower scale is u
    */
    p_colours[0][m_a] = p_colours[3-r][m_a] = Flow::Counter();
    p_colours[1][m_a] = p_colours[2+r][m_a] = Flow::Counter();
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
      pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq'->sqsq', set scale u "<<u<<"\n";
  }
  else {
    /*
      0-----+ +-----2
            | |
	    | |  t
            | |
      1-----+ +-----3

      shower scale is s
    */
    p_colours[0][m_a]   = p_colours[1][m_p]   = Flow::Counter();
    p_colours[2+r][m_a] = p_colours[3-r][m_p] = Flow::Counter();
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
      pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb'->sqsqb', set scale s "<<s<<"\n";
  }
  if (swap) SwapInOrder();
  return 1; 
}

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1q2_sQ1LsQ2R>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsQuark() && flavours[1].IsQuark() &&
	flavours[2].IsSquark() && flavours[3].IsSquark() &&
	flavours[0]!=flavours[1] &&
	((SuperPartner(flavours[2],flavours[0])  && 
	  SuperPartner(flavours[3],flavours[1])) ||
	 (SuperPartner(flavours[3],flavours[0])  && 
	  SuperPartner(flavours[2],flavours[1]))) &&
	((abs(flavours[2].Kfcode())>2000000 && abs(flavours[3].Kfcode())<2000000) ||
	 (abs(flavours[2].Kfcode())<2000000 && abs(flavours[3].Kfcode())>2000000))) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1q2_sQ1LsQ2R(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1q2_sQ1LsQ2R::XS_q1q2_sQ1LsQ2R(const size_t nin,const size_t nout, 
			       const ATOOLS::Flavour *fl,
			       XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=fl[1].IsAnti();
  m_g3=std::abs(p_model->Constant("g_3"));
  m_mgluino2=sqr(Flavour(kf_Gluino).Mass());
  m_msq32=sqr(fl[2].Mass());
  m_msq42=sqr(fl[3].Mass());
  m_nstrong=4;
}

double XS_q1q2_sQ1LsQ2R::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double pT2=(u*t-m_msq32*m_msq42)/s;
  return sqr(m_g3*m_g3)*2.0/9.0*s*pT2/sqr(t-m_mgluino2);
}

bool XS_q1q2_sQ1LsQ2R::SetColours(double s,double t,double u) 
{ 
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  int r = !SuperPartner(p_flavours[2],p_flavours[0]);
  if (m_a==m_p) {
    /*
      0-----\   /-----2, if fl[0]==~fl[2]
             \ /
	      X  t
             / \
      1-----/   \-----3, if fl[1]==~fl[3]

      shower scale is u
    */
    p_colours[0][m_a] = p_colours[3-r][m_a] = Flow::Counter();
    p_colours[1][m_a] = p_colours[2+r][m_a] = Flow::Counter();
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
      pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq'->sqsq', set scale u "<<u<<"\n";
  }
  else {
    /*
      0-----+ +-----2
            | |
	    | |  t
            | |
      1-----+ +-----3

      shower scale is s
    */
    p_colours[0][m_a]   = p_colours[1][m_p]   = Flow::Counter();
    p_colours[2+r][m_a] = p_colours[3-r][m_p] = Flow::Counter();
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
      pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb'->sqsqb', set scale s "<<s<<"\n";
  }
  if (swap) SwapInOrder();
  return 1; 
}

//

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1qbar2_sQ1sQbar2>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsQuark() && flavours[1].IsQuark() &&
	flavours[2].IsSquark() && flavours[3].IsSquark() &&
	flavours[0]!=flavours[1] && 
	((flavours[0].IsAnti() && !flavours[1].IsAnti()) ||
	 (!flavours[0].IsAnti() && flavours[1].IsAnti()) )
	&&
	((SuperPartner(flavours[2],flavours[0])  && 
	  SuperPartner(flavours[3],flavours[1])) ||
	 (SuperPartner(flavours[3],flavours[0])  && 
	  SuperPartner(flavours[2],flavours[1]))) &&
	((abs(flavours[2].Kfcode())>2000000 && abs(flavours[3].Kfcode())>2000000) ||
	 (abs(flavours[2].Kfcode())<2000000 && abs(flavours[3].Kfcode())<2000000))) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1qbar2_sQ1sQbar2(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1qbar2_sQ1sQbar2::XS_q1qbar2_sQ1sQbar2(const size_t nin,const size_t nout, 
			       const ATOOLS::Flavour *fl,
			       XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=fl[1].IsAnti();
  m_g3=std::abs(p_model->Constant("g_3"));
  m_mgluino2=sqr(Flavour(kf_Gluino).Mass());
  m_msq32=sqr(fl[2].Mass());
  m_msq42=sqr(fl[3].Mass());
  m_nstrong=4;
}

double XS_q1qbar2_sQ1sQbar2::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double pT2=(u*t-m_msq32*m_msq42)/s;
  return sqr(m_g3*m_g3)*2.0/9.0*s*pT2/sqr(t-m_mgluino2);
}

bool XS_q1qbar2_sQ1sQbar2::SetColours(double s,double t,double u) 
{ 
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  int r = !SuperPartner(p_flavours[2],p_flavours[0]);
  /*
    0-----+ +-----2
          | |
          | |  t
          | |
    1-----+ +-----3

      shower scale is s
  */
  p_colours[0][m_a]   = p_colours[1][m_p]   = Flow::Counter();
  p_colours[2+r][m_a] = p_colours[3-r][m_p] = Flow::Counter();
  m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
    pow(s*t*u,1.0/3.0);
  msg_Debugging()<<"xs: qqb'->sqsqb', set scale s "<<s<<"\n";
  if (swap) SwapInOrder();
  return 1; 
}

//

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1qbar1_sQ2sQbar2>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model,
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsQuark() && flavours[1]==flavours[0].Bar() &&
	flavours[2].IsSquark() && flavours[3]==flavours[2].Bar() &&
	!SuperPartner(flavours[2],flavours[0])) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1qbar1_sQ2sQbar2(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1qbar1_sQ2sQbar2::XS_q1qbar1_sQ2sQbar2(const size_t nin,const size_t nout, 
				       const ATOOLS::Flavour *fl,
				       XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  m_g3=std::abs(p_model->Constant("g_3"));
  m_msquark2=sqr(p_flavours[2].Mass());
  m_nstrong=4;
}

double XS_q1qbar1_sQ2sQbar2::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double pT2= (u*t-sqr(m_msquark2))/s;
  return sqr(m_g3*m_g3)*4.0/9.0*pT2/s; 
}

bool XS_q1qbar1_sQ2sQbar2::SetColours(double s,double t,double u) 
{ 
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  int r = !(p_flavours[0].IsAnti() == p_flavours[2].IsAnti());
  /*
    0\         /2, if fl[0].IsAnti()==fl[2].IsAnti()
      \   s   /
       ======= 
      /       \
    1/         \3, if fl[0].IsAnti()==fl[2].IsAnti()

    shower scale is t
  */
  p_colours[0][m_a] = p_colours[2+r][m_a] = Flow::Counter();
  p_colours[1][m_p] = p_colours[3-r][m_p] = Flow::Counter();

  m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
  msg_Debugging()<<"xs: qqb->sq'sqb', set scale t "<<t<<"\n";
  if (swap) SwapInOrder();
  return 1; 
}

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1q1_sQ1sQ1>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsQuark() && flavours[2].IsSquark() && 
	flavours[1]==flavours[0] && flavours[2]==flavours[3] &&
	SuperPartner(flavours[2],flavours[0])) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1q1_sQ1sQ1(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1q1_sQ1sQ1::XS_q1q1_sQ1sQ1(const size_t nin,const size_t nout, 
			       const ATOOLS::Flavour *fl,
			       XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_g3=std::abs(p_model->Constant("g_3"));
  m_msquark2=sqr(p_flavours[2].Mass());
  m_mgluino2=sqr(Flavour(kf_Gluino).Mass());
  m_nstrong=4;
}

double XS_q1q1_sQ1sQ1::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  return sqr(m_g3*m_g3)/9.0 * m_mgluino2*s *
    (1./sqr(t-m_mgluino2) + 1./sqr(u-m_mgluino2) - 2./3./(t-m_mgluino2)/(u-m_mgluino2));
}

bool XS_q1q1_sQ1sQ1::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double fac(m_mgluino2*s);
  double Mt(fac/sqr(t-m_mgluino2)); 
  double Mu(fac/sqr(u-m_mgluino2)); 
  if (Mt > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq->sqsq, set scale u "<<u<<"\n";
    /*
      0----\   /----2
            \ /
             X  t
            / \
      1----/   \----3

      shower scale is u
    */
    p_colours[3][m_a] = p_colours[0][m_a] = Flow::Counter();
    p_colours[2][m_a] = p_colours[1][m_a] = Flow::Counter();
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq->sqsq, set scale t "<<t<<"\n";
    /*
      0----\   /----2
            \ /
             =  u
            / \
      1----/   \----3

      shower scale is t
    */
    p_colours[2][m_a] = p_colours[0][m_a] = Flow::Counter();
    p_colours[3][m_a] = p_colours[1][m_a] = Flow::Counter();
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1q1_sQ1LsQ1R>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsQuark() && flavours[2].IsSquark() && 
	flavours[1]==flavours[0] &&
	SuperPartner(flavours[2],flavours[0]) &&
	SuperPartner(flavours[3],flavours[0]) && 
	((abs(flavours[2].Kfcode())<2000000 && abs(flavours[3].Kfcode())>2000000) ||
	 (abs(flavours[2].Kfcode())>2000000 && abs(flavours[3].Kfcode())<2000000))) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1q1_sQ1LsQ1R(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1q1_sQ1LsQ1R::XS_q1q1_sQ1LsQ1R(const size_t nin,const size_t nout, 
			       const ATOOLS::Flavour *fl,
			       XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_g3=std::abs(p_model->Constant("g_3"));
  m_msq32=sqr(p_flavours[2].Mass());
  m_msq42=sqr(p_flavours[3].Mass());
  m_mgluino2=sqr(Flavour(kf_Gluino).Mass());
  m_nstrong=4;
}

double XS_q1q1_sQ1LsQ1R::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double pT2(u*t-m_msq32*m_msq42);
  return sqr(m_g3*m_g3)*2./9. * m_mgluino2*pT2 *
    (1./sqr(t-m_mgluino2) + 1./sqr(u-m_mgluino2));
}

bool XS_q1q1_sQ1LsQ1R::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double pT2(u*t-m_msq32*m_msq42);
  double fac(m_mgluino2*pT2);
  double Mt(fac/sqr(t-m_mgluino2)); 
  double Mu(fac/sqr(u-m_mgluino2)); 
  if (Mt > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq->sqsq, set scale u "<<u<<"\n";
    /*
      0----\   /----2
            \ /
             X  t
            / \
      1----/   \----3

      shower scale is u
    */
    p_colours[3][m_a] = p_colours[0][m_a] = Flow::Counter();
    p_colours[2][m_a] = p_colours[1][m_a] = Flow::Counter();
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qq->sqsq, set scale t "<<t<<"\n";
    /*
      0----\   /----2
            \ /
             =  u
            / \
      1----/   \----3

      shower scale is t
    */
    p_colours[2][m_a] = p_colours[0][m_a] = Flow::Counter();
    p_colours[3][m_a] = p_colours[1][m_a] = Flow::Counter();
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1qbar1_sQ1sQbar1>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsQuark()  && flavours[1]==flavours[0].Bar() &&
	flavours[2].IsSquark() && flavours[2]==flavours[3].Bar() &&
	((SuperPartner(flavours[2],flavours[0])  && 
	  SuperPartner(flavours[3],flavours[1])) ||
	 (SuperPartner(flavours[3],flavours[0])  && 
	  SuperPartner(flavours[2],flavours[1])))) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1qbar1_sQ1sQbar1(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1qbar1_sQ1sQbar1::XS_q1qbar1_sQ1sQbar1(const size_t nin,const size_t nout, 
					   const ATOOLS::Flavour *fl,
					   XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  m_r=!SuperPartner(fl[2],fl[0]);
  m_g3=std::abs(p_model->Constant("g_3"));
  m_mgluino2=sqr(Flavour(kf_Gluino).Mass());
  m_msquark2=sqr(fl[2].Mass());
  m_nstrong=4;
}

double XS_q1qbar1_sQ1sQbar1::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double pT2 = (u*t-sqr(m_msquark2))/s;
  return sqr(m_g3*m_g3)*2.0/9.0*s*pT2*(1./sqr(t-m_mgluino2)+2./sqr(s)-2./3/(s*(t-m_mgluino2)));
}


bool XS_q1qbar1_sQ1sQbar1::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double fac = (u*t-sqr(m_msquark2));
  double Mt(fac/sqr(t-m_mgluino2)); 
  double Ms(fac*2./sqr(s)); 
  if (Ms >  (Mt+Ms) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->sqsqb, set scale t "<<t<<"\n";
    /*
      0\         /2, if fl[0]==fl[2]
        \   s   /
         =======
        /       \
      1/         \3, if fl[0]==fl[2]

      shower scale is t
    */
    p_colours[0][m_a] = p_colours[2+m_r][m_a] = Flow::Counter();	
    p_colours[1][m_p] = p_colours[3-m_r][m_p] = Flow::Counter();
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->sqsqb, set scale s "<<s<<"\n";
    /*
      0----+ +----2
           | |
           | | t
           | |
      1----+ +----3

      shower scale is s
    */
    p_colours[0][m_a]   = p_colours[1][m_p]   = Flow::Counter();	
    p_colours[2+m_r][m_a] = p_colours[3-m_r][m_p] = Flow::Counter();
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1qbar1_GluinoGluino>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsQuark() && flavours[1]==flavours[0].Bar() &&
	flavours[2].IsGluino() && flavours[3].IsGluino()) { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1qbar1_GluinoGluino(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1qbar1_GluinoGluino::XS_q1qbar1_GluinoGluino(const size_t nin,const size_t nout, 
				       const ATOOLS::Flavour *fl,
				       XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  int flav = abs(fl[0].Kfcode());
  m_msqL2=sqr(Flavour((kf_code)(1000000+flav)).Mass());
  m_msqR2=sqr(Flavour((kf_code)(2000000+flav)).Mass());
  m_g3=std::abs(p_model->Constant("g_3"));
  m_nstrong=4;
}

double XS_q1qbar1_GluinoGluino::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double t3(t-m_mgluino2), u4(u-m_mgluino2);
  double deltaL = m_msqL2 - m_mgluino2;
  double deltaR = m_msqR2 - m_mgluino2;
  double pT2 = (u*t-sqr(m_mgluino2))/s;
  double CL = 2.*pT2/s*((sqr(u4)-sqr(deltaL))+(sqr(t3)-sqr(deltaL)) - sqr(s)/9.)/((u4-deltaL)*(t3-deltaL))
    + sqr(deltaL)*(1./sqr(t3-deltaL) + 1./sqr(u4-deltaL) - 1./9.*sqr(1./(t3-deltaL) - 1./(u4-deltaL)));
  double CR = 2.*pT2/s*((sqr(u4)-sqr(deltaR))+(sqr(t3)-sqr(deltaR)) - sqr(s)/9.)/((u4-deltaR)*(t3-deltaR))
    + sqr(deltaR)*(1./sqr(t3-deltaR) + 1./sqr(u4-deltaR) - 1./9.*sqr(1./(t3-deltaR) - 1./(u4-deltaR)));
  return sqr(m_g3*m_g3)/3.*(CL+CR); 
}


bool XS_q1qbar1_GluinoGluino::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double deltaL = m_msqL2 - m_mgluino2;
  double deltaR = m_msqR2 - m_mgluino2;
  double pT2 = (u*t-sqr(m_mgluino2))/s;
  double t3(t-m_mgluino2), u4(u-m_mgluino2);
  double Mt = 2.*pT2/s*(((sqr(u4)-sqr(deltaL)))/((u4-deltaL)*(t3-deltaL)) +
			((sqr(u4)-sqr(deltaR)))/((u4-deltaR)*(t3-deltaR)))
    + sqr(deltaL)/sqr(t3-deltaL) + sqr(deltaR)/sqr(t3-deltaR);
  double Mu = 2.*pT2/s*(((sqr(t3)-sqr(deltaL)))/((u4-deltaL)*(t3-deltaL)) +
			((sqr(t3)-sqr(deltaR)))/((u4-deltaR)*(t3-deltaR)))
    + sqr(deltaL)/sqr(u4-deltaL) + sqr(deltaR)/sqr(u4-deltaR);
  p_colours[0][m_a] = Flow::Counter();
  p_colours[1][m_p] = Flow::Counter();
  if (Mt > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->GluinoGluino, set scale s/t "<<s<<"/"<<t<<"\n";
    /*
      0------+====2
             |
             | t
             |
      1------+====3

      shower scale is s / t
    */
    p_colours[2][m_a] = p_colours[0][m_a];
    p_colours[3][m_p] = p_colours[1][m_p];
    p_colours[2][m_p] = p_colours[3][m_a] = Flow::Counter();
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: qqb->GluinoGluino, set scale s/u "<<s<<"/"<<u<<"\n";
    /*
      0----\ +-==2
            \|/
             | u
            /|\
      1----/ +-==3

      shower scale is s / u
    */
    p_colours[3][m_a] = p_colours[0][m_a];
    p_colours[2][m_p] = p_colours[1][m_p];
    p_colours[3][m_p] = p_colours[2][m_a] = Flow::Counter();
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_gg_sQ1sQbar1>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsGluon() && flavours[1].IsGluon() && 
	flavours[2].IsSquark() && flavours[3]==flavours[2].Bar()) { 
      if (nqcd==2 && nqed==0) {
	return new XS_gg_sQ1sQbar1(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_gg_sQ1sQbar1::XS_gg_sQ1sQbar1(const size_t nin,const size_t nout, 
				 const ATOOLS::Flavour *fl,
				 XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_r=fl[2].IsAnti();
  m_g3=std::abs(p_model->Constant("g_3"));
  m_msquark2=sqr(p_flavours[2].Mass());
  m_nstrong=4;
}

double XS_gg_sQ1sQbar1::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double t3(t-m_msquark2), u4(u-m_msquark2);
  double pT2= (u*t-sqr(m_msquark2))/s;
  return sqr(m_g3*m_g3)*3./16.*(sqr(s*pT2)+sqr(m_msquark2*s))*
    (sqr(u4)+sqr(t3)-sqr(s)/9.)/sqr(s*t3*u4); 
}


bool XS_gg_sQ1sQbar1::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  double t3(t-m_msquark2), u4(u-m_msquark2);
  double pT2= (u*t-sqr(m_msquark2))/s;
  double fac((sqr(s*pT2)+sqr(m_msquark2*s))/sqr(s*t3*u4));
  double Mt(fac*sqr(u4));
  double Mu(fac*sqr(t3));
  p_colours[0][0] = Flow::Counter();
  p_colours[0][1] = Flow::Counter();
  if (Mt*(1-m_r) +Mu*m_r > (Mt+Mu) * ran.Get()) {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: gg->sqsqb, set scale t/s "<<t<<"/"<<s<<"\n";
    /*
      0====+------2
           |
           |  t
           |
      1====+------3

      shower scale is s / t
    */
    p_colours[2+m_r][0] = p_colours[0][0];
    p_colours[3-m_r][1] = p_colours[1][1] = Flow::Counter();
    p_colours[1][0] = p_colours[0][1];
  }
  else {
    m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
    msg_Debugging()<<"xs: gg->sqsqb, set scale u/s "<<u<<"/"<<s<<"\n";
    /*
      0==-+ /----2
         \|/
          |  u
         /|\
      1==-+ \----3

      shower scale is u / s
    */
    p_colours[2+m_r][0] = p_colours[1][0] = Flow::Counter();
    p_colours[3-m_r][1] = p_colours[0][1];
    p_colours[1][1] = p_colours[0][0];
  }
  if (swap) SwapInOrder();
  return true;
}

namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_q1g_sQ1Gluino>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (((flavours[0].IsQuark() && flavours[1].IsGluon()) && 
	 ((SuperPartner(flavours[2],flavours[0]) && flavours[3].IsGluino()) || 
	  (SuperPartner(flavours[3],flavours[0]) && flavours[2].IsGluino()))) ||
	((flavours[1].IsQuark() && flavours[0].IsGluon()) && 
	 ((SuperPartner(flavours[2],flavours[1]) && flavours[3].IsGluino()) || 
	  (SuperPartner(flavours[3],flavours[1]) && flavours[2].IsGluino()))))  { 
      if (nqcd==2 && nqed==0) {
	return new XS_q1g_sQ1Gluino(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_q1g_sQ1Gluino::XS_q1g_sQ1Gluino(const size_t nin,const size_t nout, 
				   const ATOOLS::Flavour *fl,
				   XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_iniq=0;
  m_swaput=0;
  if (fl[1].IsQuark()){
    m_iniq=1;
    m_swaput=1;
  }
  m_finq=2;
  if (fl[3].IsSquark()) {
    m_finq=3;
    if (m_swaput) m_swaput=0;
    else m_swaput=1;
  }
  m_a=fl[m_iniq].IsAnti();
  m_p=1-m_a;
  m_msquark2=sqr(p_flavours[m_finq].Mass());
  m_mgluino2=sqr(Flavour(kf_Gluino).Mass());
  m_g3=std::abs(p_model->Constant("g_3"));
  m_nstrong=4;
}

double XS_q1g_sQ1Gluino::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  if (m_swaput) std::swap<double>(t,u);
  
  double u4(u-m_msquark2), t3(t-m_mgluino2);
  return sqr(m_g3*m_g3)/4.*(-u4-2.*(m_msquark2-m_mgluino2)*(1.+m_mgluino2/t3+m_msquark2/u4))*
    (sqr(u4)+sqr(s)-sqr(t3/3.))/(s*t3*u4); 
}

bool XS_q1g_sQ1Gluino::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  if (m_swaput) std::swap<double>(t,u);
  
  double u4(u-m_msquark2), t3(t-m_mgluino2);
  double fac((-u4-2.*(m_msquark2-m_mgluino2)*(1.+m_mgluino2/t3+m_msquark2/u4)));
  double Ms(fac*sqr(u4)/(s*t3*u4));
  double Mu(fac*sqr(s)/(s*t3*u4));
  p_colours[m_iniq][m_a] = Flow::Counter();
  p_colours[m_finq][m_a] = Flow::Counter();
  if (Mu > (Ms+Mu) * ran.Get()) {
    /*
      1====+----2, if fl[2].IsSquark() 
           |
           |  u
           |
      0----+====3, if fl[0].IsQuark()

      shower scale is t/u
    */
    p_colours[5-m_finq][m_a] = p_colours[m_iniq][m_a];
    p_colours[5-m_finq][m_p] = p_colours[1-m_iniq][m_p] = Flow::Counter();
    p_colours[1-m_iniq][m_a] = p_colours[m_finq][m_a];
    if (dabs(t)>dabs(u)) {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
	pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->sqGluino, set scale t "<<t<<"\n";
    }
    else {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
	pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->sqGluino, set scale u "<<u<<"\n";
    }
  }
  else {
    /*
      0\        /2, if fl[0].IsQuark && fl[2].IsSquark() 
        \      /
        | +--+ | 
        //    \\
      1//      \\3

      shower scale is s/t
    */
    p_colours[5-m_finq][m_p] = p_colours[m_finq][m_a];
    p_colours[1-m_iniq][m_a] = p_colours[5-m_finq][m_a] = Flow::Counter();
    p_colours[1-m_iniq][m_p] = p_colours[m_iniq][m_a];
    if (dabs(t)>s) {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
	pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->sqGluino, set scale t "<<t<<"\n";
    }
    else {
      m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = 
	pow(s*t*u,1.0/3.0);
      msg_Debugging()<<"xs: qg->sqGluino, set scale s "<<s<<"\n";
    }
  }
  if (swap) SwapInOrder();
  return true;
}
      
namespace EXTRAXS {

  template <> 
  Single_Process *Single_Process::GetProcess<XS_gg_GluinoGluino>
  (const size_t nin,const size_t nout,
   const ATOOLS::Flavour *flavours,XS_Model_Base *const model, 
   const size_t nqed, const size_t nqcd)
  {
    if (!model->IncludesModel("MSSM")) return NULL;
    if (flavours[0].IsGluon()  && flavours[1].IsGluon() &&
	flavours[2].IsGluino() && flavours[3].IsGluino()) { 
      if (nqcd==2 && nqed==0) {
	return new XS_gg_GluinoGluino(nin,nout,flavours,model); 
      }
    }
    return NULL;
  }
}

XS_gg_GluinoGluino::XS_gg_GluinoGluino(const size_t nin,const size_t nout, 
				 const ATOOLS::Flavour *fl,
				 XS_Model_Base *const model): 
  Single_Process(nin,nout,fl,model) 
{
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_g3=std::abs(p_model->Constant("g_3"));
  m_mgluino2=sqr(Flavour(kf_Gluino).Mass());
  m_nstrong=4;
}

double XS_gg_GluinoGluino::operator()(double s,double t,double u) 
{
  if (s<m_threshold) return 0.;
  double u4(u-m_mgluino2), t3(t-m_mgluino2);
  double pT2= (u*t-sqr(m_mgluino2))/s;
  return sqr(m_g3*m_g3)*9.0/8.0*u4*t3/2.*
    (sqr(u4)+sqr(t3)+4.*m_mgluino2*sqr(s)*pT2/(u4*t3))*
    (1./(sqr(s)*sqr(t3)) + 1./(sqr(s)*sqr(u4)) + 1./(sqr(u4)*sqr(t3)));
}
  
bool XS_gg_GluinoGluino::SetColours(double s, double t, double u) 
{
  bool swap=m_swaped;
  RestoreInOrder();
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    -1.0/(1.0/s+1.0/t+1.0/u);
  m_scale[PHASIC::stp::sfs] = m_scale[PHASIC::stp::sis] = pow(s*t*u,1.0/3.0);
  msg_Debugging()<<"xs: gg->GluinoGluino, set scale s "<<s<<"\n";
  p_colours[0][0] = Flow::Counter();
  p_colours[1][1] = Flow::Counter();
  
  double u4(u-m_mgluino2);
  double t3(t-m_mgluino2);
  double pT2= (u*t-sqr(m_mgluino2))/s;
  
  double fac(u4*t3*(sqr(u4)+sqr(t3)+4.*m_mgluino2*sqr(s)*pT2/(u4*t3)));
  double Mst(fac/(sqr(s)*sqr(t3)));
  double Msu(fac/(sqr(s)*sqr(u4)));
  double Mut(fac/(sqr(u4)*sqr(t3)));
  double rr = ran.Get() * (Mst+Msu+Mut);
  if (rr-Mst < 0.) {
    /*
      0====++====2
           ||
           ||  t
           ||
      1====++====3

      shower scale is s/t/u
    */
    p_colours[2][0] = p_colours[0][0];
    p_colours[3][1] = p_colours[1][1];
    p_colours[0][1] = p_colours[1][0] = Flow::Counter();
    p_colours[2][1] = p_colours[3][0] = Flow::Counter();
  }
  else {
    if (rr-Mst-Msu < 0.) {
      /*
	0====+\---==3
             ||\ /
             || X u
             ||/ \
	1====+/---==2
	   
	shower scale is s/t/u
      */
      p_colours[3][0] = p_colours[0][0];
      p_colours[2][1] = p_colours[1][1];
      p_colours[0][1] = p_colours[1][0] = Flow::Counter();
      p_colours[3][1] = p_colours[2][0] = Flow::Counter();
    }
    else {
      /*
	0\\       //3
          \\  s  //
          | ===== |
          //     \\
	1//       \\2
	   
	shower scale is s/t/u
      */
      p_colours[2][0] = p_colours[0][0];
      p_colours[3][1] = p_colours[0][1] = Flow::Counter();
      p_colours[2][1] = p_colours[1][1];
      p_colours[3][0] = p_colours[1][0] = Flow::Counter();
    }
  }
  if (swap) SwapInOrder();
  return true;
}

