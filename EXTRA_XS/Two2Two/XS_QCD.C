#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flow.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Process/Process_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"

#define PropID(i,j) ((1<<i)|(1<<j))

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;


namespace EXTRAXS {
  class XS_q1q2_q1q2 : public ME2_Base {
  private:

    int    m_a, m_p, m_r;
    double m_m12, m_m22, m_g;
    
  public:

    XS_q1q2_q1q2(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}

DECLARE_TREEME2_GETTER(XS_q1q2_q1q2,"1XS_q1q2_q1q2")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_q1q2_q1q2>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo ||
      pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsQuark() && fl[1].IsQuark() &&
      fl[0]!=fl[1] &&
      ((fl[2]==fl[0] && fl[3]==fl[1]) ||
       (fl[3]==fl[0] && fl[2]==fl[1]))) {
    if (pi.m_maxcpl[0]==2 && pi.m_maxcpl[1]==0 &&
	pi.m_mincpl[0]==2 && pi.m_mincpl[1]==0) {
      return new XS_q1q2_q1q2(pi,fl);
    }
  }
  return NULL;
}

XS_q1q2_q1q2::XS_q1q2_q1q2(const Process_Info& pi, const Flavour_Vector& fl):
  ME2_Base(pi,fl) 
{
  DEBUG_FUNC(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_r=!(fl[0] == fl[2]);
  m_a=fl[0].IsAnti();
  m_p=fl[1].IsAnti();
  m_g=sqrt(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_S"));
  m_m12=sqr(m_flavs[0].Mass());
  m_m22=sqr(m_flavs[1].Mass());
  m_oew=0; m_oqcd=2;
  if (!m_r) {
    m_cfls[PropID(0,2)].push_back(kf_gluon);
    m_cfls[PropID(1,3)].push_back(kf_gluon);
  }
  else {
    m_cfls[PropID(1,2)].push_back(kf_gluon);
    m_cfls[PropID(0,3)].push_back(kf_gluon);
  }
}

double XS_q1q2_q1q2::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2+m_r]).Abs2();
  double u=(mom[0]-mom[3-m_r]).Abs2();
  //if (s<m_threshold) return 0.;
  double Mt(sqr(u-m_m12-m_m22)+sqr(s-m_m12-m_m22)+2.0*t*(m_m12+m_m22));
  return sqr(m_g*m_g)*CouplingFactor(2,0)*4.0/9.0*Mt/(t*t);
}

bool XS_q1q2_q1q2::SetColours(const Vec4D_Vector& mom) 
{ 
  if (m_a==m_p) {
    /*
    
      0-----\   /-----2, if fl[0]==fl[2]
             \ /
	      X  t
             / \
      1-----/   \-----3, if fl[1]==fl[3]

    */
    p_colours[0][m_a] = p_colours[3-m_r][m_a] = Flow::Counter();
    p_colours[1][m_a] = p_colours[2+m_r][m_a] = Flow::Counter();
    //msg_Debugging()<<"xs: qq'->qq', set scale u "<<u<<"\n";
  }
  else {
    /*
    
      0-----+ +-----2
            | |
	    | |  t
            | |
      1-----+ +-----3

    */
    p_colours[0][m_a]   = p_colours[1][m_p]   = Flow::Counter();
    p_colours[2+m_r][m_a] = p_colours[3-m_r][m_p] = Flow::Counter();
    //msg_Debugging()<<"xs: qqb'->qqb', set scale s "<<s<<"\n";
  }
  return 1; 
}


namespace EXTRAXS {
  class XS_q1qbar1_q2qbar2 : public ME2_Base {
  private:

    int    m_a, m_p, m_r;
    double m_m12, m_m32, m_g;

  public:

    XS_q1qbar1_q2qbar2(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}

DECLARE_TREEME2_GETTER(XS_q1qbar1_q2qbar2,"1XS_q1qbar1_q2qbar2")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_q1qbar1_q2qbar2>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo ||
      pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsQuark() && fl[1]==fl[0].Bar() &&
      fl[2].IsQuark() && fl[3]==fl[2].Bar() &&
      fl[0]!=fl[2] && fl[0]!=fl[3]) {
    if (pi.m_maxcpl[0]==2 && pi.m_maxcpl[1]==0 &&
	pi.m_mincpl[0]==2 && pi.m_mincpl[1]==0) {
      return new XS_q1qbar1_q2qbar2(pi,fl); 
    }
  }
  return NULL;
}

XS_q1qbar1_q2qbar2::XS_q1qbar1_q2qbar2(const Process_Info& pi,
                                       const Flavour_Vector& fl):
  ME2_Base(pi,fl)
{
  DEBUG_FUNC(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_r=!(fl[0].IsAnti()==fl[2].IsAnti());
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  m_g=sqrt(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_S"));
  m_m12=sqr(m_flavs[0].Mass());
  m_m32=sqr(m_flavs[2].Mass());
  m_oew=0; m_oqcd=2;
  m_cfls[PropID(0,1)].push_back(kf_gluon);
  m_cfls[PropID(2,3)].push_back(kf_gluon);
}

double XS_q1qbar1_q2qbar2::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  //if (s<m_threshold) return 0.;
  double Ms(sqr(u-m_m12-m_m32)+sqr(t-m_m12-m_m32)+2.0*s*(m_m12+m_m32));
  return sqr(m_g*m_g)*CouplingFactor(2,0)*4.0/9.0*Ms/(s*s); 
}

bool XS_q1qbar1_q2qbar2::SetColours(const Vec4D_Vector& mom) 
{ 
  /*
  
    0\         /2, if fl[0].IsAnti()==fl[2].IsAnti()
    \   s   /
    ======= 
    /       \
    1/         \3, if fl[0].IsAnti()==fl[2].IsAnti()

  */
  p_colours[0][m_a] = p_colours[2+m_r][m_a] = Flow::Counter();
  p_colours[1][m_p] = p_colours[3-m_r][m_p] = Flow::Counter();

  //msg_Debugging()<<"xs: qqb->q'qb', set scale t "<<t<<"\n";
  return 1; 
}


namespace EXTRAXS {
  class XS_q1q1_q1q1 : public ME2_Base {
  private:

    int    m_a;
    double m_m12, m_g;

  public:

    XS_q1q1_q1q1(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}

DECLARE_TREEME2_GETTER(XS_q1q1_q1q1,"1XS_q1q1_q1q1")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_q1q1_q1q1>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo ||
      pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsQuark() && fl[1]==fl[0] &&
      fl[2]==fl[0] && fl[3]==fl[0]) { 
    if (pi.m_maxcpl[0]==2 && pi.m_maxcpl[1]==0 &&
	pi.m_mincpl[0]==2 && pi.m_mincpl[1]==0) {
      return new XS_q1q1_q1q1(pi,fl); 
    }
  }
  return NULL;
}

XS_q1q1_q1q1::XS_q1q1_q1q1(const Process_Info& pi, const Flavour_Vector& fl): 
  ME2_Base(pi,fl) 
{
  DEBUG_FUNC(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_g=sqrt(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_S"));
  m_m12=sqr(m_flavs[0].Mass());
  m_oew=0; m_oqcd=2;
  m_cfls[PropID(0,2)].push_back(kf_gluon);
  m_cfls[PropID(0,3)].push_back(kf_gluon);
  m_cfls[PropID(1,2)].push_back(kf_gluon);
  m_cfls[PropID(1,3)].push_back(kf_gluon);
}

double XS_q1q1_q1q1::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  //if (s<m_threshold) return 0.;
  double Mt(sqr(u-2.0*m_m12)+sqr(s-2.0*m_m12)+4.0*t*m_m12); 
  double Mu(sqr(t-2.0*m_m12)+sqr(s-2.0*m_m12)+4.0*u*m_m12); 
  double Mtu(-2.0/3.0*(s*s-8.0*(s-2.0*m_m12)*m_m12-4.0*m_m12*m_m12));
  return sqr(m_g*m_g)*CouplingFactor(2,0)*4.0/9.0
         *(Mt/(t*t)+Mu/(u*u)+Mtu/(t*u))/2.0;
}

bool XS_q1q1_q1q1::SetColours(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  double Mt(sqr(u-2.0*m_m12)+sqr(s-2.0*m_m12)+4.0*t*m_m12); 
  double Mu(sqr(t-2.0*m_m12)+sqr(s-2.0*m_m12)+4.0*u*m_m12); 
  if (Mt > (Mt+Mu) * ran->Get()) {
    msg_Debugging()<<"xs: qq->qq, set scale u "<<u<<"\n";
    /*
    
      0----\   /----2
            \ /
             X  t
            / \
      1----/   \----3

    */
    p_colours[3][m_a] = p_colours[0][m_a] = Flow::Counter();
    p_colours[2][m_a] = p_colours[1][m_a] = Flow::Counter();
  }
  else {
    msg_Debugging()<<"xs: qq->qq, set scale t "<<t<<"\n";
    /*

      0----\   /----2
            \ /
             =  u
            / \
      1----/   \----3

    */
    p_colours[2][m_a] = p_colours[0][m_a] = Flow::Counter();
    p_colours[3][m_a] = p_colours[1][m_a] = Flow::Counter();
  }
  return true;
}


namespace EXTRAXS {
  class XS_q1qbar1_q1qbar1 : public ME2_Base {
  private:

    int    m_a, m_p, m_r;
    double m_m12, m_g;

  public:

    XS_q1qbar1_q1qbar1(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };
}

DECLARE_TREEME2_GETTER(XS_q1qbar1_q1qbar1,"1XS_q1qbar1_q1qbar1")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_q1qbar1_q1qbar1>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo ||
      pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsQuark() && fl[1]==fl[0].Bar() &&
      ((fl[2]==fl[0] && fl[3]==fl[1]) ||
       (fl[3]==fl[0] && fl[2]==fl[1]))) { 
    if (pi.m_maxcpl[0]==2 && pi.m_maxcpl[1]==0 &&
	pi.m_mincpl[0]==2 && pi.m_mincpl[1]==0) {
      return new XS_q1qbar1_q1qbar1(pi,fl); 
    }
  }
  return NULL;
}

XS_q1qbar1_q1qbar1::XS_q1qbar1_q1qbar1(const Process_Info& pi,
                                       const Flavour_Vector& fl):
  ME2_Base(pi,fl) 
{
  DEBUG_FUNC(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  m_r=(fl[0]!=fl[2]);
  m_g=sqrt(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_S"));
  m_m12=sqr(m_flavs[0].Mass());
  m_oew=0; m_oqcd=2;
  m_cfls[PropID(0,1)].push_back(kf_gluon);
  m_cfls[PropID(2,3)].push_back(kf_gluon);
  if (!m_r) {
    m_cfls[PropID(0,2)].push_back(kf_gluon);
    m_cfls[PropID(1,3)].push_back(kf_gluon);
  }
  else {
    m_cfls[PropID(1,2)].push_back(kf_gluon);
    m_cfls[PropID(0,3)].push_back(kf_gluon);
  }
}

double XS_q1qbar1_q1qbar1::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2+m_r]).Abs2();
  double u=(mom[0]-mom[3-m_r]).Abs2();
  //if (s<m_threshold) return 0.;
  double Mt(sqr(s-2.0*m_m12)+sqr(u-2.0*m_m12)+4.0*t*m_m12); 
  double Ms(sqr(t-2.0*m_m12)+sqr(u-2.0*m_m12)+4.0*s*m_m12); 
  double Mts(-2.0/3.0*(u*u-8.0*(u-2.0*m_m12)*m_m12-4.0*m_m12*m_m12));
  return sqr(m_g*m_g)*CouplingFactor(2,0)*4.0/9.0*(Mt/(t*t)+Ms/(s*s)+Mts/(t*s));
}


bool XS_q1qbar1_q1qbar1::SetColours(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  double Mt(sqr(s-2.0*m_m12)+sqr(u-2.0*m_m12)+4.0*t*m_m12); 
  double Ms(sqr(t-2.0*m_m12)+sqr(u-2.0*m_m12)+4.0*s*m_m12); 
  if (Ms >  (Mt+Ms) * ran->Get()) {
    msg_Debugging()<<"xs: qqb->qqb, set scale t "<<t<<"\n";
    /*
    
      0\         /2, if fl[0]==fl[2]
        \   s   /
         =======
        /       \
      1/         \3, if fl[0]==fl[2]

    */
    p_colours[0][m_a] = p_colours[2+m_r][m_a] = Flow::Counter();	
    p_colours[1][m_p] = p_colours[3-m_r][m_p] = Flow::Counter();
  }
  else {
    msg_Debugging()<<"xs: qqb->qqb, set scale s "<<s<<"\n";
    /*

      0----+ +----2
           | |
           | | t
           | |
      1----+ +----3

    */
    p_colours[0][m_a]   = p_colours[1][m_p]   = Flow::Counter();	
    p_colours[2+m_r][m_a] = p_colours[3-m_r][m_p] = Flow::Counter();
  }
  return true;
}


namespace EXTRAXS {
  class XS_q1qbar1_gg : public ME2_Base {
  private:

    int    m_a, m_p;
    double m_m12, m_g;

  public:

    XS_q1qbar1_gg(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
    bool SetColours();

  };
}

DECLARE_TREEME2_GETTER(XS_q1qbar1_gg,"1XS_q1qbar1_gg")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_q1qbar1_gg>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo ||
      pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsQuark() && fl[1]==fl[0].Bar() &&
      fl[2].IsGluon() && fl[3].IsGluon()) { 
    if (pi.m_maxcpl[0]==2 && pi.m_maxcpl[1]==0 &&
	pi.m_mincpl[0]==2 && pi.m_mincpl[1]==0) {
      return new XS_q1qbar1_gg(pi,fl); 
    }
  }
  return NULL;
}

XS_q1qbar1_gg::XS_q1qbar1_gg(const Process_Info& pi, const Flavour_Vector& fl): 
  ME2_Base(pi,fl) 
{
  DEBUG_FUNC(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_a=fl[0].IsAnti();
  m_p=1-m_a;
  m_m12=sqr(m_flavs[0].Mass());
  m_g=sqrt(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_S"));
  m_oew=0; m_oqcd=2;
  m_cfls[PropID(0,1)].push_back(kf_gluon);
  m_cfls[PropID(2,3)].push_back(kf_gluon);
  m_cfls[PropID(0,2)].push_back(fl[0].Bar());
  m_cfls[PropID(0,3)].push_back(fl[0].Bar());
  m_cfls[PropID(1,2)].push_back(fl[1].Bar());
  m_cfls[PropID(1,3)].push_back(fl[1].Bar());
}

double XS_q1qbar1_gg::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  //if (s<m_threshold) return 0.;
  double tp(t-m_m12), up(u-m_m12);
  double Mt(32.0/27.0*(tp*up-m_m12*(4.0*(m_m12+tp)+m_m12*tp/s))/(tp*tp));
  double Mu(32.0/27.0*(up*tp-m_m12*(4.0*(m_m12+up)+m_m12*up/s))/(up*up));
  double Ms(-8.0/3.0*(tp*tp+up*up+4.0*m_m12*s)/(s*s));
  return sqr(m_g*m_g)*CouplingFactor(2,0)*(Mt+Mu+Ms)/2.0; 
}


bool XS_q1qbar1_gg::SetColours(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  double tp(t-m_m12), up(u-m_m12);
  double Mt(32.0/27.0*(tp*up-m_m12*(4.0*(m_m12+tp)+m_m12*tp/s))/(tp*tp));
  double Mu(32.0/27.0*(up*tp-m_m12*(4.0*(m_m12+up)+m_m12*up/s))/(up*up));
  p_colours[0][m_a] = Flow::Counter();
  p_colours[1][m_p] = Flow::Counter();
  if (Mt > (Mt+Mu) * ran->Get()) {
    msg_Debugging()<<"xs: qqb->gg, set scale s/t "<<s<<"/"<<t<<"\n";
    /*
    
      0------+====2
             |
             | t
             |
      1------+====3

    */
    p_colours[2][m_a] = p_colours[0][m_a];
    p_colours[3][m_p] = p_colours[1][m_p];
    p_colours[2][m_p] = p_colours[3][m_a] = Flow::Counter();
  }
  else {
    msg_Debugging()<<"xs: qqb->gg, set scale s/u "<<s<<"/"<<u<<"\n";
    /*

      0----\ +-==2
            \|/
             | u
            /|\
      1----/ +-==3

    */
    p_colours[3][m_a] = p_colours[0][m_a];
    p_colours[2][m_p] = p_colours[1][m_p];
    p_colours[3][m_p] = p_colours[2][m_a] = Flow::Counter();
  }
  return true;
}


namespace EXTRAXS {
  class XS_gg_q1qbar1 : public ME2_Base {
  private:

    int    m_r;
    double m_m32, m_g;

  public:

    XS_gg_q1qbar1(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}

DECLARE_TREEME2_GETTER(XS_gg_q1qbar1,"1XS_gg_q1qbar1")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_gg_q1qbar1>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo ||
      pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsGluon() && fl[1].IsGluon() && 
      fl[2].IsQuark() && fl[3]==fl[2].Bar()) { 
    if (pi.m_maxcpl[0]==2 && pi.m_maxcpl[1]==0 &&
	pi.m_mincpl[0]==2 && pi.m_mincpl[1]==0) {
      return new XS_gg_q1qbar1(pi,fl); 
    }
  }
  return NULL;
}

XS_gg_q1qbar1::XS_gg_q1qbar1(const Process_Info& pi, const Flavour_Vector& fl): 
  ME2_Base(pi,fl) 
{
  DEBUG_FUNC(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_r=fl[2].IsAnti();
  m_g=sqrt(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_S"));
  m_m32=sqr(m_flavs[2].Mass());
  m_oew=0; m_oqcd=2;
  m_cfls[PropID(0,1)].push_back(kf_gluon);
  m_cfls[PropID(2,3)].push_back(kf_gluon);
  m_cfls[PropID(0,2)].push_back(fl[2]);
  m_cfls[PropID(1,2)].push_back(fl[2]);
  m_cfls[PropID(0,3)].push_back(fl[3]);
  m_cfls[PropID(1,3)].push_back(fl[3]);
}

double XS_gg_q1qbar1::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  //if (s<m_threshold) return 0.;
  double tp(t-m_m32), up(u-m_m32);
  double Mt(1.0/6.0*(tp*up-m_m32*(4.0*(m_m32+tp)+m_m32*tp/s))/(tp*tp));
  double Mu(1.0/6.0*(up*tp-m_m32*(4.0*(m_m32+up)+m_m32*up/s))/(up*up));
  double Ms(-3.0/8.0*(tp*tp+up*up+4.0*m_m32*s)/(s*s));
  return sqr(m_g*m_g)*CouplingFactor(2,0)*(Mt+Mu+Ms); 
}


bool XS_gg_q1qbar1::SetColours(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  double tp(t-m_m32), up(u-m_m32);
  double Mt(1.0/6.0*(tp*up-m_m32*(4.0*(m_m32+tp)+m_m32*tp/s))/(tp*tp));
  double Mu(1.0/6.0*(up*tp-m_m32*(4.0*(m_m32+up)+m_m32*up/s))/(up*up));
  p_colours[0][0] = Flow::Counter();
  p_colours[0][1] = Flow::Counter();
  if (Mt*(1-m_r) +Mu*m_r > (Mt+Mu) * ran->Get()) {
    msg_Debugging()<<"xs: gg->qqb, set scale t/s "<<t<<"/"<<s<<"\n";
    /*
    
      0====+------2
           |
           |  t
           |
      1====+------3

    */
    p_colours[2+m_r][0] = p_colours[0][0];
    p_colours[3-m_r][1] = p_colours[1][1] = Flow::Counter();
    p_colours[1][0] = p_colours[0][1];
  }
  else {
    msg_Debugging()<<"xs: gg->qqb, set scale u/s "<<u<<"/"<<s<<"\n";
    /*

      0==-+ /----2
         \|/
          |  u
         /|\
      1==-+ \----3

    */
    p_colours[2+m_r][0] = p_colours[1][0] = Flow::Counter();
    p_colours[3-m_r][1] = p_colours[0][1];
    p_colours[1][1] = p_colours[0][0];
  }
  return true;
}


namespace EXTRAXS {
  class XS_q1g_q1g : public ME2_Base {
  private:
  
    int    m_a, m_p, m_finq, m_iniq, m_swaput;
    double m_mq2, m_g;

  public:

    XS_q1g_q1g(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}

DECLARE_TREEME2_GETTER(XS_q1g_q1g,"1XS_q1g_q1g")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_q1g_q1g>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo ||
      pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (((fl[0].IsQuark() && fl[1].IsGluon()) && 
       ((fl[2]==fl[0] && fl[3].IsGluon()) || 
        (fl[3]==fl[0] && fl[2].IsGluon()))) ||
      ((fl[1].IsQuark() && fl[0].IsGluon()) && 
       ((fl[2]==fl[1] && fl[3].IsGluon()) || 
        (fl[3]==fl[1] && fl[2].IsGluon()))))  { 
    if (pi.m_maxcpl[0]==2 && pi.m_maxcpl[1]==0 &&
	pi.m_mincpl[0]==2 && pi.m_mincpl[1]==0) {
      return new XS_q1g_q1g(pi,fl); 
    }
  }
  return NULL;
}

XS_q1g_q1g::XS_q1g_q1g(const Process_Info& pi, const Flavour_Vector& fl): 
  ME2_Base(pi,fl) 
{
  DEBUG_FUNC(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_iniq=0;
  m_swaput=0;
  if (fl[1].IsQuark()){
    m_iniq=1;
    m_swaput=1;
  }
  m_finq=2;
  if (fl[3].IsQuark()) {
    m_finq=3;
    if (m_swaput) m_swaput=0;
    else m_swaput=1;
  }
  m_a=fl[m_iniq].IsAnti();
  m_p=1-m_a;
  m_mq2=sqr(m_flavs[m_iniq].Mass());
  m_g=sqrt(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_S"));
  m_oew=0; m_oqcd=2;
  m_cfls[PropID(0,1)].push_back(fl[m_iniq].Bar());
  m_cfls[PropID(2,3)].push_back(fl[m_finq]);
  m_cfls[PropID((1-m_iniq),m_finq)].push_back(fl[m_finq]);
  m_cfls[PropID((5-m_finq),m_finq)].push_back(fl[m_finq]);
  m_cfls[PropID(m_iniq,(5-m_finq))].push_back(fl[m_iniq].Bar());
  m_cfls[PropID(m_iniq,(1-m_iniq))].push_back(fl[m_iniq].Bar());
}

double XS_q1g_q1g::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  //if (s<m_threshold) return 0.;
  if (m_swaput) std::swap<double>(t,u);
  double sp(s-m_mq2), up(u-m_mq2);
  double Ms(4.0/9.0*(sp*up-m_mq2*(4.0*(m_mq2+sp)+m_mq2*sp/t))/(sp*sp));
  double Mu(4.0/9.0*(up*sp-m_mq2*(4.0*(m_mq2+up)+m_mq2*up/t))/(up*up));
  double Mt(-(sp*sp+up*up+4.0*m_mq2*t)/(t*t));
  return -sqr(m_g*m_g)*CouplingFactor(2,0)*(Ms+Mu+Mt); 
}

bool XS_q1g_q1g::SetColours(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  if (m_swaput) std::swap<double>(t,u);
  double sp(s-m_mq2), up(u-m_mq2);
  double Ms(4.0/9.0*(sp*up-m_mq2*(4.0*(m_mq2+sp)+m_mq2*sp/t))/(sp*sp));
  double Mu(4.0/9.0*(up*sp-m_mq2*(4.0*(m_mq2+up)+m_mq2*up/t))/(up*up));
  p_colours[m_iniq][m_a] = Flow::Counter();
  p_colours[m_finq][m_a] = Flow::Counter();
  if (Mu > (Ms+Mu) * ran->Get()) {
    /*
    
      1====+----2, if fl[2].IsQuark() 
           |
           |  u
           |
      0----+====3, if fl[0].IsQuark()

    */
    p_colours[5-m_finq][m_a] = p_colours[m_iniq][m_a];
    p_colours[5-m_finq][m_p] = p_colours[1-m_iniq][m_p] = Flow::Counter();
    p_colours[1-m_iniq][m_a] = p_colours[m_finq][m_a];
    if (dabs(t)>dabs(u)) {
      msg_Debugging()<<"xs: qg->qg, set scale t "<<t<<"\n";
    }
    else {
      msg_Debugging()<<"xs: qg->qg, set scale u "<<u<<"\n";
    }
  }
  else {
    /*
    
      0\        /2, if fl[0].IsQuark && fl[2].IsQuark() 
        \      /
        | +--+ | 
        //    \\
      1//      \\3

    */
    p_colours[5-m_finq][m_p] = p_colours[m_finq][m_a];
    p_colours[1-m_iniq][m_a] = p_colours[5-m_finq][m_a] = Flow::Counter();
    p_colours[1-m_iniq][m_p] = p_colours[m_iniq][m_a];
    if (dabs(t)>s) {
      msg_Debugging()<<"xs: qg->qg, set scale t "<<t<<"\n";
    }
    else {
      msg_Debugging()<<"xs: qg->qg, set scale s "<<s<<"\n";
    }
  }
  return true;
}


namespace EXTRAXS {
  class XS_gg_gg : public ME2_Base {
  private:
  
    double m_g;

  public:

    XS_gg_gg(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}

DECLARE_TREEME2_GETTER(XS_gg_gg,"1XS_gg_gg")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,XS_gg_gg>::
operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype!=nlo_type::lo ||
      pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].IsGluon() && fl[1].IsGluon() &&
      fl[2].IsGluon() && fl[3].IsGluon()) { 
    if (pi.m_maxcpl[0]==2 && pi.m_maxcpl[1]==0 &&
	pi.m_mincpl[0]==2 && pi.m_mincpl[1]==0) {
      return new XS_gg_gg(pi,fl); 
    }
  }
  return NULL;
}

XS_gg_gg::XS_gg_gg(const Process_Info& pi, const Flavour_Vector& fl): 
  ME2_Base(pi,fl) 
{
  DEBUG_FUNC(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  m_g=sqrt(4.*M_PI*MODEL::s_model->ScalarConstant("alpha_S"));
  m_oew=0; m_oqcd=2;
  m_cfls[PropID(0,1)].push_back(kf_gluon);
  m_cfls[PropID(0,2)].push_back(kf_gluon);
  m_cfls[PropID(0,3)].push_back(kf_gluon);
  m_cfls[PropID(1,2)].push_back(kf_gluon);
  m_cfls[PropID(1,3)].push_back(kf_gluon);
  m_cfls[PropID(2,3)].push_back(kf_gluon);
}

double XS_gg_gg::operator()(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  //if (s<m_threshold) return 0.;
  double Ms(1.0-t*u/(s*s));
  double Mt(1.0-s*u/(t*t));
  double Mu(1.0-s*t/(u*u));
  return sqr(m_g*m_g)*CouplingFactor(2,0)*9.0/2.0*(Ms+Mt+Mu)/2.0;
}
  
bool XS_gg_gg::SetColours(const Vec4D_Vector& mom) 
{
  double s=(mom[0]+mom[1]).Abs2();
  double t=(mom[0]-mom[2]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();
  p_colours[0][0] = Flow::Counter();
  p_colours[1][1] = Flow::Counter();
  double Mu(1.0+t*t/(u*s)-s*t/(u*u)-t*u/(s*s));
  double Ms(1.0+s*s/(t*u)-s*t/(u*u)-u*s/(t*t));
  double Mt(1.0+u*u/(s*t)-u*s/(t*t)-t*u/(s*s));
  double rr = ran->Get() * (Ms+Mt+Mu);
  if (rr-Mt < 0.) {
    /*
    
      0====++====2
           ||
           ||  t
           ||
      1====++====3

    */
    p_colours[2][0] = p_colours[0][0];
    p_colours[3][1] = p_colours[1][1];
    p_colours[0][1] = p_colours[1][0] = Flow::Counter();
    p_colours[2][1] = p_colours[3][0] = Flow::Counter();
  }
  else {
    if (rr-Mu-Mt < 0.) {
      /*

	0====+\---==3
             ||\ /
             || X u
             ||/ \
	1====+/---==2
	   
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
	   
      */
      p_colours[2][0] = p_colours[0][0];
      p_colours[3][1] = p_colours[0][1] = Flow::Counter();
      p_colours[2][1] = p_colours[1][1];
      p_colours[3][0] = p_colours[1][0] = Flow::Counter();
    }
  }
  return true;
}
