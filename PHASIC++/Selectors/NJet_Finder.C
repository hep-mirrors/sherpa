#ifndef PHASIC_Selectors_NJet_Finder_h
#define PHASIC_Selectors_NJet_Finder_h

#include "ATOOLS/Phys/Particle_List.H"
#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Math/Poincare.H"

namespace PHASIC {
  class NJet_Finder : public Selector_Base {
    double m_pt2min,m_et2min,m_delta_r,m_r2min,m_etamax;
    int    m_exp, m_type;
    double m_ene, m_s, m_sprime; 
    int    m_nstrong;

    double ** p_ktij;
    int    *  p_imap;
    double *  p_kis;
    std::vector<double> m_jetpt;


    double DEta12(ATOOLS::Vec4D &,ATOOLS::Vec4D &);
    double CosDPhi12(ATOOLS::Vec4D &,ATOOLS::Vec4D &);
    double DPhi12(ATOOLS::Vec4D &,ATOOLS::Vec4D &);
    double R2(ATOOLS::Vec4D &,ATOOLS::Vec4D &);
    double Y12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &) const;
    double DCos12(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &) const;

    void   AddToJetlist(const ATOOLS::Vec4D &);

    void   ConstructJets(ATOOLS::Vec4D * p, int n);
  public:
    NJet_Finder(int nin, int nout,ATOOLS::Flavour * fl, 
		double ptmin, double etmin, double dr, int exp, double etamax, int nn);

    ~NJet_Finder();


    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    bool   Trigger(const ATOOLS::Vec4D_Vector &);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &,ATOOLS::NLO_subevtlist *const);

    void   BuildCuts(Cut_Data *);
  };
}

#endif

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"


using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

NJet_Finder::NJet_Finder(int nin, int nout,ATOOLS::Flavour * fl, 
			 double ptmin, double etmin, double dr, int exp, double etamax, int nn) : 
  Selector_Base("NJetfinder"), m_pt2min(sqr(ptmin)), m_et2min(sqr(etmin)), 
  m_delta_r(dr), m_etamax(etamax), m_exp(exp)
{
  m_fl         = fl;
  m_ene        = rpa.gen.Ecms()/2.;
  m_sprime     = m_s = sqr(2.*m_ene); 
  m_smin       = Max(m_pt2min,m_et2min);
  m_smax       = m_s;

  m_r2min      = sqr(m_delta_r);

  m_nin        = nin;
  m_nout       = nout;
  m_n          = nn;

  if (m_nin==2) {
    int instrong(0);
    for (int j=0;j<m_nin;j++) { if (m_fl[j].Resummed()) instrong++; }
    if (instrong==0) m_type = 1;
    if (instrong==1) m_type = 2;
    if (instrong==2) m_type = 4;
  }

  p_kis  = new double[m_nout];
  p_imap = new int[m_nout];
  p_ktij = new double*[m_nout];
  for (int i=0;i<m_nout;++i) p_ktij[i] = new double[m_nout];
  for (int i=0;i<m_nout;++i) p_imap[i] = i;
  
  m_nstrong = 0;
  for (int i=m_nin;i<m_nout+m_nin;i++) {
    if (fl[i].Resummed()) m_nstrong++;
  }

  m_sel_log    = new Selector_Log(m_name);
}


NJet_Finder::~NJet_Finder() {
      for (int i=0;i<m_nout;++i) delete [] p_ktij[i];
      delete [] p_ktij;
      delete [] p_imap;
      delete [] p_kis;
}


void NJet_Finder::AddToJetlist(const Vec4D &jet) {
  if (dabs(jet.Eta())<m_etamax) {
    if (jet.EPerp2()>=m_et2min&&jet.PPerp2()>=m_pt2min) m_jetpt.push_back(jet.PPerp2());
  }
}

bool NJet_Finder::NoJetTrigger(const Vec4D_Vector &p)
{
  double s=(p[0]+p[1]).Abs2();
  return (s>m_smin*4.);
}

bool NJet_Finder::Trigger(const Vec4D_Vector &p)
{
  if (m_n<1) return true;

  // create copy
  m_jetpt.clear();
  int n=0;
  Vec4D * moms = new Vec4D[m_nout];
  for (int i=m_nin;i<m_nout+m_nin;i++) {
    if (m_fl[i].Resummed()) {
      moms[n]=p[i];
      n++;
    }
  }

  // cluster
  for (int i=0;i<n;++i) p_imap[i] = i;
  ConstructJets(moms,n);

  delete [] moms;
  if (n<m_n) return 0;

  bool trigger(true);
  if (m_jetpt.size()<size_t(m_n)) trigger=false;

  return (1-m_sel_log->Hit(1-trigger));
}

bool NJet_Finder::JetTrigger(const Vec4D_Vector &p,NLO_subevtlist *const subs)
{
  if (m_n<1) return true;

  // create copy
  m_jetpt.clear();
  int n=0, ns=subs->back()->m_n-m_nin;
  Vec4D * moms = new Vec4D[ns];
  for (int i=m_nin;i<ns+m_nin;i++) {
    if (subs->back()->p_fl[i].Resummed()) {
      moms[n]=p[i];
      n++;
    }
  }

  // cluster
  for (int i=0;i<n;++i) p_imap[i] = i;

  ConstructJets(moms,n);

  delete [] moms;
  if (n<m_n) return 0;

  bool trigger(true);
  if (m_jetpt.size()<size_t(m_n)) trigger=false;
  
  return (1-m_sel_log->Hit(1-trigger));
}

void NJet_Finder::ConstructJets(Vec4D * p, int n)
{
  if (n==0) return;
  if (n==1) {
    AddToJetlist(p[0]);
    return;
  }

  //cal first matrix
  int ii=0, jj=0;
  double dmin=(m_type>1)?p[0].PPerp2():sqr(p[0][0]);
  dmin=pow(dmin,m_exp);
  {
    
    for (int i=0;i<n;++i) {
      double di = (m_type>1)?p[i].PPerp2():sqr(p[i][0]);
      p_ktij[i][i] = di = pow(di,m_exp);
      if (di<dmin) { dmin=di; ii=i; jj=i;}
      for (int j=0;j<i;++j) {
	double dj  = p_ktij[j][j]; 
	double dij = p_ktij[i][j] = Min(di,dj)*R2(p[i],p[j]) /m_r2min;
	if (dij<dmin) {dmin=dij; ii=i; jj=j;}
      }
    }
  }
  // recalc matrix
  while (n>0) {
    if (ii!=jj) {
      // combine precluster
      p[p_imap[jj]]+=p[p_imap[ii]];
    }
    else {
      // add to jet list
      AddToJetlist(p[p_imap[ii]]);
    }

    --n;

    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];

    if (n==1) {
      int jjx=p_imap[jj];
      p_ktij[jjx][jjx] = pow((m_type>1)?p[jjx].PPerp2():sqr(p[jjx][0]),m_exp);
      break;
    }
    // update map (remove precluster)
    {

    
    // update matrix (only what is necessary)
    int jjx=p_imap[jj];
    p_ktij[jjx][jjx] = pow((m_type>1)?p[jjx].PPerp2():sqr(p[jjx][0]),m_exp);
    for (int j=0;j<jj;++j)   p_ktij[jjx][p_imap[j]] = 
			       Min(p_ktij[jjx][jjx],p_ktij[p_imap[j]][p_imap[j]])
			       *R2(p[jjx],p[p_imap[j]])/m_r2min;
    for (int i=jj+1;i<n;++i) p_ktij[p_imap[i]][jjx] = 
			       Min(p_ktij[jjx][jjx],p_ktij[p_imap[i]][p_imap[i]])
			       *R2(p[p_imap[i]],p[jjx])/m_r2min;
    }
    // redetermine rmin and dmin
    ii=jj=0;
    {

    dmin=p_ktij[p_imap[0]][p_imap[0]];
    for (int i=0;i<n;++i) {
      int ix=p_imap[i];
      double di = p_ktij[ix][ix];
      if (di<dmin) { dmin=di; ii=jj=i;}
      for (int j=0;j<i;++j) {
	int jx=p_imap[j];
	double dij = p_ktij[ix][jx];
	if (dij<dmin) { dmin=dij; ii=i; jj=j;}
      }
    }
    }
  }

  // add remaining preclusters to jetlist
  for (int i=0;i<n;++i) {
    AddToJetlist(p[p_imap[i]]);
  }
}


void NJet_Finder::BuildCuts(Cut_Data * cuts) 
{
  return;
}


/*----------------------------------------------------------------------------------

  Utilities

  ----------------------------------------------------------------------------------*/
double NJet_Finder::DEta12(Vec4D & p1,Vec4D & p2)
{
  // eta1,2 = -log(tan(theta_1,2)/2)   
  double c1=p1[3]/Vec3D(p1).Abs();
  double c2=p2[3]/Vec3D(p2).Abs();
  return  0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
}

double NJet_Finder::CosDPhi12(Vec4D & p1,Vec4D & p2)
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return (p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2);
}

double NJet_Finder::DPhi12(Vec4D & p1,Vec4D & p2)
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
}


double NJet_Finder::R2(Vec4D &p1, Vec4D &p2)
{
  return sqr(p1.Y()-p2.Y()) + sqr(DPhi12(p1,p2));
}

double NJet_Finder::Y12(const Vec4D & p1, const Vec4D & p2) const
{
  return 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2))/m_sprime;
}

double NJet_Finder::DCos12(const Vec4D & p1,const Vec4D & p2) const
{
  double s  = p1[1]*p2[1]+p1[2]*p2[2]+p1[3]*p2[3];
  double b1 = p1[1]*p1[1]+p1[2]*p1[2]+p1[3]*p1[3];
  double b2 = p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3];
  return s/sqrt(b1*b2);
  //  return Vec3D(p1)*Vec3D(p2)/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

namespace PHASIC{

DECLARE_ND_GETTER(NJet_Finder_Getter,"NJetFinder",Selector_Base,Selector_Key,true);

Selector_Base *NJet_Finder_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<4) THROW(critical_error,"Invalid syntax");
  
  int exp(1);
  if (key.front().size()>=5) exp=ToType<int>(key[0][4]);
  double etamax(100.);
  if (key.front().size()>=6) etamax=ToType<double>(key[0][5]);

  NJet_Finder *jf(new NJet_Finder(key.p_proc->NIn(),key.p_proc->NOut(),
				  (Flavour*)&key.p_proc->Process()->Flavours().front(),
				  ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1])),
				  ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2])),
				  ToType<double>(key[0][3]),exp,etamax,
				  ToType<int>(key[0][0])));
  jf->SetProcess(key.p_proc);
  return jf;
}

void NJet_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"NJetFinder n ptmin etmin dr [exp=1] [etamax=100]"; 
}

}
