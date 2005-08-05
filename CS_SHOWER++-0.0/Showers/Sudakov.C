#include "Sudakov.H"
#include "QCD_Splitting_Functions.H"
#include "Run_Parameter.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

Sudakov::Sudakov() : 
  m_k0sq(1.), m_asmax(0.5), p_as(NULL)
{
  for (int i=1;i<7;++i) {
    Flavour fl = Flavour(kf::code(i));
    if (fl.IsOn() && fl.Strong()) {
      Add(new q_qg_FF(fl));
      Add(new q_gq_FF(fl));
      Add(new q_qg_FF(fl.Bar()));
      Add(new q_gq_FF(fl.Bar()));
      if (fl.PSMass()<100.) Add(new g_qq_FF(fl));
    }
  }
  Add(new g_gg_FF());
  p_as    = (MODEL::Running_AlphaS*)(rpa.gen.GetScalarFunction("alpha_S"));
  m_asmax = (*p_as)(m_k0sq);
}

Sudakov::~Sudakov() {}

bool Sudakov::Dice(Parton * split,Parton * spect,const double kt_ext) {
  m_flavs[0] = split->GetFlavour();
  if (m_flavs[0].IsQuark()) m_nspect = 2.;
                       else m_nspect = 1.;  
  m_type     = cstp::none;
  double s   = (split->Momentum()+spect->Momentum()).Abs2();
  if (4.*m_k0sq>s) return false;
  //cout<<"##################################################"<<endl
  //   <<"Try this : "<<endl<<(*split)<<(*spect);
  if (split->GetType()==pst::FS && spect->GetType()==pst::FS) {
    m_type=cstp::FF;
    m_deltaz = sqrt(1.-4.*m_k0sq/s);
    m_zmin   = 0.5*(1.-m_deltaz);
    m_zmax   = 0.5*(1.+m_deltaz);
    m_scale  = s/4.;
  }
  if (m_type==cstp::none) {
    msg.Error()<<"Error in Sudakov::Dice : No type for splitter. "<<endl<<(*split);
    abort();
  }
  if (OverIntegrated(m_zmin,m_zmax,m_scale)<=0.) {
    msg.Error()<<"Error in Sudakov::Dice : "<<endl
	       <<"   Integral for SF's<0 : {"<<m_zmin<<","<<m_zmax<<","<<m_scale<<"}"<<endl<<(*split);
    abort();
  }
  m_kperp2       = split->KtStart();
  double ktveto2 = split->KtVeto(); 
  while (m_kperp2>=Max(m_k0sq,kt_ext)) {
    ProduceT();
    SelectOne();
    m_z = Z();
    m_y = m_kperp2/(s*m_z*(1.-m_z));
    //cout<<m_flavs[0]<<" -> "<<GetFlavourB()<<" + "<<GetFlavourC()
    //	<<" : { kt^2 = "<<m_kperp2<<", z = "<<m_z<<"} "
    //	<<m_z*(1.-m_z)<<" > "<<m_kperp2/s<<endl;
    //cout<<"Test : "<<m_kperp2<<" "<<ktveto2<<endl;
    if (Veto(s,ktveto2)) break;
  }
  if (m_kperp2<Max(m_k0sq,kt_ext)) return false;
  split->SetTest(m_kperp2,m_z,m_y);
  //cout<<split<<" : "<<m_kperp2<<","<<m_z<<","<<m_y<<" vs "
  //   <<split->KtTest()<<","<<split->ZTest()<<","<<split->YTest()<<endl;
  return true;
}


double Sudakov::OverIntegrated(const double zmin,const double zmax,const double scale) {
  m_lastint = 0.;
  for (m_splitter=m_splittings.begin();m_splitter!=m_splittings.end();m_splitter++) {
    if ((*m_splitter)->GetFlavourA()==m_flavs[0] &&
	(*m_splitter)->GetType()==m_type) {
      m_lastint += (*m_splitter)->OverIntegrated(zmin,zmax,scale);
      //if (m_lastint>0. && m_lastint <0.) cout<<(*this);    
    }
    else (*m_splitter)->SetLast(0.);
  }
  //cout<<"   Try ("<<m_flavs[0]<<"): "<<zmin<<" "<<zmax<<" "<<scale<<" "<<m_lastint<<endl;
  return m_lastint;  
}

void Sudakov::ProduceT() {
  m_kperp2 *= exp(log(ran.Get())*2.*m_nspect*M_PI/m_asmax/m_lastint);
}

bool Sudakov::Veto(double s,double ktveto2) {
  if (!KinCheck(s,ktveto2)) return false;
  if (!Splitting())         return false;
  if (!Coupling())          return false;
  return true;
}

bool Sudakov::KinCheck(double s,double ktveto2) {
  if (m_kperp2>s/4. || m_y<0. || m_y>1.) return false;
  if (m_kperp2>ktveto2)                  return false;
  if (m_z*(1.-m_z) < m_kperp2/s)         return false;
  return true;
}

bool Sudakov::Splitting() {
  double wt = RejectionWeight(m_z,m_y);
  //cout<<"   Spl-Weight("<<sqrt(m_kperp2)<<","<<m_z<<","<<m_y<<") = "<<wt<<endl;
  if (ran.Get()>wt) return false; 
  return true;
}

bool Sudakov::Coupling() {
  double as = (*p_as)(m_kperp2);
  //cout<<"   Cpl-Weight("<<sqrt(m_kperp2)<<") :"<<as<<"/"<<m_asmax<<" = "<<as/m_asmax<<endl;
  if (ran.Get()>as/m_asmax) return false;
  return true;
}
