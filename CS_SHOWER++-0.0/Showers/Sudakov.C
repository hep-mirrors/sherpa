#include "Sudakov.H"
#include "QCD_Splitting_Functions.H"
#include "Run_Parameter.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

Sudakov::Sudakov(PDF::ISR_Handler *isr) : 
  m_k0sq(1.), m_asmax(0.5), p_as(NULL)
{
  int hadron = rpa.gen.Beam1().Strong()==1?0:1;
  for (int i=1;i<7;++i) {
    Flavour fl = Flavour(kf::code(i));
    if (fl.IsOn() && fl.Strong()) {
      //the FF case
      Add(new q_qg_FF(fl));
      Add(new q_gq_FF(fl));
      Add(new q_qg_FF(fl.Bar()));
      Add(new q_gq_FF(fl.Bar()));
      if (fl.PSMass()<100.) Add(new g_qq_FF(fl));
      if (isr->On()) {
	//the FI case
	Add(new q_qg_FI(fl,isr->PDF(hadron)));
	Add(new q_qg_FI(fl.Bar(),isr->PDF(hadron)));
	Add(new q_gq_FI(fl,isr->PDF(hadron)));
	Add(new q_gq_FI(fl.Bar(),isr->PDF(hadron)));
	Add(new g_qq_FI(fl,isr->PDF(hadron)));
	//the IF case
	Add(new q_qg_IF(fl,isr->PDF(hadron)));
	Add(new q_qg_IF(fl.Bar(),isr->PDF(hadron)));
	Add(new q_gq_IF(fl,isr->PDF(hadron)));
	Add(new q_gq_IF(fl.Bar(),isr->PDF(hadron)));
	Add(new g_qq_IF(fl,isr->PDF(hadron)));
      }
    }
  }
  Add(new g_gg_FF());
  if (isr->On()) {
    Add(new g_gg_FI(isr->PDF(hadron)));
    Add(new g_gg_IF(isr->PDF(hadron)));
  }

  p_as    = (MODEL::Running_AlphaS*)(rpa.gen.GetScalarFunction("alpha_S"));
  m_asmax = (*p_as)(m_k0sq);
}

Sudakov::~Sudakov() {}

bool Sudakov::Dice(Parton * split,Parton * spect,const double kt_ext) {
  
  m_flavs[0] = split->GetFlavour();
  if (m_flavs[0].IsQuark()) m_nspect = 2.;
                       else m_nspect = 1.;  
  m_type     = cstp::none;
  
  double s=0,t=0;
  cout<<"####################### In  Sudakov::Dice ###########################"<<endl;
  if (split->GetType()==pst::FS && spect->GetType()==pst::FS) {

    s = (split->Momentum()+spect->Momentum()).Abs2();    

    //std::cout<<"FF ---------------------------------------------------------------"<<std::endl
    //	     <<"Try this : k0sq = "<<4.*m_k0sq<<", s = "<<s<<endl;
    if (4.*m_k0sq>s) return false;
    
    m_type=cstp::FF;
    m_deltaz = sqrt(1.-4.*m_k0sq/s);
    m_zmin   = 0.5*(1.-m_deltaz);
    m_zmax   = 0.5*(1.+m_deltaz);
    m_scale  = s/4.;
  }
  
  if (split->GetType()==pst::FS && spect->GetType()==pst::IS) { 
    //for test purposes
    //if (m_flavs[0]==Flavour(kf::gluon)) return false;
    
    t = 2.*(1-spect->Xbj())/spect->Xbj()*split->Momentum()*spect->Momentum();  
    //    std::cout<<"FI ---------------------------------------------------------------"<<std::endl;
    //	     <<"Try this : k0sq = "<<4.*m_k0sq<<", t = "<<t<<" from "<<spect->Xbj()<<endl;

    if (4.*m_k0sq>t) return false;
    
    m_type=cstp::FI;
    m_deltaz = sqrt(1.-4.*m_k0sq/t);
    m_zmin   = 0.5*(1.-m_deltaz);
    m_zmax   = 0.5*(1.+m_deltaz);
    m_scale  = t/4.;
  }
  
  if (split->GetType()==pst::IS && spect->GetType()==pst::FS) { 
    //for test purposes
    //if (m_flavs[0]==Flavour(kf::gluon)) return false;
    
    t = 2.*(1-spect->Xbj())/spect->Xbj()*split->Momentum()*spect->Momentum();  
    //    std::cout<<"FI ---------------------------------------------------------------"<<std::endl;
    //	     <<"Try this : k0sq = "<<4.*m_k0sq<<", t = "<<t<<" from "<<spect->Xbj()<<endl;

    if (4.*m_k0sq>t) return false;
    
    m_type=cstp::IF;
    m_deltaz = sqrt(1.-4.*m_k0sq/t);
    m_zmin   = 0.5*(1.-m_deltaz);
    m_zmax   = 0.5*(1.+m_deltaz);
    m_scale  = t/4.;
  }
  
  if (m_type==cstp::none) {
    std::cout<<split->GetFlavour()<<" "<<spect->GetFlavour()<<std::endl;
    msg.Error()<<"Error in Sudakov::Dice : No type for splitter. "<<endl<<(*split);
    abort();
  }
  if (OverIntegrated(m_zmin,m_zmax,m_scale,spect->Xbj())<=0.) {
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
    //FF
    if (split->GetType()==pst::FS && spect->GetType()==pst::FS) {
      //std::cout<<"FF ---------------------------------------------------------------"<<std::endl;
      m_y = m_kperp2/(s*m_z*(1.-m_z));
      /*
	cout<<m_flavs[0]<<" -> "<<GetFlavourB()<<" + "<<GetFlavourC()
	<<" : { kt^2 = "<<m_kperp2<<", z = "<<m_z<<"} "
	<<m_z*(1.-m_z)<<" > "<<m_kperp2/s<<endl;
	cout<<"Test : "<<m_kperp2<<" "<<ktveto2<<endl;
      */
      if (m_kperp2<Max(m_k0sq,kt_ext))  return false;
      if (Veto(s,0.,ktveto2)) break;
     }
    //FI
    if (split->GetType()==pst::FS && spect->GetType()==pst::IS) {
      //std::cout<<"FI ---------------------------------------------------------------"<<std::endl;
      double ta = 2.*split->Momentum()*spect->Momentum(); 
      m_y = 1./(1.+ta*m_z*(1.-m_z)/m_kperp2);
      /*
	cout<<m_flavs[0]<<" -> "<<GetFlavourB()<<" + "<<GetFlavourC()
        <<" : { kt^2 = "<<m_kperp2<<", z = "<<m_z<<"} "
        <<m_z*(1.-m_z)<<" > "<<m_kperp2/t<<endl;
	cout<<"Test : "<<m_kperp2<<" "<<ktveto2<<endl;
      */
      if (m_kperp2<Max(m_k0sq,kt_ext))  return false;
      if (Veto(t,spect->Xbj(),ktveto2)) break;
    }
  }
  split->SetTest(m_kperp2,m_z,m_y);
  cout<<"Succeed :"<<split->GetFlavour()<<" : "<<m_kperp2<<","<<m_z<<","<<m_y<<" vs "
      <<split->KtTest()<<","<<split->ZTest()<<","<<split->YTest()<<endl;
  return true;
}


double Sudakov::OverIntegrated(const double zmin,const double zmax,const double scale,const double xbj) {
  m_lastint = 0.;
  for (m_splitter=m_splittings.begin();m_splitter!=m_splittings.end();m_splitter++) {
    if ((*m_splitter)->GetFlavourA()==m_flavs[0] &&
	(*m_splitter)->GetType()==m_type) {
      m_lastint += (*m_splitter)->OverIntegrated(zmin,zmax,scale,xbj);
      if (m_lastint>0. && m_lastint <0.) cout<<(*this);    
    }
    else (*m_splitter)->SetLast(0.);
  }
  return m_lastint;  
}

void Sudakov::ProduceT() {
  m_kperp2 *= exp(log(ran.Get())*2.*m_nspect*M_PI/m_asmax/m_lastint);
}

bool Sudakov::Veto(double s,double x,double ktveto2) {
  
  if (!KinCheck(s,x,ktveto2)) return false;
  if (!Splitting(x))         return false;
  if (!Coupling())          return false;
  return true;
}

bool Sudakov::KinCheck(double s,double x, double ktveto2) {
 
  if (m_type==cstp::FF) {
    if (m_kperp2>s/4. || m_y<0. || m_y>1.)    return false;
    if (m_z*(1.-m_z) < m_kperp2/s)            return false;
  }
  if (m_type==cstp::FI) {
    if (m_kperp2>s || m_y<0. || m_y>1.)       return false;
    if (x>1.)                               return false;
    if (m_z*(1.-m_z) < m_kperp2/(s*(1.-x))) return false;
  }
  if (m_kperp2>ktveto2)                       return false;
  
  return true;
}

bool Sudakov::Splitting(double x) {
  
  double wt = 0.;
  
  if (m_type==cstp::FF) wt = RejectionWeight(m_z,m_y);
  if (m_type==cstp::FI) wt = RejectionWeight(m_z,m_y,x,m_scale);
  if (wt>1.) std::cout<<" ERROR : In Sudakov::Splitting, weight is larger than 1: "<<wt<<std::endl;
  //cout<<"   Spl-Weight("<<sqrt(m_kperp2)<<","<<m_z<<","<<m_y<<") = "<<wt<<endl;
  if (ran.Get()>wt) return false; 
  return true;
}

bool Sudakov::Coupling() {
  double wt = (*p_as)(m_kperp2)/m_asmax;
  /*
  if (wt>1.) {
    std::cout<<" ERROR : In Sudakov::Coupling, weight is larger than 1: "<<wt<<std::endl;
    cout<<"   Cpl-Weight("<<sqrt(m_kperp2)<<") :"<<(*p_as)(m_kperp2)<<"/"<<m_asmax<<" = "<<wt<<endl;
  }
  */
  if (ran.Get()>wt) return false;
  return true;
}
