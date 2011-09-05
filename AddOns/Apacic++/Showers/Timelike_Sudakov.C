#include "AddOns/Apacic++/Showers/Timelike_Sudakov.H"

#include "AddOns/Apacic++/Showers/QCD_Splitting_Functions.H"
#include "AddOns/Apacic++/Showers/QED_Splitting_Functions.H"
#include "AddOns/Apacic++/Showers/SUSY_QCD_Splitting_Functions.H"
#include "AddOns/Apacic++/Showers/Timelike_Kinematics.H"
#include "AddOns/Apacic++/Showers/Sudakov_Tools.H"
#include "AddOns/Apacic++/Main/Knot.H"
#include "ATOOLS/Math/Random.H"

#include <iomanip>

using namespace APACIC;
using namespace ATOOLS;

Timelike_Sudakov::Timelike_Sudakov(Timelike_Kinematics *const kin,
				   MODEL::Model_Base *const model) :
  Splitting_Group(p_rms), p_tools(new Sudakov_Tools(model)), 
  p_kin(kin), p_rms(NULL),
  m_pt2min(1.0), m_pt2max(sqr(rpa->gen.Ecms())), m_t0(4.0),
  m_pt2min_qed(0.0025), m_t0_qed(0.01), m_rbmax(1.0), m_fcs(0.0) {}

Timelike_Sudakov::~Timelike_Sudakov()
{
  if (p_tools) delete p_tools;
}

void Timelike_Sudakov::Init(const double fmed)
{
  // std::cout<<"Init Sudakov : "<<fmed<<std::endl;
  // for static couplings, set first argument to 0.
  p_tools->CalculateMaxCouplings
    (m_cpl_scheme,m_pt2min*m_rscalefac,m_pt2max*m_rscalefac);
  for (int i=1;i<17;++i) {
    if (i==7) i=11;
    Flavour fl = Flavour((kf_code)(i));
    if (fl.IsOn()) {
      if (fl.Strong()) {
	Add(new q_qg(p_rms,fl,p_tools,fmed));
	Add(new q_qg(p_rms,fl.Bar(),p_tools,fmed));
	if (fl.Mass()<100.) Add(new g_qq(p_rms,fl,p_tools,fmed));
      }
      if (!(fl.Charge()==0) && (m_direct_photons)) {
	Add(new f_fp(p_rms,fl,p_tools));
	Add(new f_fp(p_rms,fl.Bar(),p_tools));
	if (fl.Mass()<100.) Add(new p_ff(p_rms,fl,p_tools));
      }
    }
  }
  Add(new g_gg(p_rms,p_tools,fmed));
  
  //susy splitting functions
  if (MODEL::s_model->Name()==std::string("MSSM")) {
    for (short int i=1;i<3;i++) {
      for (short int j=1;j<7;j++) {
	Flavour fl = Flavour((kf_code)(i*1000000 + j));
	if (fl.IsOn()) {
	  Add(new SQuark__SQuark_Gluon(p_rms,fl,p_tools));
	  Add(new SQuark__SQuark_Gluon(p_rms,fl.Bar(),p_tools));
	  Add(new Gluino__Gluino_Gluon(p_rms,p_tools));
        }
      }
    }
  }
  PrintStat();
}

void Timelike_Sudakov::Add(Splitting_Function *spl) 
{
  Splitting_Vector::iterator sit(m_splittings.begin());
  for (;sit!=m_splittings.end();++sit) {
    if ((*sit)->GetFlA()==spl->GetFlA()) {
      (*sit)->Add(spl);
      return;
    }
  }
  m_splittings.push_back(new Splitting_Group(p_rms,spl));
  p_selected=spl;
}

double Timelike_Sudakov::CrudeInt(double zmin, double zmax) 
{
  Splitting_Vector::iterator sit(m_splittings.begin());
  for (;sit!=m_splittings.end();++sit)
    if ((*sit)->GetFlA()==m_inflav) { 
      p_selected=*sit; 
      break;
    }
  if (sit==m_splittings.end()) return m_lastint=-1.0;
  return m_lastint=p_selected->CrudeInt(zmin,zmax)*m_rbmax;
}        

void Timelike_Sudakov::SelectOne() 
{ 
  p_selected->SelectOne(); 
}

void Timelike_Sudakov::AcceptEmission(const Knot *const d) 
{
  ++m_nem;
  p_lemkn=d;
}

void Timelike_Sudakov::AcceptBranch(const Knot *const mo) 
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"):"<<std::endl;
  Knot *d1(mo->left), *d2(mo->right);
  if (!mo->part->Flav().Strong()) {
    d1->maxpt2 = d2->maxpt2 = mo->maxpt2;
    d1->thcrit = d2->thcrit = mo->thcrit;
    msg_Debugging()<<"  accept ew "<<mo->kn_no<<"->("<<d1->kn_no<<","
		   <<d2->kn_no<<"), set pt2 = "<<d1->maxpt2
		   <<", th = "<<d1->thcrit<<" from "<<mo->sthcrit<<"\n";
    return;
  }
  d1->maxpt2 = d2->maxpt2 = mo->smaxpt2;
  d1->thcrit = d2->thcrit = mo->sthcrit;
  msg_Debugging()<<"  accept "<<mo->kn_no<<"->("<<d1->kn_no
		 <<","<<d2->kn_no<<"), set pt2 = "<<d1->maxpt2
		 <<", th = "<<d1->thcrit<<" from "<<mo->sthcrit<<"\n";
}

bool Timelike_Sudakov::Dice(Knot *const mother, Knot *const granny) 
{
  switch (m_evolution_scheme) {
  case 0: return DiceT(mother,granny);
  }
  THROW(fatal_error,"Invalid evolution scheme.");
  return false;
}

bool Timelike_Sudakov::DiceT(Knot *const mother, Knot *const granny) 
{
  msg_Debugging()<<METHOD<<"(): m_nem = "<<m_nem
		 <<" vs. m_maxem = "<<m_maxem<<"\n";
  if (m_nem>=m_maxem) return false;
  m_inflav=mother->part->Flav(); 
  double t0=m_inflav.Strong()?m_t0:m_t0_qed;
  m_oldt=m_t=mother->t;                  
  m_E2=mother->E2;
  m_shower=mother->shower;
  if (m_t-t0<rpa->gen.Accu()) return false;
  double z0;
  while (m_t>t0) {
    z0=Max(0.5*(1.-sqrt(1.-t0/(m_t-mother->tout))),1.e-6);
    CrudeInt(z0,1.-z0);
    if (m_lastint<0.0) return false;
    if (m_mass_scheme&2) m_t-=mother->tout;
    ProduceT(m_t,mother->tout);
    if (m_t<t0) return 0;
    if (m_mass_scheme&2) m_t+=mother->tout;
    // determine estimate for energy
    if (granny) m_E2=0.25*sqr(m_t+granny->E2)/granny->E2;
    if (m_t<m_E2) {
      SelectOne();
      m_z=GetZ();
      if (!Veto(mother,m_t)) {
	if (m_azimuthal_correlation) m_phi=p_selected->GetPhi(m_z);
	else m_phi=UniformPhi();
	mother->zs=mother->z=m_z;
	mother->t=m_t;
	mother->phi=m_phi;
	return true;
      }
    }
  }
  return false; 
}

void Timelike_Sudakov::ProduceT(double ta, double m2) 
{
  m_t*=exp(2.0*M_PI*log(ran.Get())/m_lastint);
}

bool Timelike_Sudakov::MEWeight(Knot *const me)
{
  msg_Debugging()<<METHOD<<"(): "<<*me;
  if (me->left==NULL || me->decay==NULL || 
      me->decay->kn_no!=me->kn_no) return true;
  if (!(me->part->Flav().Strong() &&
	me->left->part->Flav().Strong() &&
	me->right->part->Flav().Strong())) return true;
  double t(m_t=me->part->Momentum().Abs2());
  double z(me->left->part->Momentum()[0]/me->part->Momentum()[0]);
  Flavour flb(me->left->part->Flav()), flc(me->right->part->Flav());
  m_E2=sqr(me->part->Momentum()[0]); 
  if (me->prev) m_E2=0.25*sqr(t+me->prev->E2)/me->prev->E2;
  m_tb=me->left->tout; m_tc=me->right->tout;
  switch (m_kt_scheme) {
  case 0: // durham scheme
    m_pt2=0.5*Min(z*z*m_E2-m_tb,(1.0-z)*(1.0-z)*m_E2-m_tc)*
      (1.0-cos(p_kin->GetOpeningAngle(z,m_E2,t,m_tb,m_tc)));
    break;
  case 1: // relative transverse momentum in lc kinematics
    m_pt2=p_kin->GetRelativeKT2(z,m_E2,t,m_tb,m_tc);
    break;
  case 2: { // case 1 w/ zero daughter masses
    double zlc(p_kin->LightConeZ(z,m_E2,t,m_tb,m_tc));
    m_pt2=zlc*(1.0-zlc)*t;
    break;
  }
  case 3: // pseudo transverse momentum
    m_pt2=z*(1.0-z)*t-(1.0-z)*m_tb-z*m_tc;
    break;
  case 4: {// durham scheme with light cone kinematics
    double zlc(p_kin->LightConeZ(z,m_E2,t,m_tb,m_tc));
    double f1(zlc/(1.0-zlc)), f2((1.0-zlc)/zlc);
    double kt2(p_kin->GetRelativeKT2(z,m_E2,t,m_tb,m_tc));
    double c1(flb.IsQuark()?1.0:f1);
    double c2(flc.IsQuark()?1.0:f2);
    m_pt2=0.5*Min(c1,c2)*(f1*m_tc+f2*m_tb+kt2/(zlc*(1.0-zlc)));
    break;
  }
  default:
    THROW(fatal_error,"No kt definition.");
  }
  msg_Debugging()<<" -> pt = "<<sqrt(m_pt2)<<", as = "
		 <<GetCoupling(m_rscalefac*m_pt2)
		 <<" vs. "<<me->asme<<"\n";
  double asc(0.0);
  switch (m_cpl_scheme) {
  case 2 : asc=GetCoupling(m_rscalefac*0.25*t); break;
  default: asc=GetCoupling(m_rscalefac*m_pt2); break;
  }
  if (/*4.0*m_pt2<sqr(me->qjv) ||*/ asc>me->asme) return false;
  return asc/me->asme>ran.Get();
}

bool Timelike_Sudakov::Veto(Knot *const mo,double t) 
{
  // branch okay and mother timelike
  m_last_veto=0;
  // kincheck for first timelike from spacelike
//   if (mo->tmax!=0 && m_t>mo->tmax) return true;
  if (m_E2<t) return true;
 
 //  sum m1 + m2 < sqrt(ta)
  m_tb=sqr(p_ms->Mass(GetFlB()));
  m_tc=sqr(p_ms->Mass(GetFlC()));

  if (mo->shower==3) {
    m_tb=t;
    t=mo->tmo;
  }

  m_last_veto=1;
  if (t<sqr(sqrt(m_tb)+sqrt(m_tc))) return true;

  double z(p_kin->GetZ(m_z,t,m_tb,m_tc));
  if (z<0.0 || z>1.0) return true;

  switch (m_kt_scheme) {
  case 0: // durham scheme
    m_pt2=0.5*Min(z*z*m_E2-m_tb,(1.0-z)*(1.0-z)*m_E2-m_tc)*
      (1.0-cos(p_kin->GetOpeningAngle(z,m_E2,t,m_tb,m_tc)));
    break;
  case 1: // relative transverse momentum in lc kinematics
    m_pt2=p_kin->GetRelativeKT2(z,m_E2,t,m_tb,m_tc);
    break;
  case 2: { // case 1 w/ zero daughter masses
    double zlc(p_kin->LightConeZ(z,m_E2,t,m_tb,m_tc));
    m_pt2=zlc*(1.0-zlc)*t;
    break;
  }
  case 3: // pseudo transverse momentum
    m_pt2=z*(1.0-z)*t-(1.0-z)*m_tb-z*m_tc;
    break;
  case 4: {// durham scheme with light cone kinematics
    double zlc(p_kin->LightConeZ(z,m_E2,t,m_tb,m_tc));
    double f1(zlc/(1.0-zlc)), f2((1.0-zlc)/zlc);
    double kt2(p_kin->GetRelativeKT2(z,m_E2,t,m_tb,m_tc));
    double c1(GetFlB().IsQuark()?1.0:f1);
    double c2(GetFlC().IsQuark()?1.0:f2);
    m_pt2=0.5*Min(c1,c2)*(f1*m_tc+f2*m_tb+kt2/(zlc*(1.0-zlc)));
    break;
  }
  default:
    THROW(fatal_error,"No kt definition.");
  }
  if (m_inflav.Strong()) {
    if (m_pt2<m_pt2min || m_pt2<mo->minpt2) return true;
  }
  else {
    if (m_pt2<m_pt2min_qed) return true;
  }
  // timelike daughters
  m_last_veto=2;
  double wb(z*z*m_E2), wc((1.-z)*(1.-z)*m_E2);
  if (m_tb>wb || m_tc>wc) return true;

  // z-range and splitting function
  m_last_veto=3;
  if (MassVeto(t,m_E2,z)) return true;

  // 2. alphaS
  m_last_veto=4;
  if (CplVeto(t)) return true;

  // 3. angular ordering
  m_last_veto=5;
  if (OrderingVeto(mo,t,m_E2,z)) return true; 

  // 5. ME
  m_last_veto=7;
  if (MEVeto(mo,t)) return true;
  m_last_veto=-1;
  return false;
}

double sql(const double &s,const double &s1,const double &s2)
{
  return (sqr(s-s1-s2)-4.0*s1*s2)/(4.0*s);
}

bool Timelike_Sudakov::MassVeto(double t, double E2,double z) 
{
  double psw(1.0);
  if (m_shower==3) {
    // additional phasespace weight
    psw=m_oldt/m_t*sqrt(sql(m_t,m_tb,m_tc)/sql(m_oldt,m_tb,m_tc));
  }
  if (GetWeight(m_z,m_pt2,m_mass_scheme&1)*psw<ran.Get()) 
    return true;
  if ((m_width_scheme>0) && (sqr(m_inflav.Width())>0.)) {
    if (m_width_scheme==1 && m_pt2<sqr(m_inflav.Width())) return true;
    else if (m_pt2/(m_pt2+sqr(m_inflav.Width()))<ran.Get()) return true;
  }
  return false;
}

bool Timelike_Sudakov::CplVeto(double t) 
{
  if (m_fcs>0.0) return 
    GetCoupling(m_rscalefac*m_fcs)/GetCoupling()<ran.Get();   
  switch (m_cpl_scheme) {
  case 0 : return false;
  case 2 : 
    return GetCoupling(m_rscalefac*0.25*t)/GetCoupling()<ran.Get();   
  default : 
    return GetCoupling(m_rscalefac*m_pt2)/GetCoupling()<ran.Get();   
  }
  return true;
}

bool Timelike_Sudakov::OrderingVeto(Knot * mo,double t, double E2, double z) 
{
  double th(p_kin->GetOpeningAngle(z,E2,t,m_tb,m_tc));
  mo->sthcrit=th;
  mo->smaxpt2=m_pt2;
  msg_Debugging()<<"ts("<<mo->kn_no
		 <<"): thcrit = "<<mo->thcrit<<", th = "<<th<<std::endl;
  if (!m_inflav.Strong()) return false;
  switch (m_ordering_scheme) {
  case 0 : return false;
  case 2 : 
    if (m_pt2<mo->maxpt2) return false;
    return true;
  case 1 :
  default :
    if (mo->shower==3) {
      if (th>mo->thcrit || 3.14>mo->thcrit) return false; 
      return true;
    }
    if (th<mo->thcrit || 3.14<mo->thcrit) return false; 
    return true;
  }
  return true;
}

bool Timelike_Sudakov::MEVeto(Knot * mo,double t) 
{
  //the QED emission me veto
  if (m_qed_mecorr_scheme==1) {
    bool qed_veto=true;
    Knot *gr(mo->prev);
    if (gr->t<0) qed_veto = false;
    if (gr->prev && qed_veto) qed_veto = false;
    if (!m_inflav.IsFermion() && qed_veto) qed_veto = false;
    // determine which is the current twig of the tree:
    if (qed_veto) {
      Knot *twig(mo);
      while (gr->prev) {
	twig=gr;
	gr=gr->prev;
      }
     /*
       if "first" branch perform ME - Correction for photon radiation
      (t',z') *      
             / \      x_i = 2 E_i / sqrt(t')
      (t,z) *   \
           / \   \
          1   3   2
   */
      double mass123(gr->t), mass13(t);
      if (m_ordering_scheme==0) {
	mass123/=gr->z*(1.-gr->z);
	mass13/=m_z*(1.-m_z);
      }
      double x2(1.-mass13/mass123), x1(m_z*(2.-x2)), x3(2.-x1-x2); 
      double ds_ps((1.-x1)/x3*(1.+sqr(x1/(2.-x2))) + 
		   (1.-x2)/x3*(1.+sqr(x2/(2.-x1))));
      double ds_me(sqr(x1)+sqr(x2)), ratio(ds_me/ds_ps);
      qed_veto =  ratio<ran.Get()?true:false;
      if (qed_veto) return true;
    }
  }
    
  //the QCD emission me veto
  if (!m_inflav.Strong()) return false;
  Knot *gr(mo->prev);
  if (gr->t<0) return false;
  if (m_mecorr_scheme == 0 || gr==NULL) return false;
  if (m_mecorr_scheme == 2 && gr->prev) return false;
  if (m_mecorr_scheme == 1 && m_pt2<mo->maxpt2) return false;
  if (!m_inflav.IsQuark()) return false;
  // determine which is the current twig of the tree:
  Knot *twig(mo);
  while (gr->prev) {
    twig=gr;
    gr=gr->prev;
  }
  /*
     if "first" branch perform ME - Correction for gluon radiation
     (t',z') *      
            / \      x_i = 2 E_i / sqrt(t')
     (t,z) *   \
          / \   \
         1   3   2
  */
  double mass123(gr->t), mass13(t);
  if (m_ordering_scheme==0) {
    mass123/=gr->z*(1.-gr->z);
    mass13/=m_z*(1.-m_z);
  }
  double x2(1.-mass13/mass123), x1(m_z*(2.-x2)), x3(2.-x1-x2); // quark or antiquark !!! 
  double ds_ps((1.-x1)/x3*(1.+sqr(x1/(2.-x2))) + 
              (1.-x2)/x3*(1.+sqr(x2/(2.-x1))));
  double ds_me(sqr(x1)+sqr(x2)), ratio(ds_me/ds_ps);
  return ratio<ran.Get()?true:false;
}

const Simple_Polarisation_Info &Timelike_Sudakov::GetPolB()
{
  if (m_azimuthal_correlation) m_pol_b=p_selected->GetPolB(m_z,m_phi);
  return m_pol_b;
}

const Simple_Polarisation_Info &Timelike_Sudakov::GetPolC()
{
  if (m_azimuthal_correlation) 
    m_pol_c=p_selected->GetPolC(m_z,m_phi,m_pol_b.Angle());
  return m_pol_c;
}

double Timelike_Sudakov::UniformPhi() const 
{ 
  return 2.0*M_PI*ATOOLS::ran.Get(); 
}

void Timelike_Sudakov::SetPT2Min(const double &pt2)      
{
  m_pt2min=pt2; 
  m_t0=4.0*m_pt2min;
}

#ifdef CHECK_SPLITTINGS
struct Spl_Data {
  int masses;
  double ta;
  std::vector<double> values;
};

void Timelike_Sudakov::CheckSplittings() 
{
  const int dsize=7;
  std::vector<Spl_Data> data(dsize);
  data[0].masses=0;   data[0].ta=0;           // z values
  data[1].masses=0;   data[1].ta=sqr(8.);    // weights
  data[2].masses=1;   data[2].ta=sqr(8.);
  data[3].masses=0;   data[3].ta=sqr(12.);    // weights
  data[4].masses=1;   data[4].ta=sqr(12.);
  data[5].masses=0;   data[5].ta=sqr(50.);    // weights
  data[6].masses=1;   data[6].ta=sqr(50.);
  for (double z=0;z<=1.;z+=10.e-2) {
    data[0].values.push_back(z);
    for (int j=1;j<dsize;++j) {
      data[j].values.push_back(0.);
    }
  }  

  for (size_t i(0);i<m_splittings.size();++i) {
    Splitting_Function *sf(m_splittings[i]);
    sf->PrintStat();
    // calculation
    for (size_t i=1; i<data.size();++i) {
      for (size_t j=0;j<data[0].values.size();++j) {
	double z=data[0].values[j];
	int masses=data[i].masses;
	double ta=data[i].ta;
	m_tb=sqr(p_ms->Mass(sf->GetFlB()));
	m_tc=sqr(p_ms->Mass(sf->GetFlC()));
	data[i].values[j]=sf->GetWeight(m_z,m_pt2,masses);
      }
    }
    // output 
  }
  THROW(normal_exit,"Finished check.");
}
#endif

