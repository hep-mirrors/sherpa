#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "AHADIC++/Tools/Splitting_Tools.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AHADIC;
using namespace ATOOLS;

size_t Dipole::s_cnt=0;
const Vec3D Splitting_Tools::s_ex(Vec3D(1.,0.,0.));
const Vec3D Splitting_Tools::s_ey(Vec3D(0.,1.,0.));
const Vec3D Splitting_Tools::s_ez(Vec3D(0.,0.,1.));

PTOrder::code AHADIC::DefinePTOrder(const int & ptorder) {
  switch (ptorder) {
  case 3:
    return PTOrder::total;
  case 2:
    return PTOrder::gluon_split;
  case 1:
    return PTOrder::gluon_emit;
  case 0:
  default:
    break;
  }
  return PTOrder::none;
}


Splitting_Tools::
Splitting_Tools(const leading::code & lead,const PTOrder::code & ptorder,
		const ZForm::code & zform,MODEL::Strong_Coupling * as,
		const bool & analyse) :
  m_leading(lead), m_ptorder(ptorder),
  m_fourquarks(false), m_analyse(true),
  m_masstreatment(int(hadpars->Get(std::string("Mass_treatment")))), 
  p_as(as), p_kernels(new Splitting_Functions(zform,m_masstreatment)), 
  p_options(NULL), p_spect(NULL), p_split(NULL), p_out1(NULL), p_out2(NULL),
  m_mmin_2(sqr(hadpars->GetConstituents()->MinMass())),
  m_pt2max(sqr(hadpars->Get(std::string("ptmax")))), 
  m_pt2max_factor(sqr(hadpars->Get(std::string("ptmax_factor")))), 
  m_lastpt2(-1.), 
  m_tot(0),m_d(0),m_s(0),m_u(0),m_reject_y(0),m_reject_z(0)
{ 
  if (m_analyse) {
    m_histograms[std::string("Splitting_Trials")]        = new Histogram(0,0.,1000,500);
    m_histograms[std::string("Enforced_Trials")]         = new Histogram(0,0.,1000,500);
    m_histograms[std::string("Kinematics_Trials")]       = new Histogram(0,0.,1000,500);
    m_histograms[std::string("PT_Gluon_Splitting")]      = new Histogram(0,0.,10.,1000);
    m_histograms[std::string("PT_Gluon_Emission")]       = new Histogram(0,0.,10.,1000);
    m_histograms[std::string("PT_Gluon_First")]          = new Histogram(0,0.,250.,1000);
    m_histograms[std::string("PT_Gluon_Leading")]        = new Histogram(0,0.,250.,1000);
    m_histograms[std::string("PT2_Gluon_Splitting")]     = new Histogram(0,0.,100.,1000);
    m_histograms[std::string("PT2_Gluon_Emission")]      = new Histogram(0,0.,100.,1000);
    m_histograms[std::string("Z_Gluon_Splitting")]       = new Histogram(0,0.,1.,50);
    m_histograms[std::string("Z_Gluon_Emission")]        = new Histogram(0,0.,1.,50);
    m_histograms[std::string("Z_Gluon_First")]           = new Histogram(0,0.,1.,50);
    m_histograms[std::string("Z_Gluon_Leading")]         = new Histogram(0,0.,1.,50);
    m_histograms[std::string("XB_bquark_before_gsplit")] = new Histogram(0,0.,1.,50);
    m_histograms[std::string("XB_bquark_before_qsplit")] = new Histogram(0,0.,1.,50);
    m_histograms[std::string("XB_bquark_after_gsplit")]  = new Histogram(0,0.,1.,50);
    m_histograms[std::string("XB_bquark_after_qsplit")]  = new Histogram(0,0.,1.,50);
  }
}


Splitting_Tools::~Splitting_Tools() { 
  if (m_analyse) {
    Histogram * histo;
    std::string name;
    for (std::map<std::string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = 
	std::string("Fragmentation_Analysis/")+hit->first+
	std::string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
  delete p_kernels;
}

void Splitting_Tools::SetSpectatorAndSplitter(Dipole * dip) {
  p_split=p_spect=NULL;
  if (dip->Triplet()->m_flav.IsGluon() && 
      !dip->AntiTriplet()->m_flav.IsGluon()) {
    // asymmetric case: triplet splits (gluon)
    p_spect = dip->AntiTriplet(); 
    p_split = dip->Triplet();
    dip->SetSwitched(true);
  }
  else if (dip->AntiTriplet()->m_flav.IsGluon() && 
	   !dip->Triplet()->m_flav.IsGluon()) {
    // asymmetric case: antitriplet splits (gluon)
    p_split = dip->AntiTriplet(); 
    p_spect = dip->Triplet();
    dip->SetSwitched(false);
  }
  else if ((dip->Triplet()->m_info=='B' && dip->AntiTriplet()->m_info!='B') ||
	   (dip->Triplet()->m_info=='L' && dip->AntiTriplet()->m_info=='l')) {
    // asymmetric case: triplet splits 
    // (leading/beam remnant with non leading spectator)
    p_split = dip->Triplet(); 
    p_spect = dip->AntiTriplet();
    dip->SetSwitched(true);
  }
  else if ((dip->AntiTriplet()->m_info=='B' && dip->Triplet()->m_info!='B') ||
	   (dip->AntiTriplet()->m_info=='L' && dip->Triplet()->m_info=='l')) {
    // asymmetric case: triplet splits 
    // (leading/beam remnant with non leading spectator)
    p_split = dip->Triplet(); 
    p_spect = dip->AntiTriplet();
    dip->SetSwitched(true);
  }
  else {
    // symmetric case: pick splitter according to mass.
    double l1(sqr(Max(dip->Triplet()->m_mass,1.e-4)));
    double l2(sqr(Max(dip->AntiTriplet()->m_mass,1.e-4)));
    double disc = l1/(l1+l2);
    if (ran->Get()<disc) {
      p_spect = dip->AntiTriplet(); 
      p_split = dip->Triplet();
      dip->SetSwitched(true);
    }
    else {
      p_split = dip->AntiTriplet(); 
      p_spect = dip->Triplet();
      dip->SetSwitched(false);
    }
  }
}


void Splitting_Tools::SwapSpectatorAndSplitter(Dipole * dip) {
  Proto_Particle * help(p_split); p_split = p_spect; p_spect = help;
  if (dip->IsSwitched()) dip->SetSwitched(false);
  else dip->SetSwitched(true); 
}

bool Splitting_Tools::
PrepareKinematics(Dipole * dip,const bool & first,const bool & enforce) {
  if (dip->MassBar2()<4.*m_mmin_2) return false;

  m_mom1   = p_spect->m_mom;
  m_mom3   = p_split->m_mom;
  m_mom0   = m_mom1+m_mom3;
  m_Q2     = m_mom0.Abs2();
  if (m_Q2<0 || IsNan(m_Q2)) {
    msg_Error()<<"Error in "<<METHOD<<" cannot prepare kinematics for "
	       <<"   "<<m_mom1<<" + "<<m_mom3<<" from \n";
    dip->Output();
  }
  m_Q      = sqrt(m_Q2);
  m_cms    = Poincare(m_mom0);
  m_nperp  = Vec4D(0.0,cross(Vec3D(m_mom3),Vec3D(m_mom1)));
  if (m_nperp.PSpat2()<=1.0e-6) {
    msg_Debugging()<<"Set fixed n_perp\n";
    m_nperp=Vec4D(0.0,1.0,1.0,0.0);
    Poincare zrot(m_mom3,Vec4D::ZVEC);
    zrot.RotateBack(m_nperp);
  }
  m_nperp *= 1.0/m_nperp.PSpat();
  Vec4D help3(m_mom3);
  m_cms.Boost(help3);
  m_lperp  = Vec4D(0.0,cross(Vec3D(help3),Vec3D(m_nperp)));
  m_lperp *= 1.0/m_lperp.PSpat();

  m_m1    = p_spect->m_flav.HadMass();
  m_m1_2  = sqr(m_m1);
  m_m23_2 = sqr(p_split->m_flav.HadMass());

  m_leadsplit = p_split->m_info=='L' || p_split->m_info=='B';
  m_glusplit  = p_split->m_flav.IsGluon();

  switch (m_ptorder) {
  case PTOrder::total:
    m_kt2max = p_split->m_kt2max;
    break;
  case PTOrder::gluon_split:
    m_kt2max = m_glusplit?p_split->m_kt2max:m_Q2;
    break;
  case PTOrder::gluon_emit:
    m_kt2max = m_glusplit?m_Q2:p_split->m_kt2max;
    break;
  case PTOrder::none:
  case PTOrder::flat:
  default:
    m_kt2max = m_Q2;
  }
  return true;
}

bool Splitting_Tools::SelectFlavour(const bool & vetodiquark)
{
  if (m_glusplit) {
    m_flav = Flavour(kf_none);
    double maxmass((m_Q-m_m1)/2.);   ///sqrt(2.));
    double sumwt(0.);
    
    for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
      if (vetodiquark && fdit->first.IsDiQuark()) continue;
      if (fdit->second->popweight>0. && fdit->second->massmin<maxmass) {
	sumwt += fdit->second->popweight;
      }
    }
    if (sumwt<=0) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   no flavour can be picked for gluon splitting: Q = "
		 <<m_Q<<", "<<"m_1 = "<<m_m1<<" "
		 <<"--> mass < "<<(sqrt((m_Q2-m_m1_2)/2.))<<"."<<std::endl;
      return false;
    }
    sumwt *= ran->Get();
    for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
      if (vetodiquark && fdit->first.IsDiQuark()) continue;
      if (fdit->second->popweight>0. && fdit->second->massmin<maxmass) {
	sumwt -= fdit->second->popweight;
      }
      if (sumwt<0.) {
	m_flav = fdit->first;
	break;
      }
    }
    if (m_flav.IsDiQuark()) m_flav = m_flav.Bar(); 
    if      (m_flav==Flavour(kf_d)) m_d++;
    else if (m_flav==Flavour(kf_u)) m_u++;
    else if (m_flav==Flavour(kf_s)) m_s++;
    m_tot++;
    
    m_m2  = m_m3 = m_flav.HadMass();
  }
  else {
    m_flav = Flavour(kf_gluon);
    m_m2   = m_flav.HadMass();
    m_m3   = p_split->m_flav.HadMass();
  }
  m_m2_2   = sqr(m_m2); 
  m_m3_2   = sqr(m_m3); 
  m_Qt2    = m_Q2-m_m1_2-m_m2_2-m_m3_2;
  return (m_flav!=Flavour(kf_none));
}


bool Splitting_Tools::
DetermineSplitting(Dipole * dip1,const bool & first,const bool & vetodiquark) {
  if (m_Q < m_m1) return false;
  int trials(0);
  while ((trials++)<1000) {
    if (ProduceKinematics(first,vetodiquark)) {
      if (m_analyse) {
	m_histograms[std::string("Splitting_Trials")]->Insert(trials);
      }
      //msg_Out()<<METHOD<<" yields kt = "<<m_kt<<" < "<<sqrt(m_kt2max)<<".\n";
      return true;
    }
  }
  return false;
}


bool Splitting_Tools::
ProduceKinematics(const bool & first,const bool & vetodiquark) {
  if (!SelectFlavour(vetodiquark)) return false;
  double s23min(m_m2_2+m_m3_2);
  double s23max(m_Qt2); //s23max(sqr(m_m2+m_m3)+m_pt2max*m_pt2max_factor);
  double ymax(1.-2.*m_m1*(m_Q-m_m1)/m_Qt2);
  if (ymax>(s23max-s23min)/m_Qt2) ymax = (s23max-s23min)/m_Qt2;
  double ymin(2.*m_m2*m_m3/m_Qt2);
  if (IsZero(ymin/ymax)) {
    ymin = sqr(m_m2+m_m3)/m_Qt2;
    if (ymin>ymax) ymin = ymax/10.;
  }
  double zmin,zmax;
  double OneMinExpo(-1),rand;
  double masscor(m_glusplit?1.:Max(m_m3_2/m_pt2max,1.));
  double exparg(m_pt2max_factor*m_pt2max/masscor);

  int trials(0);
  while ((trials++)<1000) {
    rand  = ATOOLS::ran->Get();
    if (OneMinExpo<=0.)  
      m_y = ymin * pow(ymax/ymin,rand);
    else 
      m_y = pow(pow(ymin,OneMinExpo) + 
		rand*(pow(ymax,OneMinExpo)-pow(ymin,OneMinExpo)),
		1./OneMinExpo); 
    m_s23 = m_y*m_Qt2+s23min;
    FixZRange(zmin,zmax);
    m_z   = p_kernels->SelectZ(zmin,zmax,1.,m_glusplit,m_leadsplit);
    m_kt2 = m_z*(1.-m_z)*m_s23-(1.-m_z)*m_m3_2-m_z*m_m2_2;
    double weight = m_kt2<0.||m_kt2>m_kt2max?0.:1.;
    weight *= (*p_as)(m_kt2,false)/p_as->MaxValue() * 
      (m_glusplit?m_kt2/m_s23<ran->Get():1.) * 
      (m_kt2>m_pt2max?exp(-sqr((m_kt2-m_pt2max)/exparg)):1.);
    if (weight<ran->Get()) continue;
    m_phi = 2.0*M_PI*ATOOLS::ran->Get();
    m_kt  = sqrt(m_kt2);
    if (ConstructKinematics()) {
      m_lastpt2 = m_kt2;
      break;
    }
    else {
      /*
	msg_Error()<<METHOD<<" "<<p_split->m_flav<<" -> "<<m_flav
	<<" ["<<p_split->m_info<<p_spect->m_info<<"]:"
	<<"crash for Q^2 = "<<m_Q2<<" --> "
	<<"s_23 = "<<m_s23<<", kt^2 = "<<m_kt2<<", z = "<<m_z
	<<" in ["<<zmin<<", "<<zmax<<"] from "<<m_m3_2
	<<"--> kt = "<<m_mom3.PPerp(m_mom2)<<"\n";
	exit(1);
      */
      return false;
    }
  }
  if (m_analyse) {
    m_histograms.find(std::string("Kinematics_Trials"))->second->Insert(trials);
  }
  return (trials<1000); 
}

bool Splitting_Tools::ConstructKinematics() {
  if (IsNan(m_kt) || IsNan(m_z) || IsNan(m_y)) {
    msg_Error()<<"Error in "<<METHOD
    	       <<"(kt = "<<m_kt<<", z = "<<m_z<<", y = "<<m_y<<").\n";
    return false;
  }
  double po(sqr(m_Q2-m_m23_2-m_m1_2)- 4.*m_m23_2*m_m1_2);
  double pn(sqr(m_Q2-m_s23-m_m1_2)  - 4.*m_s23*m_m1_2);
  if (po<0. || pn<0.) {
    msg_Error()<<"Error in "<<METHOD<<"(po,n = "<<po<<", "<<pn<<").\n";
    return false;
  }
  po = sqrt(po);
  pn = sqrt(pn);
  double gamma = (m_Q2-m_s23-m_m1_2+pn)/2.;

  double pref1 = (m_Q2-m_m23_2+m_m1_2)/(2.*m_Q2);
  double pref2 = (m_Q2-m_s23+m_m1_2)/(2.*m_Q2);
  Vec4D k1 = pn/po*m_mom1+(pref2-pn/po*pref1)*m_mom0;
  Vec4D k2 = m_mom0 - k1;
  Vec4D k3 = m_kt*sin(m_phi)*m_lperp;
  m_cms.BoostBack(k3);
  k3              += m_kt*cos(m_phi)*m_nperp;
  k3              += 
    m_z/pn*(gamma*k2-m_s23*k1)+
    (m_m3_2+m_kt2)/(m_z*pn)*(k1-m_m1_2/gamma*k2);

  m_mom1 = k1;
  m_mom3 = k3;
  m_mom2 = m_mom0-m_mom1-m_mom3;
  if (IsNan(m_mom2.Abs2())) {
    msg_Error()<<"Error in "<<METHOD<<"(kt2 = "<<m_kt2<<").\n"
	       <<m_mom1<<"+"<<m_mom2<<"+"<<m_mom3<<".\n";
    return false;
  }
  if (m_analyse) AnalyseKinematics(m_mom3,m_mom2,m_mom1);
  return true;
}


bool Splitting_Tools::
FixZRange(double & zmin,double & zmax) {
  //double m2_2(IsZero(m_m2_2/m_s23)?m_mmin_2:m_m2_2), m2(sqrt(m2_2));
  double m2_2(m_m2_2), m2(sqrt(m2_2));
  double disc = sqr(m_s23-m2_2-m_m3_2)-sqr(2.*m2*m_m3);
  double mean = m_s23+m_m3_2-m2_2;
  if (disc<0.) {
    zmin = m_m3_2/m_s23;
    zmax = 1.-zmin;
  }
  else {
    if (IsZero(m_m2_2/m_s23)) {
      zmin = (m_m3_2+m_m2_2)/m_s23;
      zmax = 1.-sqrt(m_mmin_2/m_s23); 
    }
    else {
      disc  = sqrt(disc);
      zmin  = (mean - disc)/(2.*m_s23);
      zmax  = Min((mean + disc)/(2.*m_s23),1.-0.001);
    }
    //if ((p_split->m_info=='L' || p_split->m_info=='B') &&
    //	zmax>0.5) zmin = Max(0.5,zmin); 
  }
  if (IsNan(zmin) || IsNan(zmax)) {
    zmin = 0.01;
    zmax = 1.-zmin;
  }
  if (zmin > zmax) {
    double temp = zmin;
    zmin = zmax;
    zmax = temp;
  }
  zmin = Max(0., zmin);
  zmax = Min(1., zmax);
  return true;
}

void Splitting_Tools::AftermathOfSplitting(Dipole * dip1) {
  double s12((m_mom1+m_mom2).Abs2());
  double s13((m_mom1+m_mom3).Abs2());
  double s23((m_mom2+m_mom3).Abs2()); 
  bool swap(false);
  if (p_split->m_flav.IsGluon()) {
    p_spect->m_mom    = m_mom1;
    if (p_spect->m_info=='L') {
      swap = s12/(s12+s13)>0.5; // ran->Get());
    }
    else if (p_spect->m_flav==Flavour(kf_gluon)) {
      swap = s12/(s12+s13)<0.5; // ran->Get());
    }
    p_out1 = new Proto_Particle(m_flav.Bar(),swap?m_mom3:m_mom2,'l');
    p_out2 = new Proto_Particle(m_flav,swap?m_mom2:m_mom3,'l');
    SetInfoTagsForOutgoings();
    p_out1->p_partner = p_out2;
    p_out2->p_partner = p_out1;
    p_out1->m_kt2max = p_out2->m_kt2max = Max(s23,m_Qt2);
    delete p_split;
  }
  else {
    p_split->m_mom = m_mom3;
    p_spect->m_mom = m_mom1;
    p_out1 = new Proto_Particle(Flavour(kf_gluon),m_mom2,'l');
    p_out2 = 0;
    p_out1->p_partner = p_split;
    p_split->m_kt2max = p_out1->m_kt2max = m_lastpt2;
  }
}

void Splitting_Tools::SetInfoTagsForOutgoings() const {
  if (p_split->m_info=='B') { 
    p_out1->m_info = p_out2->m_info = p_split->m_info;
  }
  else {
    switch (m_leading) {
    case leading::quarks_and_gluons:
      if (p_split->m_flav.IsGluon()) { 
	if (p_split->m_info=='L') {
	  if (m_z<0.5) 
	    p_out1->m_info=p_split->m_info;
	  else 
	    p_out2->m_info=p_split->m_info;
	}
      }
      break;
    case leading::quarks_and_gluons2:
      if (p_split->m_flav.IsGluon() && 
	  (p_split->m_info=='L' || p_split->m_info=='B')) {
	p_out1->m_info=p_split->m_info;
	p_out2->m_info=p_split->m_info;
      }
      break;
    case leading::only_quarks:
    case leading::none:
    default:
      break;
    }
  }
}






void Splitting_Tools::
AnalyseKinematics(const Vec4D & q1,const Vec4D & q2,const Vec4D & q3) {
  double Ebeam = rpa->gen.Ecms()/2.;
  if (p_split->m_flav.IsGluon()) {
    m_histograms[std::string("PT_Gluon_Splitting")]->Insert(m_kt);
    m_histograms[std::string("PT2_Gluon_Splitting")]->Insert(m_kt*m_kt);
    m_histograms[std::string("Z_Gluon_Splitting")]->Insert(m_z);
    if (p_spect->m_flav.Kfcode()==5) {
      m_histograms[std::string("XB_bquark_before_gsplit")]->Insert(p_spect->m_mom[0]/Ebeam);
      m_histograms[std::string("XB_bquark_after_gsplit")]->Insert(q1[0]/Ebeam);
    }
  }
  else {
    if (p_split->m_info=='L' && p_spect->m_info=='L' &&
	p_split->m_flav.Kfcode()==5) {
      m_histograms[std::string("PT_Gluon_First")]->Insert(m_kt);
      m_histograms[std::string("Z_Gluon_First")]->Insert(m_z);
    }
    else if (p_split->m_info=='L' &&
	     p_split->m_flav.Kfcode()==5) {
      m_histograms[std::string("PT_Gluon_Leading")]->Insert(m_kt);
      m_histograms[std::string("Z_Gluon_Leading")]->Insert(m_z);
    }
    else {
      m_histograms[std::string("PT_Gluon_Emission")]->Insert(m_kt);
      m_histograms[std::string("Z_Gluon_Emission")]->Insert(m_z);
    }
    m_histograms[std::string("PT2_Gluon_Emission")]->Insert(m_kt*m_kt);
    if (p_spect->m_flav.Kfcode()==5) {
      m_histograms[std::string("XB_bquark_before_qsplit")]->Insert(p_spect->m_mom[0]/Ebeam);
      m_histograms[std::string("XB_bquark_after_qsplit")]->Insert(q1[0]/Ebeam);
    }
  }
}
