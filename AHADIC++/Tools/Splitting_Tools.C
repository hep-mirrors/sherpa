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


Splitting_Tools::Splitting_Tools(const leading::code & lead,const PTOrder::code & ptorder,
				 const ZForm::code & zform,
				 Strong_Coupling * as,const bool & analyse) :
  m_leading(lead), m_ptorder(ptorder),
  m_fourquarks(false), m_analyse(true),
  m_masstreatment(int(hadpars.Get(std::string("Mass_treatment")))), 
  p_as(as), p_kernels(new Splitting_Functions(zform,m_masstreatment)), p_options(NULL),
  p_spect(NULL), p_split(NULL), p_out1(NULL), p_out2(NULL),
  m_mmin_2(sqr(hadpars.GetConstituents()->MinMass())),
  m_pt2max(sqr(hadpars.Get(std::string("ptmax")))), 
  m_pt2max_factor(sqr(hadpars.Get(std::string("ptmax_factor")))), 
  m_pt02(dabs(hadpars.Get(std::string("pt02")))), 
  m_pt2min(dabs(hadpars.Get(std::string("pt2min")))),
  m_lastpt2(-1.), 
  m_tot(0),m_d(0),m_s(0),m_u(0),m_reject_y(0),m_reject_z(0)
{ 
  if (m_analyse) {
    m_histograms[std::string("Splitting_Trials")]        = new Histogram(0,0.,1000,500);
    m_histograms[std::string("Enforced_Trials")]         = new Histogram(0,0.,1000,500);
    m_histograms[std::string("Kinematics_Trials")]       = new Histogram(0,0.,1000,500);
    m_histograms[std::string("PT_Gluon_Splitting")]      = new Histogram(0,0.,100.,2000);
    m_histograms[std::string("PT_Gluon_Emission")]       = new Histogram(0,0.,100.,2000);
    m_histograms[std::string("PT_Gluon_Primary")]        = new Histogram(0,0.,250.,1000);
    m_histograms[std::string("PT2_Gluon_Splitting")]     = new Histogram(0,0.,10000.,1000);
    m_histograms[std::string("PT2_Gluon_Emission")]      = new Histogram(0,0.,10000.,1000);
    m_histograms[std::string("Z_Gluon_Splitting")]       = new Histogram(0,0.,1.,50);
    m_histograms[std::string("Z_Gluon_Emission")]        = new Histogram(0,0.,1.,50);
    m_histograms[std::string("Z_Gluon_Primary")]         = new Histogram(0,0.,1.,50);
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
      name  = std::string("Fragmentation_Analysis/")+hit->first+std::string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
  //msg_Debugging()<<"Flavour production in "<<m_tot<<" gluon splittings (tools): "
  //	   <<"s_rate = "<<double(m_s)/double(m_tot)<<", "
  //	   <<"u_rate = "<<double(m_u)/double(m_tot)<<", "
  //	   <<"d_rate = "<<double(m_d)/double(m_tot)<<"."<<std::endl;
  delete p_kernels;
}





void Splitting_Tools::SetSpectatorAndSplitter(Dipole * dip) {
  p_split=p_spect=NULL;
  if (dip->Triplet()->m_flav.IsGluon() && 
      (dip->AntiTriplet()->m_flav.IsQuark() ||
       dip->AntiTriplet()->m_flav.IsDiQuark())) {
    // asymmetric case: triplet splits (gluon)
    p_spect = dip->AntiTriplet(); 
    p_split = dip->Triplet();
    dip->SetSwitched(true);
  }
  else if (dip->AntiTriplet()->m_flav.IsGluon() && 
	   (dip->Triplet()->m_flav.IsQuark() ||
	    dip->Triplet()->m_flav.IsDiQuark())) {
    // asymmetric case: antitriplet splits (gluon)
    p_split = dip->AntiTriplet(); 
    p_spect = dip->Triplet();
    dip->SetSwitched(false);
  }
  else if ((dip->Triplet()->m_flav.IsGluon() && 
	    dip->AntiTriplet()->m_flav.IsGluon())     ||
	   (!dip->Triplet()->m_flav.IsGluon() && 
	    !dip->AntiTriplet()->m_flav.IsGluon())) { 
    if ((dip->Triplet()->m_info=='L' && 
	 dip->AntiTriplet()->m_info=='L') ||
	(dip->Triplet()->m_info=='l' && 
	 dip->AntiTriplet()->m_info=='l') ||
	(dip->Triplet()->m_info=='B' && 
	 dip->AntiTriplet()->m_info=='B')) {
      // symmetric case: pick at random.
      double l1(sqr(Max(dip->Triplet()->m_mass,1.e-4)));
      double l2(sqr(Max(dip->AntiTriplet()->m_mass,1.e-4)));
      double disc = l1/(l1+l2);
      if (ran.Get()>disc) {
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
    else if ((dip->Triplet()->m_info=='L' || dip->Triplet()->m_info=='B') && 
	     dip->AntiTriplet()->m_info!='L') {
      // asymmetric case: antitriplet splits (non-leading with leading/beam remnant spectator)
      p_split = dip->AntiTriplet(); 
      p_spect = dip->Triplet();
      dip->SetSwitched(false);
    }
    else if (dip->Triplet()->m_info!='L' && 
	     (dip->AntiTriplet()->m_info=='L' || dip->AntiTriplet()->m_info=='B')) {
      // asymmetric case: antitriplet splits (non-leading with leading spectator)
      p_split = dip->Triplet(); 
      p_spect = dip->AntiTriplet();
      dip->SetSwitched(true);
    }	
  }
  if (p_split->m_kt2max<=0.) 
    p_split->m_kt2max = p_spect->m_kt2max = ATOOLS::Max(p_split->m_kt2max,p_spect->m_kt2max);
}

void Splitting_Tools::SwapSpectatorAndSplitter(Dipole * dip) {
  Proto_Particle * help(p_split); p_split = p_spect; p_spect = help;
  if (dip->IsSwitched()) dip->SetSwitched(false);
                    else dip->SetSwitched(true); 
  if (p_split->m_kt2max<=0.) 
    p_split->m_kt2max = p_spect->m_kt2max = ATOOLS::Max(p_split->m_kt2max,p_spect->m_kt2max);
}

bool Splitting_Tools::PrepareKinematics(Dipole * dip,const bool & first,const bool & enforce) {
  if (dip->MassBar2()<4.*m_mmin_2) {
    msg_Debugging()<<METHOD<<"(massmin^2 = "<<m_mmin_2<<", red.mass^2 = "<<dip->MassBar2()<<") :"<<std::endl
		   <<"   --> Dipole cannot split."<<std::endl;
    //(*dip).Output();
    return false;
  }

  m_mom1 = p_spect->m_mom;
  m_mom3 = p_split->m_mom;
  m_mom0 = m_mom1+m_mom3;
  m_Q2   = m_mom0.Abs2();
  m_Q    = sqrt(m_Q2);

  m_m1    = p_spect->m_flav.HadMass();
  m_m1_2  = sqr(m_m1);
  m_m23_2 = sqr(p_split->m_flav.HadMass());

  m_cms  = Poincare(m_mom0);
  m_cms.Boost(m_mom0);
  m_cms.Boost(m_mom1);
  m_cms.Boost(m_mom3);

  m_zrot = Poincare(m_mom1,Vec4D::ZVEC);
  m_zrot.Rotate(m_mom1);
  m_zrot.Rotate(m_mom3);

  m_glusplit = p_split->m_flav.IsGluon();

  msg_Debugging()<<METHOD<<"(glusplit = "<<m_glusplit<<", Q = "<<m_Q<<"): " 
		 <<"can in principle split."<<std::endl;
  //(*dip).Output();
  return true;
}



bool Splitting_Tools::SelectFlavour(const bool & vetodiquark)
{
  if (m_glusplit) {
    m_flav = Flavour(kf_none);
    // from QT^2>0:
    double maxmass((m_Q-m_m1)/sqrt(2.));
    double sumwt(0.);
    
    for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
      if (vetodiquark && fdit->first.IsDiQuark()) continue;
      if (fdit->second->popweight>0. && fdit->second->massmin<maxmass) {
	sumwt += fdit->second->popweight;
      }
    }
    if (sumwt<=0) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   no flavour can be picked for gluon splitting: Q = "<<m_Q<<", "
		 <<"m_1 = "<<m_m1<<" --> mass < "<<(sqrt((m_Q2-m_m1_2)/2.))<<"."<<std::endl;
      return false;
    }
    sumwt *= ran.Get();
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
    if (m_flav==Flavour(kf_d)) m_d++;
    else if (m_flav==Flavour(kf_u)) m_u++;
    else if (m_flav==Flavour(kf_s)) m_s++;
    m_tot++;
    
    //std::cout<<"-------------------------------------------------"<<std::endl;
    m_m2  = m_m3 = m_flav.HadMass();
  }
  else {
    m_flav = Flavour(kf_gluon);
    m_m2   = m_flav.HadMass();
    m_m3   = p_split->m_flav.HadMass();
  }
  m_m2_2   = sqr(m_m2); 
  m_m3_2   = sqr(m_m3); 
  m_Qt2    = m_Q2 - m_m1_2 - m_m2_2 - m_m3_2;

  m_ymin   = 2.*m_m2*m_m3/m_Qt2;
  m_ymax   = (m_Qt2-2.*m_m1*(m_Q-m_m1))/m_Qt2;
  return (m_flav!=Flavour(kf_none));
}

bool Splitting_Tools::DetermineSplitting(Dipole * dip1,const bool & first,const bool & vetodiquark) {
  int trials(0);
  msg_Debugging()<<"In "<<METHOD<<"(first = "<<first<<", veto = "<<vetodiquark<<")."<<std::endl;
  //dip1->Output();
  while (trials<1000) {
    trials++;
    if (ProduceKinematics(first,vetodiquark) &&
	ConstructKinematics()) break;
  }
  if (m_analyse) {
    Histogram * histo = (m_histograms.find(std::string("Splitting_Trials")))->second;
    histo->Insert(trials);
  }
  msg_Debugging()<<"In "<<METHOD<<"(first = "<<first<<", veto = "<<vetodiquark<<") -> succeeded."<<std::endl;
  return true;
}


bool Splitting_Tools::ProduceKinematics(const bool & first,const bool & vetodiquark) {
  msg_Debugging()<<"In "<<METHOD<<"."<<std::endl;
  if (!SelectFlavour(vetodiquark)) {
    msg_Error()<<"no flavour selected, return false."<<std::endl;
    return false;
  }
  msg_Debugging()<<"   flav = "<<m_flav<<"."<<std::endl;
  if (!FixRanges(first)) return false;
  msg_Debugging()<<"   kt2min = "<<m_kt2_min<<", "
		 <<"z in ["<<m_zmin<<", "<<m_zmax<<"]."<<std::endl;

  m_kt2 = m_z = -1.;
  int trials(0);
  while ((m_kt2<0. || m_z<0.) && trials<1000) {
    trials++;
    m_z   = p_kernels->SelectZ(m_zmin,m_zmax,m_glusplit,p_split->m_info=='L');
    if (!UpdateRanges(first)) {
      m_kt2 = -1.;
      continue;
    }
    m_kt2 = p_as->SelectPT(m_kt2_max,m_kt2_min);
    if (!KinCheck()) {
      m_kt2 = -1;
      continue;
    }
    if (p_kernels->Weight(m_kt2,m_z,p_split->m_flav.IsGluon(),
			  p_split->m_info=='L' && p_spect->m_info!='L')<ran.Get()) {
      m_kt2 = -1;
      continue;
    }
  }
  m_phi = 2.*M_PI*ran.Get();
  msg_Debugging()<<"___ "<<METHOD<<" : exit inner loop with kt2 = "<<m_kt2<<" "
	   <<"and z = "<<m_z<<" for "<<m_flav<<"."<<std::endl;

  if (m_analyse) {
    Histogram * histo = (m_histograms.find(std::string("Kinematics_Trials")))->second;
    histo->Insert(trials);
  }
  return (trials<1000);
}

bool Splitting_Tools::FixRanges(const bool & first) {
  double a(m_Qt2+m_m2_2+m_m3_2), b(m_Qt2+2.*m_m3_2);
  m_kt2_min = m_pt2min;
  if ((p_spect->m_info=='L') && !p_spect->m_flav.IsGluon() && m_m1_2>1.e-4) {
    m_kt2_min *= (first?4.:1.)*m_pt2min/m_m1_2;
  }
  if (!m_glusplit && m_m3_2>0.) {
    m_kt2_min *= (first?4.:1.)*m_pt2min/m_m3_2;
  }
  double pt2check(sqr(b)/(4.*a)-m_m3_2);
  if (m_kt2_min>pt2check) {
    m_kt2_min = m_pt2min;
    if (m_kt2_min>pt2check) {
      m_kt2_min = pt2check/4.;
    }
  }
  double c(m_kt2_min+m_m3_2);
  double disc(b*b-4.*a*c);
  if (disc<0.) {
    return false;
    msg_Error()<<"ERROR in "<<METHOD<<"(first = "<<first<<"):"<<std::endl
	       <<"   Discriminator for z range negative.  This should not happen, "
	       <<"kt2_min = "<<m_kt2_min<<" from "<<pt2check<<"."<<std::endl;
    disc = b;
  }
  m_zmin = Max(0.,(b-sqrt(disc))/(2.*a));
  m_zmax = Min(1.,(b+sqrt(disc))/(2.*a));
  //if (first && m_zmax-m_zmin>2.*m_kt2_min/m_m1_2) {
  //  m_zmin += m_kt2_min/m_m1_2;
  //  m_zmax -= m_kt2_min/m_m1_2;
  // }
  return true;
}

bool Splitting_Tools::UpdateRanges(const bool & first) {
  m_kt2_max = (sqr(m_Q-m_m1)-m_m2_2-m_m3_2)*m_z*(1.-m_z)-m_m2_2*m_z*m_z-m_m3_2*(1.-m_z)*(1.-m_z);
  msg_Debugging()<<"... "<<METHOD<<"(z = "<<m_z<<") --> kt^2_max = "<<m_kt2_max<<"."<<std::endl;
  if (m_kt2_max<0.) {
    return false;
  }
  if (m_ptorder==PTOrder::flat) {
    if (m_kt2_max>m_pt2max) m_kt2_max = sqrt(m_kt2_max*m_pt2max);
  }
  else {
    bool order((m_ptorder==PTOrder::gluon_split && p_split->m_flav==Flavour(kf_gluon)) ||
	       (m_ptorder==PTOrder::gluon_emit && p_split->m_flav!=Flavour(kf_gluon)) ||
	       (m_ptorder==PTOrder::total));
    if (!first && order) m_kt2_max = Min(m_kt2_max,p_split->m_kt2max);
    else {
      if ((p_spect->m_info=='B')) {
	if (p_split->m_info=='B') 
	  m_kt2_max = Min(m_kt2_max,m_pt2max_factor*sqrt(4.*p_spect->m_mom.PPerp2()*p_split->m_mom.PPerp2()));
	else 
	  m_kt2_max = Min(m_kt2_max,m_pt2max_factor*sqrt(4.*p_split->m_mom.PPerp2()*m_pt2min));
      }
      else {
	m_kt2_max = Min(m_kt2_max,m_pt2max_factor*p_spect->m_mom.PPerp2(p_split->m_mom));
	
      }
    }
    msg_Tracking()<<"... first: "<<p_split->m_info<<p_spect->m_info
		  <<"("<<p_split->m_flav<<" "<<p_spect->m_flav<<"): "
		  <<"kt_max = "<<sqrt(m_kt2_max)<<" from pt(kin) = "<<p_spect->m_mom.PPerp2(p_split->m_mom)<<", "
		  <<"p_min = "<<Min(p_split->m_mom.PSpat(),p_spect->m_mom.PSpat())<<" and "
		  <<"cos = "<<p_split->m_mom.CosTheta(p_spect->m_mom)
		  <<", red. mass = "<<sqrt(m_Qt2)<<"."<<std::endl;
    if (m_kt2_max>m_pt2max) m_kt2_max = sqrt(m_kt2_max*m_pt2max);
    m_kt2_min  = sqr(m_z*m_m2+(1.-m_z)*m_m3); 
    m_kt2_min += m_pt2min;
    if (!p_spect->m_flav.IsGluon() && m_m1_2>1.e-4) {
      m_kt2_min *= (first?4.:1.)*m_pt2min/m_m1_2;
    }
    if (!m_glusplit && m_m3_2>0.) {
      m_kt2_min *= (first?4.:1.)*m_pt2min/m_m3_2;
    }
    if (m_kt2_min>m_kt2_max) {
      m_kt2_min = m_kt2_max/4.;
      msg_Debugging()<<"   correct to kt^2_min,max = "<<m_kt2_min<<", "<<m_kt2_max<<"."<<std::endl;
    }
  }
  return true;
}
    









void Splitting_Tools::AftermathOfSplitting(Dipole * dip1) {
  if (p_split->m_flav.IsGluon()) {
    p_spect->m_mom    = m_mom1;
    if (p_spect->m_info=='L' && !p_spect->m_flav.IsGluon()) {
      double s12((m_mom1+m_mom2).Abs2()), s13((m_mom1+m_mom3).Abs2());
      bool swap(s12/(s12+s13)>0.5); // ran.Get());
      if ((p_spect->m_flav.IsQuark() && !p_spect->m_flav.IsAnti()) || 
	  (p_spect->m_flav.IsDiQuark() && p_spect->m_flav.IsAnti())) {
	p_out1 = new Proto_Particle(m_flav.Bar(),swap?m_mom3:m_mom2,'l');
	p_out2 = new Proto_Particle(m_flav,swap?m_mom2:m_mom3,'l');
      }
      else {
	p_out1 = new Proto_Particle(m_flav.Bar(),swap?m_mom2:m_mom3,'l');
	p_out2 = new Proto_Particle(m_flav,swap?m_mom3:m_mom2,'l');
      }
    }
    else {
      p_out1 = new Proto_Particle(m_flav.Bar(),m_mom2,'l');
      p_out2 = new Proto_Particle(m_flav,m_mom3,'l');
    }
    SetInfoTagsForOutgoings();
    p_out1->p_partner = p_out2;
    p_out2->p_partner = p_out1;
    p_out1->m_kt2max = p_out2->m_kt2max = m_lastpt2;
    delete p_split;
    if (false && !p_spect->m_flav.IsGluon()) {
      if ((p_spect->m_flav.IsQuark() && !p_spect->m_flav.IsAnti()) ||
	  (p_spect->m_flav.IsDiQuark() && p_spect->m_flav.IsAnti())) {
	if ((m_mom2+m_mom1).Abs2()>(m_mom3+m_mom1).Abs2()) {
	  p_out1->m_mom = m_mom3;
	  p_out2->m_mom = m_mom2;
	}
      }
      else {
	if ((m_mom2+m_mom1).Abs2()<(m_mom3+m_mom1).Abs2()) {
	  p_out1->m_mom = m_mom3;
	  p_out2->m_mom = m_mom2;
	}
      }
    }
    msg_Debugging()<<"___ "<<METHOD<<" yields new particles : "<<std::endl
		  <<(*p_out1)<<(*p_out2)<<std::endl;
  }
  else {
    p_split->m_mom = m_mom3;
    p_spect->m_mom = m_mom1;
    p_out1 = new Proto_Particle(ATOOLS::Flavour(kf_gluon),m_mom2,'l');
    p_out2 = 0;
    p_out1->p_partner = p_split;
    p_split->m_kt2max = p_out1->m_kt2max = m_lastpt2;
  }
}

void Splitting_Tools::SetInfoTagsForOutgoings() const {
  switch (m_leading) {
  case leading::quarks_and_gluons:
    if (p_split->m_flav.IsGluon() && p_split->m_info=='L') {
      if (m_z<0.5) 
	p_out1->m_info='L';
      else 
	p_out2->m_info='L';
    }
    break;
  case leading::quarks_and_gluons2:
    if (p_split->m_flav.IsGluon() && p_split->m_info=='L') {
      p_out1->m_info='L';
      p_out2->m_info='L';
    }
    break;
  case leading::only_quarks:
  case leading::none:
  default:
    break;
  }
}






bool Splitting_Tools::KinCheck() {
  msg_Debugging()<<METHOD<<"(kt^2 = "<<m_kt2<<", z = "<<m_z<<")."<<std::endl;
  //double mu1(m_m1_2/m_Q2),mu2(m_m2_2/m_Q2),mu3(m_m3_2/m_Q2);
  if (m_masstreatment==3 && p_split->m_flav.IsGluon()) {
    m_s23 = m_kt2;
    m_kt2 = m_Qt2*m_s23-(1.-m_z)/m_z*m_m2_2-m_z/(1.-m_z)*m_m3_2;
    if (m_kt2<0.0) return false;
    m_y   = (m_s23-m_m2_2-m_m3_2)/m_Qt2;
  }
  else { 
    m_y   = (m_kt2 + sqr(1.-m_z)*m_m2_2 + sqr(m_z)*m_m3_2)/(m_z*(1.-m_z)*m_Qt2);
    m_s23 = m_y*m_Qt2+m_m2_2+m_m3_2;
    if (m_s23<0.0) {
      msg_Debugging()<<"   ... negative s23."<<std::endl;
      return false;
    }
    if (m_masstreatment==3 && p_split->m_flav.IsGluon()) {
      if (m_kt2/m_s23*(*p_as)(m_s23)/(*p_as)(m_kt2)<ran.Get()) {
	msg_Debugging()<<"   ... reweighted as."<<std::endl;
	return false;
      }
    }
  }
  if (m_y<m_ymin || m_y>m_ymax) {
    msg_Debugging()<<"   ... y = "<<m_y<<" out of bounds ["<<m_ymin<<", "<<m_ymax<<"]."<<std::endl;
    m_reject_y++;
    return false;
  }
  double viji = sqrt(sqr(m_Qt2*m_y)-4.*m_m2_2*m_m3_2)/(m_Qt2*m_y+2.*m_m3_2);
  double vijk = sqrt(sqr(2.*m_m1_2+m_Qt2*(1.-m_y))-4.*m_m1_2)/(m_Qt2*(1.-m_y));
  double frac = (2.*m_m3_2+m_Qt2*m_y)/(2.*(m_m3_2+m_m2_2+m_Qt2*m_y));
  double zm   = frac*(1.-viji*vijk);
  double zp   = frac*(1.+viji*vijk);
  if (m_z<zm || m_z>zp) {
    msg_Debugging()<<"   ... z = "<<m_z<<" out of (updated) bounds ["<<zm<<", "<<zp<<"]."<<std::endl;
    m_reject_z++;
    return false;
  }
  p_kernels->SetScales(m_m1_2,m_m2_2,m_m3_2,m_Q2,m_y,
		       viji,vijk,zm,zp,m_y*m_Qt2/2.);

  return true;
}




bool Splitting_Tools::ConstructKinematics() {
  if (IsZero(m_Qt2))  
    msg_Error()<<"Warning in "<<METHOD<<": Qt^2 = "<<m_Qt2<<" for kt^2 = "<<m_pt2<<"."<<std::endl;

  m_lastpt2 = m_kt2;

  Vec4D q1(0.,0.,0.,0.),q2(0.,0.,0.,0.),q3(0.,0.,0.,0.);
  Vec4D nperp = Vec4D(0.,cos(m_phi)*s_ex + sin(m_phi)*s_ey);
  if (m_m1_2==0.) {
    if (m_kt2<0.) {
      msg_Error()<<"Warning in "<<METHOD<<"(massless): negative kt^2 = "<<m_kt2<<"."<<std::endl;
      return false;
    }
    m_kt   = sqrt(m_kt2);
    double fac1 = m_Q2 - ((1.-m_z)*m_m2_2+m_z*m_m3_2+m_kt2)/(m_z*(1.-m_z));
    if (fac1<0.) {
      msg_Debugging()<<"      ----> (fail0) for kt^2 = "<<m_kt2<<", z = "<<m_z<<std::endl;
      return false;
    }
    double div  = m_Q2-m_m23_2;
    if (IsZero(div)) {
      msg_Error()<<"Warning in "<<METHOD<<": div = "<<div<<" for kt^2 = "<<m_kt2<<"."<<std::endl;
      return false;
    }
    fac1 /= div;

    q1 =                                                             fac1 * m_mom1;
    q3 =      m_z * m_mom3 +    (m_kt2+m_m2_2-sqr(m_z)*m_m23_2)/(m_z*div) * m_mom1 + m_kt * nperp;
  }
  else {
    double po  = sqr(m_Q2-m_m23_2-m_m1_2)-4.*m_m23_2*m_m1_2; 
    double pn  = sqr(m_Q2-m_s23-m_m1_2)-4.0*m_s23*m_m1_2;
    if (po<0.0 || pn<0.0) {
      msg_Debugging()<<"      ----> (fail1) for kt^2 = "<<m_kt2<<", z = "<<m_z<<std::endl;
      return false;
    }
    if (IsZero(po) || IsZero(pn)) {
      msg_Error()<<"Warning in "<<METHOD<<": po, pn = "<<po<<", "<<pn
		 <<" for kt^2 = "<<m_kt2<<"."<<std::endl;
      return false;
    }
    q1 = 
      (m_Q2+m_m1_2-m_s23)/(2.*m_Q2)*m_mom0+
      sqrt(pn/po) * (m_mom1 - (m_Q2+m_m1_2-m_m23_2)/(2.*m_Q2) * m_mom0);
    
    Vec4D  q23 = m_mom0-q1;
    if (sqr(q23*q1)-m_s23*m_m1_2<0.) {
      //msg_Error()<<"Warning in "<<METHOD<<": gamma explodes: "<<(sqr(q23*q1)-m_s23*m_m1_2)
      //		 <<" for kt^2 = "<<m_kt2<<"."<<std::endl;
      return false;
    }
    double gam = q23*q1+sqrt(sqr(q23*q1)-m_s23*m_m1_2);
    if (IsZero(m_s23*m_m1_2-gam*gam) || IsZero(gam)) {
      msg_Error()<<"Warning in "<<METHOD<<": beta explodes: "<<(m_s23*m_m1_2-gam*gam)
		 <<" for kt^2 = "<<m_kt2<<" & gamma = "<<gam<<"."<<std::endl;
      return false;
    }
    double a23 = m_s23/gam, a1=m_m1_2/gam, beta=1.0/(1.0-a23*a1);
    Vec4D  l   = beta*(q23-a23*q1);
    Vec4D  n   = beta*(q1-a1*q23);
    if (IsZero(a1) || IsZero(1.-m_y) && IsZero(gam) || IsZero((1.+a23*a1)) ||
	IsZero((1./a1-(m_m2_2+m_m3_2)/gam))) {
      msg_Error()<<"Warning in "<<METHOD<<": z explodes: "<<std::endl
		 <<"  a23 = "<<a23<<", a1 = "<<a1<<", y "<<m_y<<" gamma = "<<gam
		 <<" & masses = "<<m_m2_2<<"/"<<m_m3_2<<"."<<std::endl;
      return false;
    }
    m_z =
      (m_z/a1-m_y/(1.-m_y)-2.*m_m3_2/gam/(1.+a23*a1))/
      ((1./a1-(m_m2_2+m_m3_2)/gam)/(1.0+a23*a1)-m_y/(1.0-m_y));
    m_kt2 = gam*m_y/(1.-m_y)*(1.+a23*a1)*m_z*(1.-m_z)-sqr(1.-m_z)*m_m3_2-m_z*m_z*m_m2_2;
    if (m_kt2<0.0 || m_z<0. || m_z>1. || IsNan(m_kt2) || IsNan(m_z)) {
      msg_Debugging()<<"      ----> (fail2) for kt^2 = "<<m_kt2<<", z = "<<m_z<<std::endl;
      return false;
    }
    m_kt = sqrt(m_kt2);
    //    q1 = zt*l + (mi2+ktt*ktt)/(gam*zt)*n +
    //  + ktt*cos(phi)*Vec4D(0.0,1.0,0.0,0.0) 
    //  + ktt*sin(phi)*Vec4D(0.0,0.0,1.0,0.0); 
    q3 = m_z*l + (m_m3_2+m_kt2)/(gam*m_z)*n + m_kt*nperp;
  }
  q2 = m_mom0-q3-q1;

  msg_Debugging()<<"___ "<<METHOD<<": before boosting back."<<std::endl
		 <<"___ "<<m_mom0<<" = "<<m_mom1<<" ("<<sqrt(Max(0.,m_mom1.Abs2()))<<")"
		 <<" + "<<m_mom3<<" ("<<sqrt(Max(0.,m_mom3.Abs2()))<<")"<<std::endl
		 <<"___ "<<q1<<" ("<<sqrt(Max(0.,q1.Abs2()))<<")"
		 <<" + "<<q2<<" ("<<sqrt(Max(0.,q2.Abs2()))<<")"
		 <<" + "<<q3<<" ("<<sqrt(Max(0.,q3.Abs2()))<<")."<<std::endl;

  if (m_m1_2==0.) {
    msg_Debugging()<<"___ "<<METHOD<<" summary: "<<std::endl
		  <<"     Masses: m_m2 = m_m3 = "<<m_m2<<" = "<<m_m3<<" --> s23 = "<<m_m23_2<<"."<<std::endl;
  }

  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Debugging()<<"      ----> (fail3) for kt^2 = "<<m_kt2<<", z = "<<m_z<<std::endl;
    return false;
  }
  

  //std::cout<<"====   Splitting "<<p_split->m_flav<<" with params : "
  //	   <<"kt2 = "<<m_kt2<<", z = "<<m_z
  //	   <<"      ---> q2 (cms) = "<<q2<<"."<<std::endl;

  m_zrot.RotateBack(q1);
  m_zrot.RotateBack(q2);
  m_zrot.RotateBack(q3);
  m_cms.BoostBack(q1);
  m_cms.BoostBack(q2);
  m_cms.BoostBack(q3);

  if (m_analyse) AnalyseKinematics(q1,q2,q3);

  m_mom1 = q1;
  m_mom2 = q2;
  m_mom3 = q3;
  m_cms.BoostBack(m_mom0);

  msg_Debugging()<<"___ "<<METHOD<<"(out): succeeded to build kinematics "
	   <<"for kt^2 = "<<m_kt2<<", z = "<<m_z<<" for "<<m_flav<<"."<<std::endl;

  return true;
}

void Splitting_Tools::AnalyseKinematics(const ATOOLS::Vec4D & q1,
					const ATOOLS::Vec4D & q2,
					const ATOOLS::Vec4D & q3) {
  double Ebeam = rpa.gen.Ecms()/2.;
  Histogram * histo;
  if (p_split->m_flav.IsGluon()) {
    histo = (m_histograms.find(std::string("PT_Gluon_Splitting")))->second;
    histo->Insert(m_kt);
    histo = (m_histograms.find(std::string("PT2_Gluon_Splitting")))->second;
    histo->Insert(m_kt*m_kt);
    histo = (m_histograms.find(std::string("Z_Gluon_Splitting")))->second;
    histo->Insert(m_z);
    if (p_spect->m_flav.Kfcode()==5) {
      //std::cout<<"In "<<METHOD<<" for spectator = b, splitter = g:"<<std::endl
      //	       <<"   Ebeam = "<<Ebeam<<", Eb = "<<q1[0]<<" / "<<p_spect->m_mom[0]<<"."<<std::endl;
      histo = (m_histograms.find(std::string("XB_bquark_before_gsplit")))->second;
      histo->Insert(p_spect->m_mom[0]/Ebeam);
      histo = (m_histograms.find(std::string("XB_bquark_after_gsplit")))->second;
      histo->Insert(q1[0]/Ebeam);
    }
  }
  else {
    if (p_split->m_info=='L') {
      histo = (m_histograms.find(std::string("PT_Gluon_Primary")))->second;
      histo->Insert(m_kt);
      histo = (m_histograms.find(std::string("Z_Gluon_Primary")))->second;
      histo->Insert(m_z);
    }
    else {
      histo = (m_histograms.find(std::string("PT_Gluon_Emission")))->second;
      histo->Insert(m_kt);
      histo = (m_histograms.find(std::string("Z_Gluon_Emission")))->second;
      histo->Insert(m_z);
    }
    histo = (m_histograms.find(std::string("PT2_Gluon_Emission")))->second;
    histo->Insert(m_kt*m_kt);
    if (p_spect->m_flav.Kfcode()==5) {
      //std::cout<<"In "<<METHOD<<" for spectator = b, splitter = q ("<<p_split->m_flav<<"):"<<std::endl
      //	       <<"   Ebeam = "<<Ebeam<<", Eb = "<<q1[0]<<" / "<<p_spect->m_mom[0]<<"."<<std::endl;
      histo = (m_histograms.find(std::string("XB_bquark_before_qsplit")))->second;
      histo->Insert(p_spect->m_mom[0]/Ebeam);
      histo = (m_histograms.find(std::string("XB_bquark_after_qsplit")))->second;
      histo->Insert(q1[0]/Ebeam);
    }
  }
}
