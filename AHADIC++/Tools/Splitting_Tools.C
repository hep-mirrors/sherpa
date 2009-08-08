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
  m_lastpt2(-1.), 
  p_as(as), p_kernels(new Splitting_Functions(zform,m_masstreatment)), p_options(NULL),
  p_spect(NULL), p_split(NULL), p_out1(NULL), p_out2(NULL),
  m_mmin_2(sqr(hadpars.GetConstituents()->MinMass())),
  m_pt2max(sqr(hadpars.Get(std::string("ptmax")))), 
  m_pt2max_factor(sqr(hadpars.Get(std::string("ptmax_factor")))), 
  m_pt02(dabs(hadpars.Get(std::string("pt02")))), 
  m_pt2min(dabs(hadpars.Get(std::string("pt2min")))),
  m_tot(0),m_d(0),m_s(0),m_u(0)
{ 
  if (m_analyse) {
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
  //msg_Out()<<"Flavour production in "<<m_tot<<" gluon splittings (tools): "
  //	   <<"s_rate = "<<double(m_s)/double(m_tot)<<", "
  //	   <<"u_rate = "<<double(m_u)/double(m_tot)<<", "
  //	   <<"d_rate = "<<double(m_d)/double(m_tot)<<"."<<std::endl;
  delete p_kernels;
}





double Splitting_Tools::SetSpectatorAndSplitter(Dipole * dip) {
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

  return p_split->m_kt2max;
}

double Splitting_Tools::SwapSpectatorAndSplitter(Dipole * dip) {
  Proto_Particle * help(p_split); p_split = p_spect; p_spect = help;
  if (dip->IsSwitched()) dip->SetSwitched(false);
                    else dip->SetSwitched(true); 
  if (p_split->m_kt2max<=0.) 
    p_split->m_kt2max = p_spect->m_kt2max = ATOOLS::Max(p_split->m_kt2max,p_spect->m_kt2max);

  return p_split->m_kt2max;
}


bool Splitting_Tools::PrepareKinematics(Dipole * dip,const bool & first,const bool & enforce) {
  if (dip->MassBar2()<m_mmin_2) {
    if (dip->Triplet()->m_info=='B' || dip->AntiTriplet()->m_info=='B') {
      if (msg->LevelIsTracking()) {
	msg_Out()<<"Warning in "<<METHOD<<":"<<std::endl
		 <<"   Can not decay dipole for massmin^2 = "<<m_mmin_2<<"."<<std::endl;
	(*dip).Output();
      }
    }
    return false;
  }

  m_mom1           = p_spect->m_mom;
  m_mom3           = p_split->m_mom;
  m_mom0           = m_mom1+m_mom3;
  m_Q2             = m_mom0.Abs2();
  m_Q              = sqrt(m_Q2);

  m_m1             = p_spect->m_flav.HadMass();
  m_m2             = Flavour(kf_gluon).HadMass();
  m_m3   = m_m23   = p_split->m_flav.HadMass();
  m_m1_2           = sqr(m_m1);
  m_m2_2           = sqr(m_m2);
  m_m3_2 = m_m23_2 = sqr(m_m23);

  if (m_Q-m_m1-2.*sqrt(m_mmin_2)<0.) return false;
  m_Qt2            = m_Q2 - m_m1_2 - m_m2_2 - m_m3_2;
  if (m_Qt2<0.) {
    if (dip->Triplet()->m_info=='B' || dip->AntiTriplet()->m_info=='B') {
      msg_Tracking()<<"___ "<<METHOD<<" yields Qt2 = "<<m_Qt2<<"."<<std::endl;
    }  
    return false;
  }
  m_kt2_max = KT2Max(first);
  m_kt2_min = KT2Min(m_kt2_max,p_split->m_flav.IsGluon());


  if (dip->Triplet()->m_info=='B' || dip->AntiTriplet()->m_info=='B') {
    msg_Tracking()<<"___ "<<METHOD<<": set up dipole for splitting, "
		  <<"splitter = "<<p_split->m_flav<<"."<<std::endl;
  }
  m_cms  = Poincare(m_mom0);
  m_cms.Boost(m_mom0);
  m_cms.Boost(m_mom1);
  m_cms.Boost(m_mom3);

  m_zrot = Poincare(m_mom1,Vec4D::ZVEC);
  m_zrot.Rotate(m_mom1);
  m_zrot.Rotate(m_mom3);

  return true;
}




bool Splitting_Tools::DetermineSplitting(Dipole * dip1,const bool & vetodiquark) {
  if (dip1->Triplet()->m_info=='B' || dip1->AntiTriplet()->m_info=='B') {
    msg_Tracking()<<"___ "<<METHOD<<"(in): pt^2(max) = "<<m_kt2_max<<"."<<std::endl;
  }
  int  trials(0);
  m_z = -1.;
  m_flav = Flavour(kf_none);
  do {
    if (trials>1000 || 
	!ProduceKinematics(m_pt2,m_z,m_phi,m_flav,vetodiquark)) {
      if (dip1->Triplet()->m_info=='B' || dip1->AntiTriplet()->m_info=='B') {
	msg_Tracking()<<"___ "<<METHOD<<"(out): no splitting determined, trials = "<<trials<<"."<<std::endl;
      } 
      return false;
    }
    trials++;
  } while (!ConstructKinematics(dip1));

  msg_Tracking()<<"___ "<<METHOD<<"(out): splitting succeeded."<<std::endl;
  return true;
}


bool Splitting_Tools::EnforcedSplitting(Dipole * dip) {
  m_kt2 = m_mmin_2;
  m_z   = 0.5;
  SelectFlavour(m_kt2,true);
  
  if (!KinCheck(m_kt2,m_z) || !ConstructKinematics(dip)) {
    if (dip->Triplet()->m_info=='B' || dip->AntiTriplet()->m_info=='B') {
      msg_Tracking()<<METHOD<<" did not work out."<<std::endl;
    }
    return false;
  }
  return true;
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
  //msg_Out()<<"Out "<<METHOD<<":"<<std::endl<<(*p_split)<<(*p_spect)<<(*p_out1)<<(*p_out2)<<std::endl;
}





void Splitting_Tools::AftermathOfSplitting(Dipole * dip1) {
  if (p_split->m_flav.IsGluon()) {
    p_spect->m_mom    = m_mom1;
    p_out1 = new Proto_Particle(m_flav.Bar(),m_mom2,'l');
    p_out2 = new Proto_Particle(m_flav,m_mom3,'l');
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
    msg_Tracking()<<"___ "<<METHOD<<" yields new particles : "<<std::endl
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




bool Splitting_Tools::ProduceKinematics(double & kt2_test,double & z_test,double & phi_test,
					ATOOLS::Flavour & flav,
					const bool & vetodiquark) {
  double disc(4.*m_kt2_min/m_Qt2), discp(disc); 
  if (disc<0. || disc>1.) {
    msg_Tracking()<<"___ "<<METHOD<<": Cannot split object with reduced mass = "
		  <<sqrt(m_Qt2)<<"."<<std::endl;
    return false;
  }
  double deltaz(sqrt(1.-disc)),zmin((1.-deltaz)/2.), zmax((1.+deltaz)/2.); 
  double deltazp(deltaz), zminp(zmin), zmaxp(zmax);
  kt2_test = -1.;
  int  trials(0);
  msg_Tracking()<<"___ "<<METHOD<<" : start inner loop with pt2 in ["<<m_kt2_min<<", "<<m_kt2_max<<"] "
  		<<"and z = "<<z_test<<" for "<<m_flav<<". "
  		<<"Glusplit = "<<p_split->m_flav.IsGluon()<<"."<<std::endl;
  while ((kt2_test<0. || z_test<0.) && trials<1000) {
    trials++;
    kt2_test = p_as->SelectPT(m_kt2_max,m_kt2_min);
    if (p_split->m_flav.IsGluon()) {
      if (m_masstreatment==3) {
	SelectFlavour(kt2_test,vetodiquark);
	disc = m_m2_2/m_Qt2;
	if (disc<0. || disc>1.) continue;
      }
      else {
	SelectFlavour(m_Qt2/4.,vetodiquark);
      }
    }
    z_test = p_kernels->SelectZ(zmin,zmax,p_split->m_flav.IsGluon(),p_split->m_info=='L');

    if (true) {
      discp   = 4.*kt2_test/m_Qt2;
      deltazp = sqrt(1.-disc);
      zminp   = (1.-deltaz)/2.;
      zmaxp   = (1.+deltaz)/2.;
      if (z_test<zminp || z_test>zmaxp) {
	kt2_test = -1.;
	continue;
      }
    }

    if (!KinCheck(kt2_test,z_test)) {
      kt2_test = -1;
      continue;
    }
    if (p_kernels->Weight(kt2_test,z_test,p_split->m_flav.IsGluon(),
			  p_split->m_info=='L' && p_spect->m_info!='L')<ran.Get()) {
      kt2_test = -1;
      continue;
    }
  }
  msg_Tracking()<<"___ "<<METHOD<<" : exit inner loop with kt2 = "<<kt2_test<<" "
		<<"and z = "<<z_test<<" for "<<m_flav<<"."<<std::endl;

  //if (p_split->m_flav.IsGluon()) {
  // msg_Out()<<"Out of "<<METHOD<<"(veto = "<<vetodiquark<<", mass = "<<m_masstreatment<<"): "
  //	     <<"Flavour = "<<m_flav<<" for pt^2 = "<<kt2<<",    Qt^2 = "<<m_Qt2<<"."<<std::endl;
  //}
  //std::cout<<"Out of "<<METHOD<<" with "<<trials<<" trials."<<std::endl;
  phi_test = 2.*M_PI*ran.Get();
  return (trials<1000);
}

double Splitting_Tools::KT2Max(const bool & first) {
  double kt2max = Min((p_split->m_flav.IsGluon()?m_Qt2:m_Qt2/4.),m_pt2max);
  if (m_ptorder==PTOrder::flat) return kt2max;

  if (!first &&
      ((m_ptorder==PTOrder::gluon_split && p_split->m_flav==Flavour(kf_gluon)) ||
       (m_ptorder==PTOrder::gluon_emit && p_split->m_flav!=Flavour(kf_gluon)) ||
       (m_ptorder==PTOrder::total))) {
    kt2max    = p_split->m_kt2max;
    if (p_spect->m_info=='L' && !p_spect->m_flav.IsGluon() && m_m1_2>1.e-6) 
      kt2max *= m_pt2min/m_m1_2;
  }
  else {
    if (p_spect->m_info=='B') {
      if (p_split->m_info=='B') 
	kt2max = m_pt2max_factor*sqrt(4.*p_spect->m_mom.PPerp2()*p_split->m_mom.PPerp2());
      else 
	kt2max = m_pt2max_factor*sqrt(4.*p_split->m_mom.PPerp2()*m_pt2min);
    }
    else {
      kt2max = Min(m_Qt2/4.,p_spect->m_mom.PPerp2(p_split->m_mom));
      if (p_spect->m_info=='L' && !p_spect->m_flav.IsGluon() && m_m1_2>1.e-6) 
	kt2max *= m_pt2min/Max(m_pt2min,m_m1_2);
      //kt2max = Max(m_pt2max,sqrt(m_pt2max*kt2max));
    }
  }
    
  msg_Tracking()<<METHOD<<"(first = "<<first<<") yields kt2_max = "<<kt2max<<"("<<sqr(kt2max)/m_pt2max<<")"
		<<" for:"<<std::endl<<(*p_split)<<(*p_spect);
  return kt2max;
}
    



double Splitting_Tools::KT2Min(const double & kt2max,const bool & glusplit) {
  double kt2min(m_pt2min);
  if (m_m3_2>0.) kt2min *= m_pt2min/m_m3_2;
  if (m_m1_2>0.) kt2min *= m_pt2min/m_m1_2;
  if ((kt2min>kt2max) && (kt2min<m_Qt2/4.)) kt2min = 0.;
  if (glusplit) kt2min = Max(kt2min,m_mmin_2);
  return kt2min;
}


bool Splitting_Tools::KinCheck(const double & kt2_test,const double & z_test) {
  //double mu1(m_m1_2/m_Q2),mu2(m_m2_2/m_Q2),mu3(m_m3_2/m_Q2);
  if (m_masstreatment==3 && p_split->m_flav.IsGluon()) {
    m_s23 = kt2_test;
    m_kt2 = m_Qt2*m_s23-(1.-z_test)/z_test*m_m2_2-z_test/(1.-z_test)*m_m3_2;
    if (m_kt2<0.0) return false;
    m_y   = (m_s23-m_m2_2-m_m3_2)/m_Qt2;
  }
  else { 
    m_kt2 = kt2_test;
    m_y   = (kt2_test - sqr(1.-z_test)*m_m2_2 - sqr(z_test)*m_m3_2)/(z_test*(1.-z_test)*m_Qt2);
    m_s23 = m_y*m_Qt2+m_m2_2+m_m3_2;
    if (m_s23<0.0) return false;
    if (m_masstreatment==3 && p_split->m_flav.IsGluon()) {
      if (m_kt2/m_s23*(*p_as)(m_s23)/(*p_as)(m_kt2)<ran.Get()) return false;
    }
  }
  double ym = 2.*m_m2*m_m3/m_Qt2;
  double yp = (m_Q2-2.*m_m1_2/m_Q2*(m_Q2-m_m3_2))/m_Qt2;
  if (m_y<ym || m_y>yp) return false;

  double viji = sqrt(sqr(m_Qt2*m_y)-4.*m_m2_2*m_m3_2)/(m_Qt2*m_y+2.*m_m3_2);
  double vijk = sqrt(sqr(2.*m_m1_2+m_Qt2*(1.-m_y))-4.*m_m1_2)/(m_Qt2*(1.-m_y));
  double frac = (2.*m_m3_2+m_Qt2*m_y)/(2.*(m_m3_2+m_m2_2+m_Qt2*m_y));
  double zm   = frac*(1.-viji*vijk);
  double zp   = frac*(1.+viji*vijk);
  if (z_test<zm || z_test>zp) return false;

  p_kernels->SetScales(m_m1_2,m_m2_2,m_m3_2,m_Q2,m_y,
		       viji,vijk,zm,zp,m_y*m_Qt2/2.);

  return true;
}


void Splitting_Tools::SelectFlavour(const double & maxmass,const bool & vetodiquark) {
  double sumwt(0.);
  m_flav = Flavour(kf_none);
  //std::cout<<"-------------------------------------------------"<<std::endl;
  while (m_flav==Flavour(kf_none)) {
    for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
      if (vetodiquark && fdit->first.IsDiQuark()) continue;
      if (fdit->second->popweight>0. && fdit->second->massmin<maxmass) {
	sumwt += fdit->second->popweight;
	//std::cout<<"  add "<<fdit->first<<" with pop = "<<fdit->second->popweight<<", "
	//	 <<"mass = "<<fdit->second->massmin<<" (max = "<<maxmass<<")  "
	//	 <<"--> sum = "<<sumwt<<"."<<std::endl;
      }
    }
    if (sumwt==0) {
      m_flav=(ran.Get()<0.5)?Flavour(kf_d):Flavour(kf_u);
      if (m_flav==Flavour(kf_d)) m_d++;
      else if (m_flav==Flavour(kf_u)) m_u++;
      m_m2    = m_m3    = m_flav.HadMass();
      m_m2_2  = m_m3_2  = sqr(m_m2); 
      m_Qt2   = m_Q2 - m_m1_2 - m_m2_2 - m_m3_2;
      return;
    }
    sumwt *= ran.Get();
    for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
      if (vetodiquark && fdit->first.IsDiQuark()) continue;
      if (fdit->second->popweight>0. && fdit->second->massmin<maxmass) {
	//std::cout<<"  sub "<<fdit->first<<" with pop = "<<fdit->second->popweight<<", "
	//	 <<"mass = "<<fdit->second->massmin<<"(max = "<<maxmass<<")  "
	//	 <<"<-- sum = "<<sumwt<<"."<<std::endl;
	sumwt -= fdit->second->popweight;
      }
      if (sumwt<0.) {
	m_flav = fdit->first;
	if (m_flav.IsDiQuark()) m_flav = m_flav.Bar(); 
	break;
      }
    }
    sumwt = 0.;
  }
  if (m_flav==Flavour(kf_d)) m_d++;
  else if (m_flav==Flavour(kf_u)) m_u++;
  else if (m_flav==Flavour(kf_s)) m_s++;
  m_tot++;

  //std::cout<<"-------------------------------------------------"<<std::endl;
  m_m2    = m_m3    = m_flav.HadMass();
  m_m2_2  = m_m3_2  = sqr(m_m2); 
  m_Qt2   = m_Q2 - m_m1_2 - m_m2_2 - m_m3_2;
}


bool Splitting_Tools::ConstructKinematics(const bool glusplit) {
  msg_Tracking()<<"___ "<<METHOD<<"(in) for kt^2 = "<<m_kt2<<", z = "<<m_z<<"."<<std::endl;
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
      msg_Tracking()<<"      ----> (fail0) for kt^2 = "<<m_kt2<<", z = "<<m_z<<std::endl;
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
      msg_Tracking()<<"      ----> (fail1) for kt^2 = "<<m_kt2<<", z = "<<m_z<<std::endl;
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
      msg_Tracking()<<"      ----> (fail2) for kt^2 = "<<m_kt2<<", z = "<<m_z<<std::endl;
      return false;
    }
    m_kt = sqrt(m_kt2);
    //    q1 = zt*l + (mi2+ktt*ktt)/(gam*zt)*n +
    //  + ktt*cos(phi)*Vec4D(0.0,1.0,0.0,0.0) 
    //  + ktt*sin(phi)*Vec4D(0.0,0.0,1.0,0.0); 
    q3 = m_z*l + (m_m3_2+m_kt2)/(gam*m_z)*n + m_kt*nperp;
  }
  q2 = m_mom0-q3-q1;

  if (m_m1_2==0.) {
    msg_Tracking()<<"___ "<<METHOD<<" summary: "<<std::endl
		  <<"     Masses: m_m2 = m_m3 = "<<m_m2<<" = "<<m_m3<<" --> s23 = "<<m_m23_2<<"."<<std::endl;
  }

  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<"      ----> (fail3) for kt^2 = "<<m_kt2<<", z = "<<m_z<<std::endl;
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

  msg_Tracking()<<"___ "<<METHOD<<"(out): succeeded to build kinematics "
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


/*
double Gluon_Decayer::PT2Max(Dipole * dip) const {
  double pt2max(Min(Min(dip->Triplet()->m_mom.PPerp2(dip->AntiTriplet()->m_mom),
			dip->AntiTriplet()->m_mom.PPerp2(dip->Triplet()->m_mom)),
		    dabs((dip->Triplet()->m_mom-dip->AntiTriplet()->m_mom).Abs2())));
  if (IsZero(pt2max)) pt2max = dip->MassBar2();
  if (dip->Triplet()->m_info=='B' && dip->AntiTriplet()->m_info!='B') {
    if (dip->AntiTriplet()->m_mom.PPerp2()<1.e-6) pt2max = Min(pt2max,m_pt02);
    else pt2max = Min(pt2max,dip->AntiTriplet()->m_mom.PPerp2());
  }
  else if (dip->AntiTriplet()->m_info=='B' && dip->Triplet()->m_info!='B') {
    if (dip->Triplet()->m_mom.PPerp2()<1.e-6) pt2max = Min(pt2max,m_pt02);
    else pt2max = Min(pt2max,dip->Triplet()->m_mom.PPerp2());
  }
  else if (dip->AntiTriplet()->m_info=='B' && dip->Triplet()->m_info=='B') {
    pt2max = Min(pt2max,sqrt(dip->Triplet()->m_mom.PPerp2()*dip->AntiTriplet()->m_mom.PPerp2()));
  }
  if (pt2max>m_pt2max) {
    msg_Tracking()<<"Unusually large pt2max in "<<METHOD<<" for : "<<std::endl;
    if (msg_LevelIsTracking()) dip->Output();
    msg_Tracking()<<"   pt^2 = "<<dip->Triplet()->m_mom.PPerp2(dip->AntiTriplet()->m_mom)
		  <<" vs. pt^2 = "<<dip->AntiTriplet()->m_mom.PPerp2(dip->Triplet()->m_mom)
		  <<" vs. m^2 = "<<(dip->Triplet()->m_mom-dip->AntiTriplet()->m_mom).Abs2()
		  <<" --> "<<(m_pt2max_factor*pt2max)<<"."<<std::endl;
    pt2max = m_pt2max;
  }
  return m_pt2max_factor * pt2max;
}
*/


