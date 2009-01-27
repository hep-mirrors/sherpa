#include "Dipole_Splitter.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;

size_t Dipole::s_cnt=0;
const Vec3D Dipole_Splitter::s_ex(Vec3D(1.,0.,0.));
const Vec3D Dipole_Splitter::s_ey(Vec3D(0.,1.,0.));
const Vec3D Dipole_Splitter::s_ez(Vec3D(0.,0.,1.));


Dipole_Splitter::Dipole_Splitter(Strong_Coupling * as,const double ptmax) :
  m_massreweighting(false), 
  m_leading(hadpars.Get(std::string("leading_particles"))<2), m_flat(false), m_pole(true),
  p_as(as), p_constituents(hadpars.GetConstituents()), p_options(NULL),
  p_spect(0), p_split(0), p_out1(0), p_out2(0), 
  m_mmin_2(sqr(p_constituents->MinMass()))
{ 
  m_histograms[std::string("PT_Gluon_Splitting")]      = new Histogram(0,0.,5.,50);
  m_histograms[std::string("PT_Gluon_Emission")]       = new Histogram(0,0.,5.,50);
  m_histograms[std::string("Z_Gluon_Splitting")]       = new Histogram(0,0.,1.,50);
  m_histograms[std::string("Z_Gluon_Emission")]        = new Histogram(0,0.,1.,50);
  m_histograms[std::string("ZWeight_Gluon_Splitting")] = new Histogram(0,0.,1.,50);
  m_histograms[std::string("ZWeight_Gluon_Emission")]  = new Histogram(0,0.,1.,50);
}


Dipole_Splitter::~Dipole_Splitter() {
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

bool Dipole_Splitter::SplitCluster(SP(Cluster) cluster,const double pt2max,
				   const bool pole) {
  //std::cout<<METHOD<<"("<<pt2max<<", "<<pole<<")."<<std::endl;
  if (m_leading) m_pole = pole;
  else m_pole = true;

  SP(Dipole) dip1 = new Dipole(new Proto_Particle((*cluster->GetTrip())),
			       new Proto_Particle((*cluster->GetAnti())));
#ifdef memchecker
  msg_Out()<<"### Two new Proto_Particles ("
  	   <<dip1->Triplet()<<"/"<<dip1->AntiTriplet()<<") from "<<METHOD<<"."<<std::endl;
#endif

  p_dip = dip1;

  SP(Dipole) dip2 = EmitGluon(pt2max);
  if (dip2==NULL) {
    msg_Tracking()<<"ERROR in "<<METHOD<<":"<<std::endl
		  <<"   Could not emit gluon, hence no dipole splitting."<<std::endl
		  <<"   Return false, may lead to new event."<<std::endl;
#ifdef memchecker
    msg_Out()<<"### Delete Proto_Particle ("<<dip1->Triplet()
    	     <<"/"<<dip1->Triplet()->m_flav<<") in "<<METHOD<<"."<<std::endl;
#endif
#ifdef memchecker
    msg_Out()<<"### Delete Proto_Particle ("<<dip1->AntiTriplet()
    	     <<"/"<<dip1->AntiTriplet()->m_flav<<") in "<<METHOD<<"."<<std::endl;
#endif
    return false;
  }
  dip1->SetAntiTriplet(dip2->Triplet());
  dip1->Update();
  dip2->Update();

  SP(Dipole) first(dip1), second(dip2);
  if (ran.Get()<(dip1->Mass2()-dip1->Triplet()->m_mom.Abs2())/
      (dip1->Mass2()-dip1->Triplet()->m_mom.Abs2()+
       dip2->Mass2()-dip2->AntiTriplet()->m_mom.Abs2())) {
    first  = dip2;
    second = dip1;
  }
  if (!SplitDipole(first,pt2max,m_pole) && !SplitDipole(second,pt2max,m_pole)) {
    msg_Tracking()<<"Error in "<<METHOD<<" :"<<std::endl
		  <<"   Two unsplittable dipoles emerging from :"<<std::endl
		  <<(*cluster)<<std::endl;
    return false;
  }


  dip1->SetAntiTriplet(p_out1);
  dip2->SetTriplet(p_out2);
  
  dip1->Update();
  dip2->Update();

  SP(Cluster) left  = new Cluster(dip1->Triplet(),dip1->AntiTriplet());
  SP(Cluster) right = new Cluster(dip2->Triplet(),dip2->AntiTriplet());
#ifdef memchecker
  msg_Out()<<"@@@ Two new clusters "<<left<<"/"<<right<<" from "<<METHOD<<"."<<std::endl;
#endif

  cluster->SetLeft(left);
  cluster->SetRight(right);
  left->SetPrev(cluster);
  right->SetPrev(cluster);

  Vec4D check = cluster->Momentum()-cluster->GetLeft()->Momentum()-cluster->GetRight()->Momentum();
  if (!IsZero(check.Abs2()) || !IsZero(check[0]/1.e6)) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   Four-momentum not conserved: "<<check<<" ("<<check.Abs2()<<") for "<<std::endl
	       <<"   "<<cluster->Momentum()<<"  ---> "<<std::endl
	       <<"   "<<cluster->GetLeft()->Momentum()
	       <<" + "<<cluster->GetRight()->Momentum()<<"."<<std::endl;
    //abort();
    return false;
  }
  return true;
}

SP(Dipole) Dipole_Splitter::EmitGluon(const double pt2max) {
  SetSpectatorAndSplitter();
  if (!PrepareKinematics(pt2max) || !DetermineSplitting(false)) {
    //msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
    //	       <<"   Can not decay dipole for massmin^2 = "<<m_mmin_2<<"."<<std::endl;
    //(*p_dip).Output();
    return false;
  }
  p_spect->m_mom = m_mom1;
  p_out1 = new Proto_Particle(m_flav2,m_mom2,'l');
  p_split->m_mom = m_mom3;

  return new Dipole(p_out1,p_dip->AntiTriplet());
}

bool Dipole_Splitter::SplitDipole(SP(Dipole) dip,const double pt2max,
				  const bool pole) {
  //std::cout<<METHOD<<"("<<pt2max<<", "<<pole<<")."<<std::endl;
  m_pole = pole;
  p_dip  = dip;
  SetSpectatorAndSplitter();
  if (!PrepareKinematics(pt2max) || !DetermineSplitting(true)) {
    //msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
    //	       <<"   Can not decay dipole for massmin^2 = "<<m_mmin_2<<"."<<std::endl;
    //(*p_dip).Output();
    return false;
  }
  p_spect->m_mom = m_mom1;
  p_out1 = new Proto_Particle(m_flav2,m_mom2,'l');
  p_out2 = new Proto_Particle(m_flav3,m_mom3,'l');
  return true;
}

bool Dipole_Splitter::DetermineSplitting(const bool glusplit) {
  int trials(0);
  do {
    //if (!m_pole) 
    // std::cout<<METHOD<<": isotropic from glusplit "<<glusplit<<" & "<<m_pole
    //	       <<" from "<<m_leading<<"."<<std::endl;
    //else  
    // std::cout<<METHOD<<": non-isotropic from glusplit "<<glusplit<<" & "<<m_pole
    //	       <<" from "<<m_leading<<"."<<std::endl;
    m_kt2 = p_as->SelectPT(m_kt2_max,m_pole);
    if (glusplit) SelectFlavour();
    if (4.*m_kt2*(m_Qt2+m_m3_2)>m_Qt2*m_Qt2) continue;
    m_z   = SelectZ(glusplit); 
    m_phi = 2.*M_PI*ran.Get();
    trials++;
    if (trials>100) return false;
  } while (!ConstructKinematics(glusplit));
  return true;
}

void Dipole_Splitter::SelectFlavour() {
  double maxmass(sqrt((m_Q2-m_m1_2)/4.)),sumwt(0.);
  m_flav3 = p_options->begin()->first;
  m_flav2 = m_flav3.Bar();
  for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
    if (fdit->second->popweight>0. && fdit->second->massmin<maxmass) 
      sumwt += fdit->second->popweight;
  }
  sumwt *= ran.Get();
  for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
    if (fdit->second->popweight>0. && fdit->second->massmin<maxmass)
      sumwt -= fdit->second->popweight;
    if (sumwt<=0.) {
      m_flav3 = fdit->first;
      if (m_flav3.IsDiQuark()) m_flav3 = m_flav3.Bar(); 
      m_flav2 = m_flav3.Bar();
      break;
    }
  }
  m_m2    = m_m3    = p_constituents->Mass(m_flav3);
  m_m2_2  = m_m3_2  = sqr(m_m2); 
  m_Qt2   = m_Q2 - m_m1_2 - m_m2_2 - m_m3_2;
}

double Dipole_Splitter::SelectZ(const bool glusplit) {
  if (m_flat) return ran.Get();
  Histogram* histo;
  double weight(0.), mu2(0.), ztest(-1.);
  if (glusplit) {
    //std::cout<<"In g->QQ"<<std::endl;
    while (ztest<0.) {
      ztest  = ran.Get();
      mu2    = 2.*m_m2_2*ztest*(1.-ztest)/(m_kt2+m_m2_2);
      weight = (1.-2.*ztest*(1.-ztest))*(1.+mu2/(1.-2.*ztest*(1.-ztest)));
      histo  = (m_histograms.find(std::string("ZWeight_Gluon_Splitting")))->second;
      histo->Insert(weight);
      if (weight>ran.Get()) {
	histo  = (m_histograms.find(std::string("Z_Gluon_Splitting")))->second;
	histo->Insert(ztest);
	return ztest;
      }
      ztest = -1.;
    }
  }
  else {
    double deltaz = sqrt(1.-4.*m_kt2*(m_Qt2+m_m3_2)/(m_Qt2*m_Qt2));
    double zmin   = (1.-deltaz)/2., zmax = (1.+deltaz)/2.; 
    //std::cout<<"In q->qg : "<<zmin<<" ... "<<zmax<<"."<<std::endl;
    while (ztest<0.) {
      ztest  = 1.-(1.-zmin)*pow((1.-zmax)/(1.-zmin),ran.Get());
      mu2    = (ztest*(1.-ztest)*m_m3_2)/(m_kt2+sqr(1.-ztest)*m_m3_2);
      weight = ((1.+sqr(ztest))/2.) * (1.-2.*(1.-ztest)*mu2/(1.+sqr(ztest)));
      histo  = (m_histograms.find(std::string("ZWeight_Gluon_Emission")))->second;
      histo->Insert(weight);
      if (weight>ran.Get()) {
	histo  = (m_histograms.find(std::string("Z_Gluon_Emission")))->second;
	histo->Insert(ztest);
	return ztest;
      }
      ztest = -1;
    }
  }
  return ran.Get();
}

bool Dipole_Splitter::PrepareKinematics(const double pt2max) {
  m_mom1    = p_spect->m_mom;
  m_mom2    = p_split->m_mom;
  m_mom0    = m_mom1+m_mom2;

  m_Q2      = m_mom0.Abs2();
  m_Q       = sqrt(m_Q2);
  m_cms     = Poincare(m_mom0);
  m_cms.Boost(m_mom0);
  m_cms.Boost(m_mom1);
  m_cms.Boost(m_mom2);

  m_zrot   = Poincare(m_mom1,Vec4D::ZVEC);
  m_zrot.Rotate(m_mom1);
  m_zrot.Rotate(m_mom2);
  
  m_flav2          = Flavour(kf_gluon);
  m_flav3          = p_split->m_flav;

  m_m1             = p_constituents->Mass(p_spect->m_flav);
  m_m3   = m_m23   = p_constituents->Mass(p_split->m_flav);
  m_m1_2           = sqr(m_m1);
  m_m3_2 = m_m23_2 = sqr(m_m23);
  m_m2             = p_constituents->Mass(m_flav2);
  m_m2_2           = sqr(m_m2);


  if (m_Q-m_m1-2.*sqrt(m_mmin_2)<0.) return false;
  
  if (pt2max>0.) m_kt2_max = ATOOLS::Min((m_Q2-m_m1_2)/4.,pt2max);
            else m_kt2_max = (m_Q2-m_m1_2)/4.;

  if (m_kt2_max<0.) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   No physical splitting possible for pt2max = "
	       <<pt2max<<" m_kt2_max = "<<m_kt2_max<<std::endl;
    return false;
  }
  return true;
}

bool Dipole_Splitter::ConstructKinematics(const bool glusplit) {
  Vec4D q1(0.,0.,0.,0.),q2(0.,0.,0.,0.),q3(0.,0.,0.,0.);
  Vec4D nperp = Vec4D(0.,cos(m_phi)*s_ex + sin(m_phi)*s_ey);
  if (m_m1_2==0.) {
    double kt   = sqrt(m_kt2);
    double fac1 = m_Q2 - ((1.-m_z)*m_m2_2+m_z*m_m3_2+m_kt2)/(m_z*(1.-m_z));
    if (fac1<0.) return false;
    double div  = m_Q2-m_m23_2;
    fac1 /= div;
    
    q1 =                                                             fac1 * m_mom1;
    q2 =      m_z * m_mom2 +    (m_kt2+m_m2_2-sqr(m_z)*m_m23_2)/(m_z*div) * m_mom1 + kt * nperp;
  }
  else {
    double y   =
      (m_kt2/(m_z*(1.-m_z))+(1.-m_z)/m_z*m_m2_2+m_z/(1.-m_z)*m_m3_2)/
      (m_Q2-m_m1_2-m_m2_2-m_m3_2);
    double s23 = y*(m_Q2-m_m1_2)+(1.-y)*(m_m2_2+m_m3_2);
    double po  = sqr(m_Q2-m_m23_2-m_m1_2)-4.*m_m23_2*m_m1_2; 
    double pn  = sqr(m_Q2-s23-m_m1_2)-4.0*s23*m_m1_2;
    if (po<0.0 || pn<0.0) return false;

    q1 = 
      (m_Q2+m_m1_2-s23)/(2.*m_Q2)*m_mom0+
      sqrt(pn/po) * (m_mom1 - (m_Q2+m_m1_2-m_m23_2)/(2.*m_Q2) * m_mom0);

    Vec4D  q23 = m_mom0-q1;
    double gam = q23*q1+sqrt(sqr(q23*q1)-s23*m_m1_2);
    double a23 = s23/gam, a1=m_m1_2/gam, beta=1.0/(1.0-a23*a1);
    Vec4D  l   = beta*(q23-a23*q1);
    Vec4D  n   = beta*(q1-a1*q23);
    double zt  =
      (m_z/a1-y/(1.-y)-2.*m_m2_2/gam/(1.+a23*a1))/
      ((1./a1-(m_m2_2+m_m3_2)/gam)/(1.0+a23*a1)-y/(1.0-y));
    double lt2 = gam*y/(1.-y)*(1.+a23*a1)*zt*(1.-zt)-sqr(1.-zt)*m_m2_2-zt*zt*m_m3_2;
    if (lt2<0.0) return false;

    double lt  = sqrt(lt2);
    q2 = zt*l + (m_m2_2+lt2)/(gam*zt)*n + lt*nperp;    
    if (q1[0]<0. || q2[0]<0. || q3[0]<0.) return false;
  }
  q3 = m_mom0-q2-q1;

  if (m_massreweighting && glusplit) {
    // Reweight with (pt2+pt02)/s23
    m_mass2 = (q2+q3).Abs2();    
    if (p_as->GetPTArgument()<m_mass2*ran.Get()) return false;
  }

  m_zrot.RotateBack(q1);
  m_zrot.RotateBack(q2);
  m_zrot.RotateBack(q3);
  m_cms.BoostBack(q1);
  m_cms.BoostBack(q2);
  m_cms.BoostBack(q3);

  m_mom1 = q1;
  m_mom2 = q2;
  m_mom3 = q3;

  m_cms.BoostBack(m_mom0);

  Vec4D check = (m_mom0-m_mom1-m_mom2-m_mom3);
  if (!IsZero(check.Abs2())) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   4 mom-check failed for kinematics: "
	       <<check<<" ("<<check.Abs2()<<")."<<std::endl;
    return false;
  }
  return true;
}

void Dipole_Splitter::SetSpectatorAndSplitter() {
  if (p_dip->Triplet()->m_flav.IsGluon() && 
      (p_dip->AntiTriplet()->m_flav.IsQuark() ||
       p_dip->AntiTriplet()->m_flav.IsDiQuark())) {
    p_spect = p_dip->AntiTriplet(); 
    p_split = p_dip->Triplet();
    p_dip->SetSwitched(true);
  }
  else if (p_dip->AntiTriplet()->m_flav.IsGluon() && 
	   (p_dip->Triplet()->m_flav.IsQuark() ||
	    p_dip->Triplet()->m_flav.IsDiQuark())) {
    p_split = p_dip->AntiTriplet(); 
    p_spect = p_dip->Triplet();
    p_dip->SetSwitched(false);
  }
  else if ((p_dip->Triplet()->m_flav.IsGluon() && 
	    p_dip->AntiTriplet()->m_flav.IsGluon()) ||
	   (!p_dip->Triplet()->m_flav.IsGluon() && 
	    !p_dip->AntiTriplet()->m_flav.IsGluon())) {
    if (ran.Get()>0.5) {
      p_spect = p_dip->AntiTriplet(); 
      p_split = p_dip->Triplet();
      p_dip->SetSwitched(true);
    }
    else {
      p_split = p_dip->AntiTriplet(); 
      p_spect = p_dip->Triplet();
      p_dip->SetSwitched(false);
    }
  }
}


