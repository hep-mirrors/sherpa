#include "Dipole_Splitter.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;




////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


Dipole_Splitter::Dipole_Splitter(Strong_Coupling * as) :
  m_recoils(recoils::energy_dep), p_as(as), 
  p_spect(0), p_split(0), p_out1(0), p_out2(0),
  m_pt2min(0.01), 
  p_constituents(hadpars.GetConstituents())
{ }

bool Dipole_Splitter::Split(Dipole * dip,const double pt2max) {
  PrepareDipole(dip,pt2max);
  msg_Out()<<METHOD<<" for : "<<dip->massbar2<<"/"<<dip->mass2<<std::endl
	   <<"   "<<p_spect->m_flav<<" "<<p_spect->m_mom<<" with : "
	   <<hadpars.GetConstituents()->Mass(p_spect->m_flav)<<"/"<<m_kin.W<<std::endl
	   <<"   "<<p_split->m_flav<<" "<<p_split->m_mom<<std::endl;    

  int counter(0);
  do { 
    m_kin.xt2tmp = 0.;
    DetermineSplittingMode();
    counter++;
  } while(m_kin.xt2tmp==0. && counter<10);
  if (counter>=10 && m_kin.xt2tmp==0.) {
    
  }

  msg_Out()<<"METHOD - testsplit : "
	   <<m_kin.x1/2.<<" vs. "<<m_kin.sm1<<", "
	   <<m_kin.x2/2.<<" vs. "<<m_kin.sm2<<", "
	   <<m_kin.x3/2.<<" vs. "<<m_kin.sm3<<"."<<std::endl;
  if (FixKinematics(dip)) {
    msg_Out()<<METHOD<<" yields success."<<std::endl
	     <<(*p_out1)<<(*p_out2)<<std::endl;
    return true;
  }
  msg_Out()<<"ERROR in "<<METHOD<<" FixKinematics did not work out!"<<std::endl;
  return false;
}

void Dipole_Splitter::PrepareDipole(Dipole * dip,const double pt2max) {
  if (dip->triplet->m_flav.IsGluon() && dip->antitriplet->m_flav.IsQuark()) {
    p_spect = dip->antitriplet; 
    p_split = dip->triplet;
    dip->switched = true;
  }
  else {
    p_split = dip->antitriplet; 
    p_spect = dip->triplet;
    dip->switched = false;
  }
  m_kin.s      = dip->mass2;
  m_kin.W      = sqrt(dip->mass2);
  m_kin.sm1    = p_constituents->Mass(p_spect->m_flav)/m_kin.W;
  m_kin.sm1_2  = sqr(m_kin.sm1);

  m_kin.yint   = 2.*m_kin.W;
  m_evolval    = 16.*M_PI/(m_kin.yint*p_as->MaxValue());

  if (pt2max>0.) m_kin.xt2cut = pt2max/m_kin.s;
  else m_kin.xt2max = 1.;

  m_kin.flav   = Flavour(kf::none);
}

void Dipole_Splitter::DetermineSplittingMode() {
  bool   nosplit(false);
  double xt2test,ytest;
  for (FDIter fdit=p_options->begin();fdit!=p_options->end();fdit++) {
    msg_Out()<<"   In "<<METHOD<<" for new flav = "
    	     <<fdit->first<<" ("<<fdit->second->popweight<<") with "
    	     <<fdit->second->massmin<<std::endl;
    if (fdit->second->popweight<ran.Get()) continue;
    if (p_constituents->Mass(p_spect->m_flav)+
	2.*fdit->second->massmin>m_kin.W)  continue;
    if (!BuildKinematicBounds(fdit))       continue;

    xt2test = 1.;
    nosplit = false;
    do { if (!SelectPT_Y(xt2test,ytest)) { nosplit=true; break; } 
    } while (Veto(xt2test,ytest));
    if (nosplit) {
      //msg_Out()<<"In "<<METHOD<<" found no splitting with pt larger than so far."<<std::endl;
      continue;
    }
    if (xt2test>m_kin.xt2tmp) {
      m_kin.xt2tmp = xt2test;
      m_kin.ytmp   = ytest;
      m_kin.flav   = fdit->first;
      m_kin.x1     = m_kin.smom1-sqrt(m_kin.xt2tmp)*exp(m_kin.ytmp);
      m_kin.x3     = m_kin.smom3-sqrt(m_kin.xt2tmp)*exp(-m_kin.ytmp);
      m_kin.x2     = 2.-m_kin.x1-m_kin.x3;
      msg_Out()<<"   Found an allowed splitting mode, flav = "<<m_kin.flav
	       <<" -> "<<m_kin.xt2tmp<<"."<<std::endl;
    }
  }
  msg_Out()<<"***Out "<<METHOD<<": flav = "<<m_kin.flav<<" -> "<<m_kin.xt2tmp<<"."<<std::endl;
}

bool Dipole_Splitter::BuildKinematicBounds(FDIter & fdit) {
  //msg_Out()<<"In "<<METHOD<<":"<<std::endl;
  m_kin.sm2   = m_kin.sm3   = p_constituents->Mass(fdit->first)/m_kin.W;
  m_kin.sm2_2 = m_kin.sm3_2 = sqr(m_kin.sm2);
  double sum_scaled_masses(m_kin.sm1+m_kin.sm2+m_kin.sm3);
  if (sum_scaled_masses>1.) return false;
  //msg_Out()<<"   Sum of scaled masses okay  : "<<sum_scaled_masses<<std::endl;
  double cms_mom(1.+sqr(m_kin.sm1_2-m_kin.sm3_2)-2.*(m_kin.sm1_2+m_kin.sm3_2));
  if (cms_mom<0.)           return false;
  //msg_Out()<<"   Three momentum in cms okay : "<<cms_mom<<std::endl;
  m_kin.smom1 = 1+m_kin.sm1_2-sqr(m_kin.sm2+m_kin.sm3);
  m_kin.smom3 = 1+m_kin.sm3_2-sqr(m_kin.sm1+m_kin.sm2);
  m_kin.xt1   = m_kin.smom1-2.*m_kin.sm1;
  m_kin.xt3   = m_kin.smom3-2.*m_kin.sm3;
  if (m_kin.xt1<0. || m_kin.xt3<0.)     return false;
  //msg_Out()<<"   XT(1,3) okay               : "<<m_kin.xt1<<", "<<m_kin.xt3<<std::endl;
  double xtkinmax(sqrt(0.25+m_kin.sm2_2)-1+0.5*(m_kin.smom1+m_kin.smom3));
  if (xtkinmax<0.)          return false;
  //msg_Out()<<"   XT_kinmax okay             : "<<xtkinmax<<std::endl;
  m_kin.xt2kinmax = sqr(xtkinmax);
  m_kin.xt2max    = Min(m_kin.xt2kinmax,m_kin.xt1*m_kin.xt3);
  if (m_kin.xt2max<m_kin.xt2tmp) return false;
  return true;
}

bool Dipole_Splitter::SelectPT_Y(double & xt2,double & y) {
  double random1(ran.Get()),random2(ran.Get());
  if (xt2>m_kin.xt2kinmax) xt2 = m_kin.xt2kinmax;
  //if (integ*log(random1)<(sqr(log(xt2))-sqr(log(pt2min/mass2)))) return false;
  //xt2  = exp(-sqrt(sqr(log(xt2))-log(random1)*integ));
  if (log(random1)*m_evolval<log(m_pt2min/(m_kin.s*xt2))) return false;  
  do { xt2 *= pow(random1,m_evolval); } while (xt2>m_kin.xt2cut);
  if (xt2<m_kin.xt2tmp) return false;
  double zmax = sqrt(m_kin.xt2kinmax/xt2)+sqrt(m_kin.xt2kinmax/xt2-1.);
  double xt   = sqrt(xt2);
  double zmin = Min(zmax,m_kin.xt1/xt);
  zmax        = Min(zmax,m_kin.xt3/xt);
  m_kin.ymin  = -log(zmin);
  m_kin.ymax  =  log(zmax);
  y    = -log(1./zmax+random2*(zmin-1./zmax));
  //msg_Out()<<"In "<<METHOD<<": "<<std::endl
  //	   <<"   xt2, y = "<<xt2<<", "<<y<<" from zmin,max = "<<zmin<<", "<<zmax
  //	   <<" xts = "<<m_kin.xt2kinmax<<", xt2 = "<<xt2<<","
  //	   <<" xt1,3 = "<<m_kin.xt1<<", "<<m_kin.xt3<<std::endl; 
  return true;
}

bool Dipole_Splitter::Veto(const double xt2,const double y) { 
  //msg_Out()<<"In "<<METHOD<<":"<<std::endl;
  if ((*p_as)(xt2*m_kin.s)/p_as->MaxValue()<ran.Get()) {
    //msg_Out()<<"   as-weight : "<<(*p_as)(xt2*m_kin.s)<<"/"<<p_as->MaxValue()<<"."<<std::endl;
    return true;
  }
  double xt(sqrt(xt2));
  double x1   = m_kin.smom1-xt*exp(y);
  double x3   = m_kin.smom3-xt*exp(-y);
  double x2   = 2.-x1-x3;
  if (x1/2.<m_kin.sm1 || x2/2.<m_kin.sm2 || x3/2.<m_kin.sm3) return true;
  double ss12 = 1.-x3+m_kin.sm3_2;
  double ss23 = 1.-x1+m_kin.sm1_2;
  double ss13 = 1.-x2+m_kin.sm2_2;
  if (ss12<sqr(m_kin.sm1+m_kin.sm2) ||
      ss13<sqr(m_kin.sm1+m_kin.sm3) ||
      ss23<sqr(m_kin.sm2+m_kin.sm3)) return true;

  double mecorr = sqr((1.-x3+sqr(m_kin.sm3)))+
    sqr((1.-x2+sqr(m_kin.sm2)))*xt*(exp(-m_kin.ymin)-exp(-m_kin.ymax))/(2.*m_kin.W);
  if (mecorr<ran.Get()) {
    //msg_Out()<<"   me-corr.  : "<<mecorr<<std::endl;
    return true;
  }
  msg_Out()<<METHOD<<" - energies/masses match : "
	   <<x1/2.<<" vs. "<<m_kin.sm1<<", "
	   <<x2/2.<<" vs. "<<m_kin.sm2<<", "
	   <<x3/2.<<" vs. "<<m_kin.sm3<<";"<<std::endl
	   <<"    and the s match : "
	   <<ss12<<" vs. "<<sqr(m_kin.sm1+m_kin.sm2)<<", "
	   <<ss13<<" vs. "<<sqr(m_kin.sm1+m_kin.sm3)<<", "
	   <<ss23<<" vs. "<<sqr(m_kin.sm2+m_kin.sm3)<<"."<<std::endl;
  
  return false; 
}

bool Dipole_Splitter::FixKinematics(Dipole * dip) {
  p_out1      = p_out2 = NULL; 
  Vec4D momsplit(p_split->m_mom), momspect(p_spect->m_mom), cms(momsplit+momspect);
  m_kin.sm1   = p_constituents->Mass(p_spect->m_flav)/m_kin.W;
  m_kin.sm2   = m_kin.sm3 = p_constituents->Mass(m_kin.flav)/m_kin.W;
  m_kin.sm1_2 = sqr(m_kin.sm1);
  m_kin.sm2_2 = m_kin.sm3_2 = sqr(m_kin.sm2);
 
  m_booster   = Poincare(cms);
  m_rotator   = Poincare(momsplit,Vec4D::ZVEC);
  
  double sE1  = m_kin.x1/2.;
  double sE2  = m_kin.x2/2.;
  double sE3  = m_kin.x3/2.;
  if (sE1<m_kin.sm1 || sE2<m_kin.sm2 || sE3<m_kin.sm3) {
    msg_Out()<<"ERROR in "<<METHOD<<" : "<<std::endl
	     <<"    energies/masses mismatch : "
	     <<sE1<<" vs. "<<m_kin.sm1<<", "
	     <<sE2<<" vs. "<<m_kin.sm2<<", "
	     <<sE3<<" vs. "<<m_kin.sm3<<"."<<std::endl;

    return false;
  }
  double ss12 = 1.-m_kin.x3+m_kin.sm3_2;
  double ss23 = 1.-m_kin.x1+m_kin.sm1_2;
  double ss13 = 1.-m_kin.x2+m_kin.sm2_2;
  double sp1  = sqrt(sqr(sE1)-m_kin.sm1_2);
  double sp2  = sqrt(sqr(sE2)-m_kin.sm2_2);
  double sp3  = sqrt(sqr(sE3)-m_kin.sm3_2);

  double cos13((2.*sE1*sE3+(m_kin.sm1_2+m_kin.sm3_2)-ss13)/(2.*sp1*sp3)),
    /*cos12((2.*E1*E2+sqr(m1)+sqr(m2)-s12)/(2.*p1*p2))*/;
    
  if (/*dabs(cos12)>1.01 ||*/ dabs(cos13)>1.01) {
    msg_Out()<<"ERROR in "<<METHOD<<" : |cos(angle)| = "<<cos13<<">1. from "
	     <<"2E1E3/2p1p3 = "<<(sE1*sE3/sp1/sp3)<<std::endl
	     <<"  E,p1 = "<<sE1<<", "<<sp1<<"; E,p3 = "<<sE3<<", "<<sp3<<" with "
	     <<"  sm1_2 = "<<m_kin.sm1_2<<", sm2_2 = "<<m_kin.sm2_2<<","
	     <<" sm3_2 = "<<m_kin.sm3_2<<", "<<std::endl
	     <<"   and ss12, ss13, ss23 = "
	     <<ss12<<"("<<sqr(m_kin.sm1+m_kin.sm2)<<"), "
	     <<ss13<<"("<<sqr(m_kin.sm1+m_kin.sm3)<<"), "
	     <<ss23<<"("<<sqr(m_kin.sm2+m_kin.sm3)<<"), "<<std::endl
	     <<"   check this (mass/mom ratios) : "
	     <<((sE1*sE1-sp1*sp1)/m_kin.sm1_2)<<", "
	     <<((sE2*sE2-sp2*sp2)/m_kin.sm2_2)<<", "
	     <<((sE3*sE3-sp3*sp3)/m_kin.sm3_2)<<std::endl; 
    return false;
  }
  //if (cos12>1.) cos12=1.; else if (cos12<-1.) cos12=-1.;
  if       (cos13>1.) cos13=1.; 
  else if (cos13<-1.) cos13=-1.;
  
  double theta13(acos(cos13)); //,theta12(acos(cos12));
  double psi(FixPsi(m_kin.x1,m_kin.x3,theta13)),chi(psi+theta13);
  double phi = 2.*M_PI*ran.Get();

  m_mom1 = m_kin.W*Vec4D(sE1,sp1*sin(psi)*sin(phi),-sp1*sin(psi)*cos(phi),sp1*cos(psi));
  m_mom3 = m_kin.W*Vec4D(sE3,sp3*sin(chi)*sin(phi),-sp3*sin(chi)*cos(phi),sp3*cos(chi));
  m_mom2 = Vec4D(m_kin.W,0.,0.,0.)-m_mom1-m_mom3;

  m_booster.Boost(momsplit);
  m_booster.Boost(momspect);
  m_rotator.Rotate(momsplit);
  m_rotator.Rotate(momspect);
  msg_Out()<<METHOD<<" cms-check : "<<std::endl
	   <<momsplit<<" + "<<momspect<<" = "<<m_kin.W<<" ---> "<<std::endl
	   <<m_mom1<<" + "<<m_mom2<<" + "<<m_mom3<<std::endl;

  m_rotator.RotateBack(m_mom1);
  m_rotator.RotateBack(m_mom2);
  m_rotator.RotateBack(m_mom3);
  m_booster.BoostBack(m_mom1);
  m_booster.BoostBack(m_mom2);
  m_booster.BoostBack(m_mom3);

  msg_Out()<<METHOD<<" : "<<std::endl
	   <<"  "<<p_split->m_mom<<" + "<<p_spect->m_mom<<" ---> "<<std::endl
	   <<"  "<<m_mom1<<" + "<<m_mom2<<" + "<<m_mom3<<std::endl
	   <<"  for x_{1,3} = "<<m_kin.x1<<", "<<m_kin.x3
	   <<" and masses = "<<m_kin.sm1<<", "<<m_kin.sm2<<","
	   <<" and "<<m_kin.sm3<<"."<<std::endl;

  p_out1 = new Proto_Particle(m_kin.flav.Bar(),m_mom2,'l');
  p_out2 = new Proto_Particle(m_kin.flav,m_mom3,'l');
  control::s_AHAprotoparticles+=2;
  
  return true;
}

double Dipole_Splitter::FixPsi(const double x1,const double x3,const double theta13) {
  switch (int(m_recoils)) {
  case 1:
    return sqr(x3)/(sqr(x1)+sqr(x3)) * (M_PI-theta13);
  case 2:
    if ((1.-x1)>(1.-x3)) return (M_PI-theta13);
    return 0.;
  case 99:
    msg_Error()<<"ERROR in Dipole_Splitter::FixPsi :"<<std::endl
	       <<"   Recoil strategy not known: "<<int(m_recoils)<<"."<<std::endl
	       <<"   Continue with default (spectator remains)."<<std::endl;
  case 0:
  default:
    break;
  }
  return (M_PI-theta13);
}
