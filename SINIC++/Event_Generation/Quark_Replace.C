#include "SINIC++/Event_Generation/Quark_Replace.H"
#include "SINIC++/Event_Generation/Ladder.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace SINIC;
using namespace ATOOLS;

Quark_Replace::Quark_Replace() :
  p_as(static_cast<MODEL::Strong_Coupling *>
       (MODEL::s_model->GetScalarFunction(std::string("strong_cpl"))))
{ 
}

Quark_Replace::~Quark_Replace() {
}

/////////////////////////////////////////////////////////////////////////////
// replace gluon pairs with quark pairs according to MEs
/////////////////////////////////////////////////////////////////////////////

void Quark_Replace::ReplaceWithQuarks(Ladder * ladder){
  LadderMap::iterator yiter1(ladder->GetEmissionsBegin());
  LadderMap::iterator yiter2(yiter1),ytmp(yiter1),yfirst(yiter1); yiter2++;
  LadderMap::iterator ylast(ladder->GetEmissionsEnd()); ylast--;
  TPropList::iterator piter(ladder->GetPropsBegin());
  int flav(0);
  Ladder_Particle out1,out2,newout1,newout2;
  ATOOLS::Flavour flout1,flout2;
  ATOOLS::Flavour inflav1, inflav2;
  ATOOLS::Flow    inflow1, inflow2;
  ATOOLS::Vec4D outmom1,outmom2;
  double y1,y2,pq;
  bool glgl(false),glgl2(false),replace(false),replacefirst(false);
  bool colourok(true),ggedgeok(true),qgedgeok(true),closedindex,quark(false),antiq(false);

  ladder->GetInFlows(inflow1,inflow2);
  ladder->GetInFlavours(inflav1,inflav2);
  yiter1->second.m_flav = inflav1;
  ylast->second.m_flav = inflav2;

  if (inflow1.Code(1)==0 && inflow1.Code(2)==0) closedindex = false;
  else if (inflow1.Code(1)!=inflow2.Code(2) && inflow1.Code(2)!=inflow2.Code(1)) closedindex = false;
  else closedindex = true;

  if ((inflav1.IsQuark() && !inflav1.IsAnti()) || 
      (inflav2.IsQuark() && !inflav2.IsAnti())) quark=true;
  if ((inflav1.IsQuark() && inflav1.IsAnti()) || 
      (inflav2.IsQuark() && inflav2.IsAnti())) antiq=true;


  // check if incoming colours allow for triplet exchanges
  if (inflow1.Code(1)==0 && inflow2.Code(2)==0) {
    if (inflav1.IsGluon() && inflav2.IsGluon() && ladder->Size()==2) {
      colourok = false;
    }
    else if ((inflav1.IsGluon() && inflav2.IsQuark()) || (inflav2.IsGluon() && inflav1.IsQuark())) {
      ggedgeok = false;
    }
  }
  else {
    if (inflav1.IsGluon() && inflav2.IsGluon() && ladder->Size()==2 && 
	inflow1.Code(1)!=inflow2.Code(2) && inflow1.Code(2)!=inflow2.Code(1)) {
      colourok = false;
    }
    else if (((inflav1.IsGluon() && inflav2.IsQuark()) || (inflav2.IsGluon() && inflav1.IsQuark())) && 
	inflow1.Code(1)!=inflow2.Code(2) && inflow1.Code(2)!=inflow2.Code(1)) {
      ggedgeok = false;
    }
    else if (inflav1.IsQuark() && !inflav1.IsAnti() && inflow1.Code(1)==inflow2.Code(2)) {
      qgedgeok = false;
    }
    else if (inflav1.IsQuark() && inflav1.IsAnti() && inflow1.Code(2)==inflow2.Code(1)) {
      qgedgeok = false;
    }
    else if (inflav2.IsQuark() && !inflav2.IsAnti() && inflow2.Code(1)==inflow1.Code(2)) {
      qgedgeok = false;
    }
    else if (inflav2.IsQuark() && inflav2.IsAnti() && inflow2.Code(2)==inflow1.Code(1)) {
      qgedgeok = false;
    }
  }

  msg_Debugging()<<METHOD<<"   initial ladder: "<<std::endl<<(*ladder)<<std::endl
		 <<METHOD<<"   initial flavours: "<<inflav1<<"/"<<inflav2<<"."<<std::endl;
  msg_Debugging()<<"   incoming colour indices: ( "<<inflow1.Code(1)<<" , "<<inflow1.Code(2)<<" ) and ( "
               <<inflow2.Code(1)<<" , "<<inflow2.Code(2)<<" ) ; colours o.k.: "<<colourok
               <<" , gg edge o.k.: "<<ggedgeok<<" , qg edge o.k.: "<<qgedgeok
               <<" , closed index: "<<closedindex<<std::endl;
  while (yiter1 != ylast) {
    if(piter->m_col==colour_type::octet) {
      if (colourok) {
        pq = GetTripletWeight(yiter1,piter,flav,glgl);
        if ((yiter1==yfirst || yiter2==ylast) && glgl && !ggedgeok) pq = 0.;
        if ((yiter1==yfirst || yiter2==ylast) && !glgl && !qgedgeok) pq = 0.;
        if (yiter2==ylast && glgl && replacefirst) pq = 0.;
      }
      else {
 	pq = 0.;
      }
      if(glgl) {
        if((closedindex && quark) ||
           (!closedindex && antiq) ||
           (!quark && !antiq)) {
          flout1 = ATOOLS::Flavour(flav,0);
          flout2 = ATOOLS::Flavour(flav,1);
        }
        else {
          flout1 = ATOOLS::Flavour(flav,1);
          flout2 = ATOOLS::Flavour(flav,0);
        }
      }
      else {
        flout1 = yiter2->second.m_flav;
        flout2 = yiter1->second.m_flav;
      }
      if (yiter2 != ylast) {
        piter++;
        if(piter->m_col==colour_type::octet) {
          pq *= (1.-GetTripletWeight(yiter2,piter,flav,glgl2));
        }
        piter--;
      }
    }
    else pq = 0.;
      // replace outgoing gluons with quarks
    if(ATOOLS::ran->Get()<pq){
      if (yiter1==yfirst && glgl) replacefirst=true;
      out1 = yiter1->second;
      out2 = yiter2->second;
      outmom1=out1.m_mom;
      outmom2=out2.m_mom;
      if(glgl && flout1.Mass()!=0.){
        replace = MassiveKinematics(flout1.Mass(),outmom1,outmom2,y1,y2,yiter1,yfirst,ylast);
        if(replace){
	  ytmp=yiter2;
	  ytmp++;
	  ladder->DeleteRapidity(yiter1);
	  ladder->DeleteRapidity(yiter2);
	  ladder->AddRapidity(y1);
	  ladder->AddRapidity(y2);
	  yiter1=ytmp;
	  yiter2=ytmp;
	  yiter1--;
	  yiter1--;
	  yiter2--;
	  if (ytmp==ladder->GetEmissionsEnd()) ylast=yiter2;
	}
      }
      if(replace || (flout1.Mass()==0. && flout2.Mass()==0.)) {
        out1.m_mom=outmom1;
        out2.m_mom=outmom2;
        out1.m_flav=flout1;
        out2.m_flav=flout2;
		if(out1.p_part!=NULL) {
	  out1.p_part->SetMomentum(outmom1);
	  out2.p_part->SetMomentum(outmom2);
	  out1.p_part->SetFlav(flout1);
	  out2.p_part->SetFlav(flout2);
        }
        ladder->SetParticle(yiter1,out1);
        ladder->SetParticle(yiter2,out2);
        piter->m_col=colour_type::triplet;
        yiter1++;
        if(yiter2 != ylast) yiter2++;
	piter++;
      }
    }
    if(yiter1 != ylast) yiter1++;
    if(yiter2 != ylast) yiter2++;
    piter++;
  }
  msg_Debugging()<<METHOD<<"   final ladder: "<<std::endl<<(*ladder);
}




double Quark_Replace::GetTripletWeight(LadderMap::iterator yiter1,
                 const TPropList::iterator piter, int & flav,bool & glgl){
  double mu2(0.),meg(0.),meq(0.),pq(0.);
  int nf;
  Ladder_Particle out1,out2;
  ATOOLS::Flavour flout1,flout2,flin1,flin2;
  double hats, hatt, hatu;
  ATOOLS::Flavour inflav1, inflav2;
  double qt2,qt02,deltay,y1,y2;
  LadderMap::iterator yiter2;
 
  yiter2=yiter1;
  yiter2++;
  out1   = yiter1->second;
  out2   = yiter2->second;
  flin1 = out1.m_flav;
  flin2 = out2.m_flav;
  // figure out flavour combination
  if(flin1.IsGluon() && flin2.IsGluon()){
    //calculate cross sections
    glgl=true;
    ConstructMandelstams(out1.m_mom,out2.m_mom,mu2,hats,hatt,hatu);
    if (hats>4.*pow(ATOOLS::Flavour(kf_b).HadMass(),2)) nf=5;
    else if (hats>4.*pow(ATOOLS::Flavour(kf_c).HadMass(),2)) nf=4;
    else if (hats>4.*pow(ATOOLS::Flavour(kf_s).HadMass(),2)) nf=3;
    else nf=2;
    flav = int(ATOOLS::ran->Get()*nf)+1;
    flout1=ATOOLS::Flavour(flav,false);
    flout2=ATOOLS::Flavour(flav,true);
    meg = m_me2s(flin1,flin2,flin1,flin2,hats,hatt,hatu);
    meq = nf*m_me2s(flin1,flin2,flout1,flout2,hats,hatt,hatu);
    if((meq<0.) || (meg<0.)){
      msg_Out()<<"We have a problem here:"<<std::endl;
      msg_Out()<<"gluonic ME: "<<meg<<" ;  quarkic ME:  "<<meq<<std::endl;
    }
  }
  else if ((flin1.IsGluon() && flin2.IsQuark()) || (flin2.IsGluon() && flin1.IsQuark())){
    glgl=false;
    flout1=flin2;
    flout2=flin1;
    //calculate cross sections
    qt2    = piter->m_qt2;
    qt02   = piter->m_q02;
    deltay = yiter2->first-yiter1->first;
    meg    = pow(qt02/(qt2+qt02),3.*(*p_as)(qt2)*deltay/M_PI);
    meq    = pow(qt02/(qt2+qt02),0.5);
  }
  else if (flin1.IsQuark() && flin2.IsQuark()){
    glgl=false;
    meq=0.;
    meg=1.;
  }
  else {
    glgl=false;
    msg_Out()<<"What is this?"<<std::endl;
    msg_Out()<<"in flavours "<<flin1<<", "<<flin2<<std::endl;
    exit(1);
  }
  return meq/(meq+meg);
}

bool Quark_Replace::MassiveKinematics(const double m, 
         ATOOLS::Vec4D & outmom1, ATOOLS::Vec4D & outmom2, double & y1, double & y2,
         const LadderMap::iterator yiter1, const LadderMap::iterator yfirst, 
         const LadderMap::iterator ylast){
  LadderMap::iterator yiter2,ytmp,ytmp2;
  ATOOLS::Vec4D dir0(0.,0.,0.,1.);
  bool keep;

  ATOOLS::Poincare boost(outmom1+outmom2);
  boost.Boost(outmom1);
  boost.Boost(outmom2);
  ATOOLS::Poincare rot(outmom1,dir0);
  outmom1=rot*outmom1;
  outmom2=rot*outmom2;
  outmom1[3]=sqrt(outmom1[0]*outmom1[0]-m*m);
  outmom2[3]=-sqrt(outmom2[0]*outmom2[0]-m*m);
  rot.RotateBack(outmom1);
  rot.RotateBack(outmom2);
  boost.BoostBack(outmom1);
  boost.BoostBack(outmom2);
  y1=0.5*log((outmom1[0]+outmom1[3])/(outmom1[0]-outmom1[3]));
  y2=0.5*log((outmom2[0]+outmom2[3])/(outmom2[0]-outmom2[3]));
  //check if we change the order of emissions
  yiter2=yiter1;
  yiter2++;
  if ((yiter1==yfirst) && (yiter2==ylast)) {
    keep=true;
  }
  else if (yiter1==yfirst){
    ytmp=yiter2;
    ytmp++;
    if (y2<ytmp->first) keep=true;
    else keep=false;
  }
  else if (yiter2==ylast){
    ytmp=yiter1;
    ytmp--;
    if (y1>ytmp->first) keep=true;
    else keep=false;
  }
  else{
    ytmp=yiter1;
    ytmp--;
    ytmp2=yiter2;
    ytmp2++;
    if ((y1>ytmp->first) && (y2<ytmp2->first)) keep=true;
    else keep=false;
  }
  return keep;
}

void Quark_Replace::ConstructMandelstams(const ATOOLS::Vec4D & out1,const ATOOLS::Vec4D & out2,
		     const double & mu2,double & hats, double & hatt, double & hatu)
{
  double y1(out1.Y()), y2(out2.Y());
  double ystar((y1-y2)/2.), ybar((y1+y2)/2.), coshy(cosh(ystar)), tanhy(tanh(ystar));
  hats  = (out1+out2).Abs2();
  hatt  = -hats/2.*(1.+tanhy)-mu2; 
  hatu  = -hats/2.*(1.-tanhy)-mu2;
  double pperp(sqrt(hats/ATOOLS::sqr(2.*coshy)));
  hats += 2.*mu2;
  //x1    = 2.*pperp/m_Ecms*exp(ybar)*coshy;
  //x2    = 2.*pperp/m_Ecms*exp(-ybar)*coshy;
}  
