#include "Jet_Finder.H"
#include "Run_Parameter.H"
#include "Message.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ATOOLS;


Jet_Finder::~Jet_Finder() {
  if (p_frame)   delete [] p_frame;
  if (p_value)   delete [] p_value;
}

void Jet_Finder::Init(const Vec4D * p)
{
  PROFILE_HERE;
  if (m_nin==2) {
    switch (m_type) {
    case 4 : 
      return;
    case 3 : {
      msg.Error()<<"Jet_Finder::Init : process-type "<<m_type
		 <<" not implemented yet !"<<std::endl;
      return;
    }
    case 2 : {
      msg.Error()<<"Jet_Finder::Init : process-type "<<m_type
		 <<" not implemented yet !"<<std::endl;
      return;
    }  
    case 1 : {
      m_sprime   = (p[0]+p[1]).Abs2();
      m_cms_boost = Poincare(p[0]+p[1]);
      return;
    }  
    case 0 : return;
    default :
      msg.Error()<<"This process-type is unknown!"<<std::endl;
    }
  }
}

double * Jet_Finder::ActualValue() { 
  return p_value; 
}

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Jet_Finder::Jet_Finder(double _ycut,int _type=1) : 
  m_ycut(_ycut), m_type(_type) , m_jet_alg(1), p_value(NULL), p_frame(NULL)
{
  rpa.gen.SetYcut(_ycut);

  m_name    = std::string("Jetfinder");

  m_ene     = rpa.gen.Ecms()/2;
  m_sprime  = m_s = sqr(2.*m_ene); 
  m_smin    = m_ycut * m_s;
  m_smax    = m_s;

  m_shower_pt2 = m_ycut * m_s;
  p_value   = new double[1];

  m_sel_log = new Selector_Log(m_name);
}

bool Jet_Finder::ConstructJets(Particle_List * pl, double y_res, bool final_only) 
{
  PROFILE_HERE;
  std::vector<Vec4D>   momsout;

  Vec4D   momsin[2];
  Flavour flavsin[2];
  if (!final_only) {
    for (int i=0;i<2;i++) {
      momsin[i]  = (*pl)[i]->Momentum();
      flavsin[i] = (*pl)[i]->Flav();
    }
    if ( (flavsin[0].Strong()) || (flavsin[1].Strong()) || (m_type != 1) ) {
      if (m_type==1) m_type=4;  // assume hadron hadron
    }

    Init(momsin);
    BoostInFrame(momsout);

    // delete first two
    pl->erase(pl->begin());
    pl->erase(pl->begin());
  }
  else {
    if (rpa.gen.Beam1().Strong() || rpa.gen.Beam2().Strong() || rpa.gen.Beam1().IsHadron() || rpa.gen.Beam2().IsHadron()) 
      m_type=4;
    flavsin[0]=rpa.gen.Beam1();
    flavsin[1]=rpa.gen.Beam1();

    momsin[0]=Vec4D::ZVEC;
    momsin[0]*=(rpa.gen.Ecms()*0.5);
    momsin[1]=Vec4D(momsin[0][0],-1.*Vec3D(momsin[0]));
  }

  // remove everything not to cluster and create momentum list
  for (Particle_List::iterator it=pl->begin(); it!=pl->end();) {
    //    if ((*parts)[i]->Flav().Strong()) {
    if (!(*it)->Flav().IsLepton()) {
      momsout.push_back((*it)->Momentum());
      ++it;
    }
    else {
      it=pl->erase(it);
    }
  }

  // Cluster vectors untill y_res reached!
  for (;;) {
    int j,k;
    double yij=YminKt(momsin,flavsin,momsout,j,k);
    if (yij>y_res) break;

    if (j<0) {
      //      momsin[j+2] += momsout[k]; // *AS*   ??!!!
    }
    else {
      momsout[j] += momsout[k];
    }
    for (int i=k;i<momsout.size()-1;i++) momsout[i] = momsout[i+1];
    momsout.pop_back();
    for (int i=k;i<pl->size()-1;i++) (*pl)[i] = (*pl)[i+1];
    pl->pop_back();
  }


  // create "complete new particle list"
  int j=0;
  for (Particle_List::iterator it=pl->begin(); it!=pl->end();++it,++j) {
    (*it)= new Particle(*it);
    (*it)->SetFlav(Flavour(kf::jet));
    (*it)->SetMomentum(momsout[j]);
  }

  return true;
}


bool Jet_Finder::ConstructJets(const Particle_List * parts,
			       const std::vector<int> & jets,std::vector<double> & lastys,bool final_only) 
{
  PROFILE_HERE;
  //  std::cout<<" in Jet_Finder::ConstructJets with "<<parts->size()<<" particles "<<std::endl;

  std::vector<Vec4D>   momsout;
  Vec4D   momsin[2];
  Flavour flavsin[2];
  if (!final_only) {
    for (int i=0;i<2;i++) {
      momsin[i]  = (*parts)[i]->Momentum();
      flavsin[i] = (*parts)[i]->Flav();
    }
    if ( (flavsin[0].Strong()) || (flavsin[1].Strong()) || (m_type != 1) ) {
      if (m_type==1) m_type=4;  // assume hadron hadron
    }

    for (int i=2;i<parts->size();i++) {
      //    if ((*parts)[i]->Flav().Strong()) {
      if (!(*parts)[i]->Flav().IsLepton()) {
	momsout.push_back((*parts)[i]->Momentum());
      }
    }
    Init(momsin);
    BoostInFrame(momsout);
  }
  else {
    if (rpa.gen.Beam1().Strong() || rpa.gen.Beam2().Strong() || rpa.gen.Beam1().IsHadron() || rpa.gen.Beam2().IsHadron()) 
      m_type=4;
    flavsin[0]=rpa.gen.Beam1();
    flavsin[1]=rpa.gen.Beam1();

    momsin[0]=Vec4D::ZVEC;
    momsin[0]*=(rpa.gen.Ecms()*0.5);
    momsin[1]=Vec4D(momsin[0][0],-1.*Vec3D(momsin[0]));

    for (int i=0;i<parts->size();i++) {
      //    if ((*parts)[i]->Flav().Strong()) {
      if (!(*parts)[i]->Flav().IsLepton()) {
	momsout.push_back((*parts)[i]->Momentum());
      }
    }

  }

  bool ordered = 1;
  while ((momsout.size()<=jets[lastys.size()]) && (lastys.size()<jets.size())) {
    //    std::cout<<" A: "<<momsout.size()<<" <=? "<<jets[lastys.size()]<<std::endl;
    lastys.push_back(-1.);
  }
  while (momsout.size()>jets.back()) {
    //    std::cout<<" B: "<<momsout.size()<<" >? "<<jets.back()<<std::endl;
    if (!ConstructJetSystem(momsin,flavsin,momsout,jets,lastys)) ordered = 0;
  }
  return ordered;
}

bool Jet_Finder::ConstructJetSystem(Vec4D * momsin,Flavour * flavsin,std::vector<Vec4D> & momsout,
				    std::vector<int> jets,std::vector<double> & lastys) 
{
  PROFILE_HERE;
  //  std::cout<<" in Jet_Finder::ConstructJetSystem with "<<momsout.size()<<" momenta "<<std::endl;
  int j,k;
  bool ordered = 1;
  // Calculate ymin and store for comparison
  lastys.push_back(YminKt(momsin,flavsin,momsout,j,k));
  if (lastys.size()>1) {
    if (lastys.back() > lastys[lastys.size()-2]) ordered = 0; 
  }
  // Erase previous y if not already jets.
  if ((momsout.size() > jets[0]) && (lastys.size()>1)) {
    lastys.front() = lastys.back();
    lastys.pop_back();
  }
  // Cluster vectors.
  if (j<0) {
    momsin[j+2] += momsout[k];
  }
  else {
    momsout[j] += momsout[k];
  }
  for (int i=k;i<momsout.size()-1;i++) momsout[i] = momsout[i+1];
  momsout.pop_back();

//   std::cout<<"============================================================"<<std::endl;
//   for (int i=0;i<momsout.size();++i) std::cout<<i<<" : "<<momsout[i]<<std::endl;

  return ordered;
}

double Jet_Finder::YminKt(Vec4D * momsin,Flavour * flavsin,std::vector<Vec4D> momsout,
			   int & j1,int & k1)
{
  PROFILE_HERE;
  double ymin = 2.;
  j1=-3; k1=-3;
  double pt2jk,pt2j,pt2k;
  for (int j=0;j<momsout.size();j++) {
    if (m_type>=3) {
      pt2j = (sqr(momsout[j][1]) + sqr(momsout[j][2]));
      if (pt2j < ymin*m_s) {
	ymin = pt2j/m_s;
	k1   = j;
	if (momsout[j][3]*momsin[0][3] > 0.) j1 = -2;
	                                else j1 = -1;
      } 
      for (int k=j+1;k<momsout.size();k++) {
	pt2k  = (sqr(momsout[k][1]) + sqr(momsout[k][2]));
	pt2jk = 2.*Min(pt2j,pt2k) * (Coshyp(DEta12(momsout[j],momsout[k])) - 
				     CosDPhi12(momsout[j],momsout[k]));
	if (pt2jk<ymin*m_s) {
	  ymin = pt2jk/m_s;
	  j1 = j;k1 = k;
	}
      }
    }
    else {
      for (int k=j+1;k<momsout.size();k++) {
	pt2jk  = 2.*sqr(Min(momsout[j][0],momsout[k][0]))*(1.-DCos12(momsout[j],momsout[k]));
	if (pt2jk<ymin*m_sprime) {
	  ymin = pt2jk/m_sprime;
	  j1 = j;k1 = k;
	}
      }
    }
  }

  if (j1==-3) {
    j1=0;
    k1=1;
  }
  return ymin;
}




/*---------------------------------------------------------------------

  Special form - flavours etc are known, will operate on momenta only.

  --------------------------------------------------------------------- */


Jet_Finder::Jet_Finder(int _n,Flavour * _fl,double _ycut,int _jetalg,int _type) : 
  m_ycut(_ycut), m_jet_alg(_jetalg), m_type(_type), p_value(NULL), p_frame(NULL) 
{
  rpa.gen.SetYcut(_ycut);

  m_name = std::string("Jetfinder");
  m_fl   = _fl;
  m_n    = _n;
  if (m_type==0) { m_nin = 1; m_nout = _n-1; }
          else { m_nin = 2; m_nout = _n-2; }

  
  p_frame = new Vec4D[m_nin];
  if (m_nin==1) {
    m_ene        = m_fl[0].Mass();
    m_sprime = m_s = sqr(m_ene); 
    p_frame[0] = Vec4D(m_ene,0.,0.,0.);
    m_cms_boost = Poincare(p_frame[0]);
  }
  if (m_nin==2) {
    if((m_type>=3) || (m_type==1)) {
      m_ene      = rpa.gen.Ecms()/2.;
      m_sprime   = m_s = sqr(2.*m_ene); 
	p_frame[0] = Vec4D(m_ene,0.,0., sqrt(sqr(m_ene)-sqr(m_fl[0].Mass())));
	p_frame[1] = Vec4D(m_ene,0.,0.,-sqrt(sqr(m_ene)-sqr(m_fl[1].Mass())));
    }    
    if (m_type==3) m_cms_boost = Poincare(p_frame[0]+p_frame[1]);
  }
  
  m_smin = m_ycut * m_s;
  m_smax = m_s;

  p_value   = new double[1];

  m_sel_log = new Selector_Log(m_name);
}

bool Jet_Finder::Trigger(const Vec4D * p)
{
  PROFILE_HERE;
  // create copy
  Vec4D * moms = new Vec4D[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) moms[i]=p[i];

  Init(moms);

  BoostInFrame(moms);

  int    j,k;
  bool   trigger = 1;
  double ymin;

  switch (m_jet_alg) {
  case 1 : 
    ymin = YminKt(moms,j,k); 
    if (ymin < m_ycut) {
      trigger = 0;
    } 
    break;
  default : 
    msg.Error()<<"No jet algorithm specified in Jet_Finder. return 0."<<std::endl;
    trigger = 0;
    break;
  }
  BoostBack(moms);

  delete [] moms;

  p_value[0] = ymin;
  return (1-m_sel_log->Hit(1-trigger));
}

double Jet_Finder::MTij2(Vec4D p1,Vec4D p2)
{
  double mt12_2;
  if (m_type>=2) {
    double pt1_2  = sqr(p1[1]) + sqr(p1[2]);
    double mt1_2  = pt1_2 + dabs(p1.Abs2()); 
    double pt2_2  = sqr(p2[1]) + sqr(p2[2]);
    double mt2_2  = pt2_2 + dabs(p2.Abs2()); 
    if (IsZero(pt1_2/(pt1_2+pt2_2))) {
      mt12_2        = mt2_2;
    }
    else 
      mt12_2        = 2.*Min(mt1_2,mt2_2) * (Coshyp(DEta12(p1,p2)) - CosDPhi12(p1,p2));
  }
  else {
    mt12_2        = 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2));
  }
  return mt12_2;
}

/*
double Jet_Finder::PTij2(Vec4D p1,Vec4D p2)
{
  double pt12_2;
  if (m_type>=2) {
    double pt1_2  = sqr(p1[1]) + sqr(p1[2]); 
    double pt2_2  = sqr(p2[1]) + sqr(p2[2]); 
    if (IsZero(pt1_2/(pt1_2+pt2_2))) {
      pt12_2        = pt2_2;
    }
    else 
      pt12_2        = 2.*Min(pt1_2,pt2_2) * (Coshyp(DEta12(p1,p2)) - CosDPhi12(p1,p2));
  }
  else {
    pt12_2        = 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2));
  }
  return pt12_2;
}
*/

bool Jet_Finder::TwoJets(const Vec4D & p1) 
{
  if (m_type>=2) {
    if (sqr(p1[1]) + sqr(p1[2])  < m_shower_pt2 ) return 0;
  }
  else {
    msg.Out()<<"Jet_Finder::TwoJets(Vec4D &) "<<std::endl;
    msg.Out()<<"    JETFINDER METHOD STILL NOT IMPLEMENTED for mode "<<m_type<<std::endl;
  }
  return 1;
}

bool Jet_Finder::TwoJets(const Vec4D & _p1,const Vec4D & _p2)
{
  Vec4D p1=_p1;
  Vec4D p2=_p2;

  BoostInFrame(p1);
  BoostInFrame(p2);

  if (m_type>=2) {
    double pt1_2  = sqr(p1[1]) + sqr(p1[2]); 
    double pt2_2  = sqr(p2[1]) + sqr(p2[2]); 
    if (pt1_2  < m_shower_pt2 ) return 0;
    if (pt2_2  < m_shower_pt2 ) return 0;
    double pt12_2 = 2.*Min(pt1_2,pt2_2) * (Coshyp(DEta12(p1,p2)) - CosDPhi12(p1,p2));
    if (pt12_2 < m_shower_pt2 ) return 0;
  }
  else {
    double pt12_2 = 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2));
    if (pt12_2 < m_shower_pt2 ) return 0;
  }

  BoostBack(p1);
  BoostBack(p2);
  return 1;
}


bool Jet_Finder::TwoJets(double & E2,double & z,double & costheta,bool mode)
{
  double pt12_2;
  if (mode == 1) {
    pt12_2 = -1000.;
    msg.Debugging()<<"JETFINDER METHOD STILL NOT IMPLEMENTED !!!!"<<std::endl;
  }
  else {
    pt12_2 = 2.*E2*sqr(Min(z,1.- z))*(1.-costheta);
    if (pt12_2 < m_shower_pt2 ) return 0;
  }
  return 1;
}

void Jet_Finder::BuildCuts(Cut_Data * cuts) 
{
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = m_fl[i].Mass();
    if (m_fl[i].Strong()) {                
      /* 
	 minimal energies : 
	 either   :  E^2 > kt^2 > y_cut s      
	             (hadron-hadron collisions)
	 or       :  4 min{E_i^2,E_j^2} > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	             ycut s' > ycut s_min   
	             (lepton-lepton collisions)
      */
      if (m_type>=2) cuts->energymin[i] = Max(sqrt(m_ycut * m_s),cuts->energymin[i]);
              else cuts->energymin[i] = Max(sqrt(m_ycut * m_smin/4.),cuts->energymin[i]);
      for (int j=i+1; j<m_nin+m_nout; ++j) {
	if (m_fl[j].Strong()) {                
	  /* 
	     minimal scut :
	     either   :  s_ij = 2 E_i E_j (1-cos(ij)) > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	                 m_ycut s' > ycut s_min   
	                 (lepton-lepton collisions)
	     or       :  similarly .... have to think ...
               	         (hadron-hadron collisions)
 
	  */
	  cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],m_ycut*m_smin);
	}
      }
    }
  }
}

void   Jet_Finder::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = m_fl[i].Mass();
    // strong interacting particles first.
    if (m_fl[i].Strong()) {                
      /* 
	 minimal energies : 
	 either   :  E^2 > kt^2 > y_cut s      
	             (hadron-hadron collisions)
	 or       :  4 min{E_i^2,E_j^2} > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	             ycut s' > ycut s_min   
	             (lepton-lepton collisions)
      */
      if (m_type>=2) cuts->energymin[i] = Max(sqrt(m_ycut * m_s),cuts->energymin[i]);
              else cuts->energymin[i] = Max(sqrt(m_ycut * m_smin/4.),cuts->energymin[i]);
      for (int j=i+1; j<m_nin+m_nout; ++j) {
	if (m_fl[j].Strong()) {                
	  /* 
	     minimal scut :
	     either   :  s_ij = 2 E_i E_j (1-cos(ij)) > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	                 ycut s' > ycut s_min   
	                 (lepton-lepton collisions)
	     or       :  similarly .... have to think ...
               	         (hadron-hadron collisions)
 
	  */
	  if (m_type>=2) cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],m_ycut*m_smin);
	          else cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],m_ycut*sprime);
	}
      }
    }
  }
}

double Jet_Finder::YminKt(Vec4D * p,int & j1,int & k1)
{
  PROFILE_HERE;
  double ymin = 2.;
  double pt2jk,pt2j,pt2k;
  for (int j=m_nin;j<m_n;j++) {
    if (m_fl[j].Strong()) {
      if (m_type>=3) {
	pt2j = (sqr(p[j][1]) + sqr(p[j][2]));
	if (pt2j < ymin*m_s) {
	  ymin = pt2j/m_s;
	  j1   = j;
	  if (p[j][3]*p[0][3] > 0.) k1 = 0;
	                       else k1 = 1;
	} 
	for (int k=j+1;k<m_n;k++) {
	  if (m_fl[k].Strong()) {
	    pt2k  = (sqr(p[k][1]) + sqr(p[k][2]));
	    pt2jk = 2.*Min(pt2j,pt2k) * (Coshyp(DEta12(p[j],p[k])) - CosDPhi12(p[j],p[k]));
	    if (pt2jk<ymin*m_s) {
	      ymin = pt2jk/m_s;j1 = j;k1 = k;
	    }
	  }
	}
      }
      else {
	for (int k=j+1;k<m_n;k++) {
	  if (m_fl[k].Strong()) {
	    pt2jk  = 2.*sqr(Min(p[j][0],p[k][0]))*(1.-DCos12(p[j],p[k]));
	    if (pt2jk<ymin*m_sprime) {
	      ymin = pt2jk/m_sprime;j1 = j;k1 = k;
	    }
	  }
	}
      }
    }
  }
  return ymin;
}

double Jet_Finder::DEta12(Vec4D & p1,Vec4D & p2)
{
  // eta1,2 = -log(tan(theta_1,2)/2)   
  double c1=p1[3]/Vec3D(p1).Abs();
  double c2=p2[3]/Vec3D(p2).Abs();
  return  0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
}

double Jet_Finder::CosDPhi12(Vec4D & p1,Vec4D & p2)
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return (p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2);
}

double Jet_Finder::DCos12(Vec4D & p1,Vec4D & p2)
{
  return Vec3D(p1)*Vec3D(p2)/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

double Jet_Finder::Coshyp(double x) 
{
  //  return 0.5*(std::exp(x)+std::exp(-x));
  // return 0.5*(exp(x)+exp(-x));
  return cosh(x);
}

void Jet_Finder::BoostInFrame(Vec4D * p)
{
  for (int i=0;i<m_n;i++) m_cms_boost.Boost(p[i]);
}

void Jet_Finder::BoostBack(Vec4D * p)
{
  for (int i=0;i<m_n;i++) m_cms_boost.BoostBack(p[i]);
}

void Jet_Finder::BoostInFrame(std::vector<Vec4D> p)
{
  for (int i=0;i<p.size();i++) m_cms_boost.Boost(p[i]);
}

void Jet_Finder::BoostBack(std::vector<Vec4D> p)
{
  for (int i=0;i<p.size();i++) m_cms_boost.BoostBack(p[i]);
}

void Jet_Finder::BoostInFrame(Vec4D & p)
{
  m_cms_boost.Boost(p);
}

void Jet_Finder::BoostBack(Vec4D & p)
{
  m_cms_boost.BoostBack(p);
}
