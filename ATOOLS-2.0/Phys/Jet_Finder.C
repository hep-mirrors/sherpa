#include "Jet_Finder.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;


Jet_Finder::~Jet_Finder() {
  if (p_frame)   delete [] p_frame;
  if (p_value)   delete [] p_value;
}

void Jet_Finder::Init(const AMATOOLS::Vec4D * p)
{
  if (m_nin==2) {
    switch (m_type) {
    case 4 : 
      // msg.Tracking()<<"Jet_Finder::Init : process-type hadron-hadron."<<std::endl;
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
      // msg.Tracking()<<"Jet_Finder::Init : process-type lepton-lepton."<<std::endl;
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

  General form - flavours etc are unknown, will operate on a Parton_List.

  --------------------------------------------------------------------- */

Jet_Finder::Jet_Finder(double _ycut,int _type=1) : 
  p_value(NULL), p_frame(NULL), m_ycut(_ycut), m_type(_type) , m_jet_alg(1)
{
  AORGTOOLS::msg.Debugging()<<"Initialize the Jet_Finder : "<<std::endl
			    <<"   Jetalg = "<<m_jet_alg<<", type = "<<m_type
			    <<", ycut = "<<m_ycut<<std::endl;

  m_name    = std::string("Jetfinder");

  m_ene     = AORGTOOLS::rpa.gen.Ecms()/2;
  m_sprime  = m_s = sqr(2.*m_ene); 
  m_smin    = m_ycut * m_s;
  m_smax    = m_s;

  p_value   = new double[1];

  m_sel_log = new Selector_Log(m_name);
}

bool Jet_Finder::ConstructJets(const APHYTOOLS::Parton_List * parts,
			       const std::vector<int> & jets,std::vector<double> & lastys) {
  std::vector<Vec4D>   momsout;
  Vec4D   momsin[2];
  Flavour flavsin[2];
  for (int i=0;i<2;i++) {
    momsin[i]  = (*parts)[i]->Momentum();
    flavsin[i] = (*parts)[i]->Flav();
  }
  if ( (flavsin[0].Strong()) || (flavsin[1].Strong()) || (m_type != 1) ) {
    //    msg.Out()<<"WARNING: Jet_Finder::ConstructJets : This is not completely tested yet! "<<std::endl;
    if (m_type==1) m_type=4;  // assume hadron hadron
  }

  for (int i=2;i<parts->size();i++) {
    if ((*parts)[i]->Flav().Strong()) {
      momsout.push_back((*parts)[i]->Momentum());
    }
  }
  Init(momsin);
  BoostInFrame(momsout);

  bool ordered = 1;
  //  while (momsout.size() <= jets[jets.front()-lastys.size()]) {
  while ((momsout.size()<=jets[lastys.size()]) && (lastys.size()<jets.size())) {
    //    cout<<" i="<<lastys.size()<<" moms="<<momsout.size()<<" "<<jets[lastys.size()]<<std::endl;
    lastys.push_back(-1.);
  }
  while (momsout.size()>jets.back()) {
    if (!ConstructJetSystem(momsin,flavsin,momsout,jets,lastys)) ordered = 0;
  }
  return ordered;
}

bool Jet_Finder::ConstructJetSystem(Vec4D * momsin,Flavour * flavsin,std::vector<Vec4D> & momsout,
				    std::vector<int> jets,std::vector<double> & lastys) {
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

  return ordered;
}

double Jet_Finder::YminKt(Vec4D * momsin,Flavour * flavsin,std::vector<Vec4D> momsout,
			   int & j1,int & k1)
{
  double ymin = 2.;
  j1=-3; k1=-3;
  double pt2jk,pt2j,pt2k,pt2min;
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
				     DPhi12(momsout[j],momsout[k]));
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
    std::cout<<" problem: "<<ymin<<std::endl;
    for (int k=0; k<momsout.size();++k)
      std::cout<<k<<" "<<momsout[k]<<std::endl;
    j1=0;
    k1=1;
  }
  return ymin;
}




/*---------------------------------------------------------------------

  Special form - flavours etc are known, will operate on momenta only.

  --------------------------------------------------------------------- */


Jet_Finder::Jet_Finder(int _n,Flavour * _fl,double _ycut,int _jetalg,int _type) : 
  p_value(NULL),p_frame(NULL),m_ycut(_ycut), m_jet_alg(_jetalg), m_type(_type) 
{
  AORGTOOLS::msg.Debugging()<<"Initialize the <"<<_n<<"> Jet_Finder : "<<std::endl;
  AORGTOOLS::msg.Debugging()<<"   Jetalg = "<<m_jet_alg<<", type = "<<m_type;
  AORGTOOLS::msg.Debugging()<<", ycut = "<<m_ycut<<std::endl;

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
      m_ene      = AORGTOOLS::rpa.gen.Ecms()/2.;
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

bool Jet_Finder::Trigger(const AMATOOLS::Vec4D * p)
{
  // create copy
  Vec4D * moms = new Vec4D[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) moms[i]=p[i];

  Init(moms);

  //msg.Out()<<"In Jet_Finder::Trigger. Before boosting :"<<std::endl;
  //for (int i=0;i<m_nin+m_nout;i++) msg.Out()<<"   "<<i<<" th mom : "<<p[i]<<std::endl;

  BoostInFrame(moms);

  //msg.Out()<<"In Jet_Finder::Trigger. Before loop :"<<std::endl;
  //for (int i=0;i<m_nin+m_nout;i++) msg.Out()<<"   "<<i<<" th mom : "<<moms[i]<<std::endl;

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

  //msg.Out()<<"In Jet_Finder::Trigger. After loop :"<<std::endl;
  //for (int i=0;i<m_nin+m_nout;i++) msg.Out()<<"   "<<i<<" th mom : "<<moms[i]<<std::endl;
  delete [] moms;

  p_value[0] = ymin;
  return (1-m_sel_log->Hit(1-trigger));
}

double Jet_Finder::PTij(AMATOOLS::Vec4D p1,AMATOOLS::Vec4D p2)
{
  double pt12_2;
  if (m_type>=2) {
    double pt1_2  = sqr(p1[1]) + sqr(p1[2]); 
    double pt2_2  = sqr(p2[1]) + sqr(p2[2]); 
    if (IsZero(pt1_2/(pt1_2+pt2_2))) {
      pt12_2        = pt2_2;
    }
    else 
      pt12_2        = 2.*Min(pt1_2,pt2_2) * (Coshyp(DEta12(p1,p2)) - DPhi12(p1,p2));
    //    cout<<" ptij = "<<pt12_2<<endl;
  }
  else {
    pt12_2        = 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2));
  }
  return pt12_2;
}


bool Jet_Finder::TwoJets(const AMATOOLS::Vec4D & p1) 
{
  if (m_type>=2) {
    if (sqr(p1[1]) + sqr(p1[2])  < m_s * m_ycut) return 0;
  }
  else {
    msg.Out()<<"Jet_Finder::TwoJets(Vec4D &) "<<std::endl;
    msg.Out()<<"    JETFINDER METHOD STILL NOT IMPLEMENTED for mode "<<m_type<<std::endl;
  }
  return 1;
}

bool Jet_Finder::TwoJets(const AMATOOLS::Vec4D & _p1,const AMATOOLS::Vec4D & _p2)
{
  Vec4D p1=_p1;
  Vec4D p2=_p2;

  BoostInFrame(p1);
  BoostInFrame(p2);

  if (m_type>=2) {
    double pt1_2  = sqr(p1[1]) + sqr(p1[2]); 
    double pt2_2  = sqr(p2[1]) + sqr(p2[2]); 
    if (pt1_2  < m_s * m_ycut) return 0;
    if (pt2_2  < m_s * m_ycut) return 0;
    double pt12_2 = 2.*Min(pt1_2,pt2_2) * (Coshyp(DEta12(p1,p2)) - DPhi12(p1,p2));
    if (pt12_2 < m_s * m_ycut) return 0;
  }
  else {
    double pt12_2 = 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2));
    if (pt12_2 < m_sprime * m_ycut) return 0;
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
    /*
    std::cout<<" ycut="<<m_ycut<<std::endl;
    std::cout<<" sprime="<<m_sprime<<std::endl;
    std::cout<<" z     ="<<z<<std::endl;
    std::cout<<" E2    ="<<E2<<std::endl;
    std::cout<<" costheta    ="<<costheta<<std::endl;
    */
    pt12_2 = 2.*E2*sqr(Min(z,1.- z))*(1.-costheta);
    if (pt12_2 < m_sprime * m_ycut) return 0;
  }
  return 1;
}

void Jet_Finder::BuildCuts(Cut_Data * cuts) 
{
  msg.Debugging()<<"In Jet_Finder::BuildCuts"<<std::endl;
  // Loop over final state particles.
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

double Jet_Finder::YminKt(AMATOOLS::Vec4D * p,int & j1,int & k1)
{
  double ymin = 2.;
  double pt2jk,pt2j,pt2k,pt2min;
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
	    pt2jk = 2.*Min(pt2j,pt2k) * (Coshyp(DEta12(p[j],p[k])) - DPhi12(p[j],p[k]));
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

double Jet_Finder::DEta12(AMATOOLS::Vec4D & p1,AMATOOLS::Vec4D & p2)
{
  // eta1,2 = -log(tan(theta_1,2)/2)   
  // => eta1 - eta2 
  //  = -log(tan1/2)+log(tan2/2) = log(tan2/tan1) = log(pt2/pt1 * pl1/pl2)
  return log(sqrt( (sqr(p2[1])+sqr(p2[2]))/(sqr(p1[1])+sqr(p1[2])) ) * dabs(p1[3]/p2[3]));
}

double Jet_Finder::DPhi12(AMATOOLS::Vec4D & p1,AMATOOLS::Vec4D & p2)
{
  // cos(phi1-phi2) = cos(phi1) cos(phi2) + sin(phi1) sin(phi2)
  //                = (p1_x p2_x + p1_y p2_y)/(p1_z * p2_z)
  return (p1[1]*p2[1]+p1[2]*p2[2])/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

double Jet_Finder::DCos12(AMATOOLS::Vec4D & p1,AMATOOLS::Vec4D & p2)
{
  return Vec3D(p1)*Vec3D(p2)/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

double Jet_Finder::Coshyp(double x) 
{
  //  return 0.5*(std::exp(x)+std::exp(-x));
  // return 0.5*(exp(x)+exp(-x));
  return cosh(x);
}

void Jet_Finder::BoostInFrame(AMATOOLS::Vec4D * p)
{
  for (int i=0;i<m_n;i++) m_cms_boost.Boost(p[i]);
}

void Jet_Finder::BoostBack(AMATOOLS::Vec4D * p)
{
  for (int i=0;i<m_n;i++) m_cms_boost.BoostBack(p[i]);
}

void Jet_Finder::BoostInFrame(std::vector<AMATOOLS::Vec4D> p)
{
  for (int i=0;i<p.size();i++) m_cms_boost.Boost(p[i]);
}

void Jet_Finder::BoostBack(std::vector<AMATOOLS::Vec4D> p)
{
  for (int i=0;i<p.size();i++) m_cms_boost.BoostBack(p[i]);
}

void Jet_Finder::BoostInFrame(AMATOOLS::Vec4D & p)
{
  m_cms_boost.Boost(p);
}

void Jet_Finder::BoostBack(AMATOOLS::Vec4D & p)
{
  m_cms_boost.BoostBack(p);
}
