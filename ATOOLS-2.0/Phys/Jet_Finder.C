#include "Jet_Finder.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;


Jet_Finder::~Jet_Finder() {
  //  if (sel_log) delete [] sel_log;
  if (frame)   delete [] frame;
}

void Jet_Finder::Init(const AMATOOLS::vec4d * p)
{
  if (nin==2) {
    switch (type) {
    case 4 : 
      // msg.Tracking()<<"Jet_Finder::Init : process-type hadron-hadron."<<std::endl;
      return;
    case 3 : {
      msg.Error()<<"Jet_Finder::Init : process-type "<<type
		 <<" not implemented yet !"<<std::endl;
      return;
    }
    case 2 : {
      msg.Error()<<"Jet_Finder::Init : process-type "<<type
		 <<" not implemented yet !"<<std::endl;
      return;
    }  
    case 1 : {
      // msg.Tracking()<<"Jet_Finder::Init : process-type lepton-lepton."<<std::endl;
      sprime   = (p[0]+p[1]).abs2();
      CMSBoost = Poincare(p[0]+p[1]);
      return;
    }  
    case 0 : return;
    default :
      msg.Error()<<"This process-type is unknown!"<<std::endl;
    }
  }
}

double * Jet_Finder::ActualValue() { 
  return value; 
}

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Parton_List.

  --------------------------------------------------------------------- */

Jet_Finder::Jet_Finder(double _ycut,int _type=1) : 
  ycut(_ycut), type(_type) , jet_alg(1)
{
  AORGTOOLS::msg.Debugging()<<"Initialize the Jet_Finder : "<<std::endl
			    <<"   Jetalg = "<<jet_alg<<", type = "<<type
			    <<", ycut = "<<ycut<<std::endl;

  name    = std::string("Jetfinder");

  E       = AORGTOOLS::rpa.gen.Ecms()/2;
  sprime  = s = sqr(2.*E); 
  smin    = ycut * s;
  smax    = s;

  value   = new double[1];

  sel_log = new Selector_Log(name);
}

bool Jet_Finder::ConstructJets(const APHYTOOLS::Parton_List * parts,
			       const std::vector<int> & jets,std::vector<double> & lastys) {
  std::vector<vec4d>   momsout;
  vec4d   momsin[2];
  Flavour flavsin[2];
  for (int i=0;i<2;i++) {
    momsin[i]  = (*parts)[i]->momentum();
    flavsin[i] = (*parts)[i]->flav();
  }
  if ( (flavsin[0].strong()) || (flavsin[1].strong()) || (type != 1) ) {
    msg.Error()<<"Jet_Finder::ConstructJets : This is not implemented yet !!!!!"<<std::endl;
    abort();
  }

  for (int i=2;i<parts->size();i++) {
    if ((*parts)[i]->flav().strong()) {
      momsout.push_back((*parts)[i]->momentum());
    }
  }
  Init(momsin);
  BoostInFrame(momsout);

  bool ordered = 1;
  //  while (momsout.size() <= jets[jets.front()-lastys.size()]) {
  while ((momsout.size() <= jets[lastys.size()])&&(lastys.size()<jets.size())) {
    //    cout<<" i="<<lastys.size()<<" moms="<<momsout.size()<<" "<<jets[lastys.size()]<<std::endl;
    lastys.push_back(-1.);
  }
  while (momsout.size() > jets.back()) {
    if (!ConstructJetSystem(momsin,flavsin,momsout,jets,lastys)) ordered = 0;
  }
  return ordered;
}

bool Jet_Finder::ConstructJetSystem(vec4d * momsin,Flavour * flavsin,std::vector<vec4d> & momsout,
				    std::vector<int> jets,std::vector<double> & lastys) {
  int j,k;
  bool ordered = 1;
  // Calculate ymin and store for comparison
  lastys.push_back(ymin_kt(momsin,flavsin,momsout,j,k));
  if (lastys.size()>1) {
    if (lastys.back() > lastys[lastys.size()-2]) ordered = 0; 
  }
  // Erase previous y if not already jets.
  if ((momsout.size() > jets[0]) && (lastys.size()>1)) {
    lastys.front() = lastys.back();
    lastys.pop_back();
  }
  // Cluster vectors.
  momsout[j] += momsout[k];
  for (int i=k;i<momsout.size()-1;i++) momsout[i] = momsout[i+1];
  momsout.pop_back();

  return ordered;
}

double Jet_Finder::ymin_kt(vec4d * momsin,Flavour * flavsin,std::vector<vec4d> momsout,
			   int & j1,int & k1)
{
  double ymin = 2.;
  j1=-1; k1=-1;
  double pt2jk,pt2j,pt2k,pt2min;
  for (int j=0;j<momsout.size();j++) {
    if (type>=3) {
      pt2j = (sqr(momsout[j][1]) + sqr(momsout[j][2]));
      if (pt2j < ymin*s) {
	ymin = pt2j/s;
	j1   = j;
	if (momsout[j][3]*momsin[0][3] > 0.) k1 = 0;
	                                else k1 = 1;
      } 
      for (int k=j+1;k<momsout.size();k++) {
	pt2k  = (sqr(momsout[k][1]) + sqr(momsout[k][2]));
	pt2jk = 2.*Min(pt2j,pt2k) * (coshyp(deta12(momsout[j],momsout[k])) - 
				     dphi12(momsout[j],momsout[k]));
	if (pt2jk<ymin*s) {
	  ymin = pt2jk/s;
	  j1 = j;k1 = k;
	}
      }
    }
    else {
      for (int k=j+1;k<momsout.size();k++) {
	pt2jk  = 2.*sqr(Min(momsout[j][0],momsout[k][0]))*(1.-dcos12(momsout[j],momsout[k]));
	if (pt2jk<ymin*sprime) {
	  ymin = pt2jk/sprime;
	  j1 = j;k1 = k;
	}
      }
    }
  }

  if (j1==-1) {
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
  ycut(_ycut), jet_alg(_jetalg), type(_type) 
{
  AORGTOOLS::msg.Debugging()<<"Initialize the <"<<_n<<"> Jet_Finder : "<<std::endl;
  AORGTOOLS::msg.Debugging()<<"   Jetalg = "<<jet_alg<<", type = "<<type;
  AORGTOOLS::msg.Debugging()<<", ycut = "<<ycut<<std::endl;

  name = std::string("Jetfinder");
  fl   = _fl;
  n    = _n;
  if (type==0) { nin = 1; nout = n-1; }
          else { nin = 2; nout = n-2; }

  
  frame = new vec4d[nin];
  if (nin==1) {
    E        = fl[0].mass();
    sprime   = s = sqr(E); 
    frame[0] = vec4d(E,0.,0.,0.);
    CMSBoost = Poincare(frame[0]);
  }
  if (nin==2) {
    if((type>=3) || (type==1)) {
      E        = AORGTOOLS::rpa.gen.Ecms()/2.;
      sprime   = s = sqr(2.*E); 
	frame[0] = vec4d(E,0.,0., sqrt(sqr(E)-sqr(fl[0].mass())));
	frame[1] = vec4d(E,0.,0.,-sqrt(sqr(E)-sqr(fl[1].mass())));
    }    
    if (type==3) CMSBoost = Poincare(frame[0]+frame[1]);
  }
  
  smin = ycut * s;
  smax = s;

  value   = new double[1];

  sel_log = new Selector_Log(name);
}

bool Jet_Finder::Trigger(const AMATOOLS::vec4d * p)
{
  // create copy
  vec4d * moms = new vec4d[nin+nout];
  for (int i=0;i<nin+nout;i++) moms[i]=p[i];

  Init(moms);

  //msg.Out()<<"In Jet_Finder::Trigger. Before boosting :"<<std::endl;
  //for (int i=0;i<nin+nout;i++) msg.Out()<<"   "<<i<<" th mom : "<<p[i]<<std::endl;

  BoostInFrame(moms);

  //msg.Out()<<"In Jet_Finder::Trigger. Before loop :"<<std::endl;
  //for (int i=0;i<nin+nout;i++) msg.Out()<<"   "<<i<<" th mom : "<<moms[i]<<std::endl;

  int    j,k;
  bool   trigger = 1;
  double ymin;

  switch (jet_alg) {
  case 1 : 
    ymin = ymin_kt(moms,j,k); 
    if (ymin < ycut) {
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
  //for (int i=0;i<nin+nout;i++) msg.Out()<<"   "<<i<<" th mom : "<<moms[i]<<std::endl;
  delete [] moms;

  value[0] = ymin;
  return (1-sel_log->Hit(1-trigger));
}

double Jet_Finder::PTij(AMATOOLS::vec4d p1,AMATOOLS::vec4d p2)
{
  double pt12_2;
  if (type>=2) {
    double pt1_2  = sqr(p1[1]) + sqr(p1[2]); 
    double pt2_2  = sqr(p2[1]) + sqr(p2[2]); 
    pt12_2        = 2.*Min(pt1_2,pt2_2) * (coshyp(deta12(p1,p2)) - dphi12(p1,p2));
  }
  else {
    pt12_2        = 2.*sqr(Min(p1[0],p2[0]))*(1.-dcos12(p1,p2));
  }
  return pt12_2;
}


bool Jet_Finder::TwoJets(const AMATOOLS::vec4d & p1) 
{
  if (type>=2) {
    if (sqr(p1[1]) + sqr(p1[2])  < s * ycut) return 0;
  }
  return 1;
}

bool Jet_Finder::TwoJets(AMATOOLS::vec4d & p1,AMATOOLS::vec4d & p2)
{
  BoostInFrame(p1);
  BoostInFrame(p2);

  if (type>=2) {
    double pt1_2  = sqr(p1[1]) + sqr(p1[2]); 
    double pt2_2  = sqr(p2[1]) + sqr(p2[2]); 
    if (pt1_2  < s * ycut) return 0;
    if (pt2_2  < s * ycut) return 0;
    double pt12_2 = 2.*Min(pt1_2,pt2_2) * (coshyp(deta12(p1,p2)) - dphi12(p1,p2));
    if (pt12_2 < s * ycut) return 0;
  }
  else {
    double pt12_2 = 2.*sqr(Min(p1[0],p2[0]))*(1.-dcos12(p1,p2));
    if (pt12_2 < sprime * ycut) return 0;
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
    std::cout<<" ycut="<<ycut<<std::endl;
    std::cout<<" sprime="<<sprime<<std::endl;
    std::cout<<" z     ="<<z<<std::endl;
    std::cout<<" E2    ="<<E2<<std::endl;
    std::cout<<" costheta    ="<<costheta<<std::endl;
    */
    pt12_2 = 2.*E2*sqr(Min(z,1.- z))*(1.-costheta);
    if (pt12_2 < sprime * ycut) return 0;
  }
  return 1;
}

void Jet_Finder::BuildCuts(Cut_Data * cuts) 
{
  msg.Debugging()<<"In Jet_Finder::BuildCuts"<<std::endl;
  // Loop over final state particles.
  for (int i=nin; i<nin+nout; ++i) {
    cuts->energymin[i] = fl[i].mass();
    // strong interacting particles first.
    if (fl[i].strong()) {                
      /* 
	 minimal energies : 
	 either   :  E^2 > kt^2 > y_cut s      
	             (hadron-hadron collisions)
	 or       :  4 min{E_i^2,E_j^2} > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	             ycut s' > ycut s_min   
	             (lepton-lepton collisions)
      */
      if (type>=2) cuts->energymin[i] = Max(sqrt(ycut * s),cuts->energymin[i]);
              else cuts->energymin[i] = Max(sqrt(ycut * smin/4.),cuts->energymin[i]);
      for (int j=i+1; j<nin+nout; ++j) {
	if (fl[j].strong()) {                
	  /* 
	     minimal scut :
	     either   :  s_ij = 2 E_i E_j (1-cos(ij)) > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	                 ycut s' > ycut s_min   
	                 (lepton-lepton collisions)
	     or       :  similarly .... have to think ...
               	         (hadron-hadron collisions)
 
	  */
	  cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],ycut*smin);
	}
      }
    }
  }
}

void   Jet_Finder::UpdateCuts(double sprime,double y,Cut_Data * cuts) {
  for (int i=nin; i<nin+nout; ++i) {
    cuts->energymin[i] = fl[i].mass();
    // strong interacting particles first.
    if (fl[i].strong()) {                
      /* 
	 minimal energies : 
	 either   :  E^2 > kt^2 > y_cut s      
	             (hadron-hadron collisions)
	 or       :  4 min{E_i^2,E_j^2} > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	             ycut s' > ycut s_min   
	             (lepton-lepton collisions)
      */
      if (type>=2) cuts->energymin[i] = Max(sqrt(ycut * s),cuts->energymin[i]);
              else cuts->energymin[i] = Max(sqrt(ycut * smin/4.),cuts->energymin[i]);
      for (int j=i+1; j<nin+nout; ++j) {
	if (fl[j].strong()) {                
	  /* 
	     minimal scut :
	     either   :  s_ij = 2 E_i E_j (1-cos(ij)) > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	                 ycut s' > ycut s_min   
	                 (lepton-lepton collisions)
	     or       :  similarly .... have to think ...
               	         (hadron-hadron collisions)
 
	  */
	  if (type>=2) cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],ycut*smin);
	          else cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],ycut*sprime);
	}
      }
    }
  }
}

double Jet_Finder::ymin_kt(AMATOOLS::vec4d * p,int & j1,int & k1)
{
  double ymin = 2.;
  double pt2jk,pt2j,pt2k,pt2min;
  for (int j=nin;j<n;j++) {
    if (fl[j].strong()) {
      if (type>=3) {
	pt2j = (sqr(p[j][1]) + sqr(p[j][2]));
	if (pt2j < ymin*s) {
	  ymin = pt2j/s;
	  j1   = j;
	  if (p[j][3]*p[0][3] > 0.) k1 = 0;
	                       else k1 = 1;
	} 
	for (int k=j+1;k<n;k++) {
	  if (fl[k].strong()) {
	    pt2k  = (sqr(p[k][1]) + sqr(p[k][2]));
	    pt2jk = 2.*Min(pt2j,pt2k) * (coshyp(deta12(p[j],p[k])) - dphi12(p[j],p[k]));
	    if (pt2jk<ymin*s) {
	      ymin = pt2jk/s;j1 = j;k1 = k;
	    }
	  }
	}
      }
      else {
	for (int k=j+1;k<n;k++) {
	  if (fl[k].strong()) {
	    pt2jk  = 2.*sqr(Min(p[j][0],p[k][0]))*(1.-dcos12(p[j],p[k]));
	    if (pt2jk<ymin*sprime) {
	      ymin = pt2jk/sprime;j1 = j;k1 = k;
	    }
	  }
	}
      }
    }
  }
  return ymin;
}

double Jet_Finder::deta12(AMATOOLS::vec4d & p1,AMATOOLS::vec4d & p2)
{
  // eta1,2 = -log(tan(theta_1,2)/2)   
  // => eta1 - eta2 
  //  = -log(tan1/2)+log(tan2/2) = log(tan2/tan1) = log(pt2/pt1 * pl1/pl2)
  return log(sqrt( (sqr(p2[1])+sqr(p2[2]))/(sqr(p1[1])+sqr(p1[2])) ) * dabs(p1[3]/p2[3]));
}

double Jet_Finder::dphi12(AMATOOLS::vec4d & p1,AMATOOLS::vec4d & p2)
{
  // cos(phi1-phi2) = cos(phi1) cos(phi2) + sin(phi1) sin(phi2)
  //                = (p1_x p2_x + p1_y p2_y)/(p1_z * p2_z)
  return (p1[1]*p2[1]+p2[1]*p2[2])/(vec3d(p1).abs()*vec3d(p2).abs());
}

double Jet_Finder::dcos12(AMATOOLS::vec4d & p1,AMATOOLS::vec4d & p2)
{
  return vec3d(p1)*vec3d(p2)/(vec3d(p1).abs()*vec3d(p2).abs());
}

double Jet_Finder::coshyp(double x) 
{
  //  return 0.5*(std::exp(x)+std::exp(-x));
  // return 0.5*(exp(x)+exp(-x));
  return cosh(x);
}

void Jet_Finder::BoostInFrame(AMATOOLS::vec4d * p)
{
  for (int i=0;i<n;i++) CMSBoost.Boost(p[i]);
}

void Jet_Finder::BoostBack(AMATOOLS::vec4d * p)
{
  for (int i=0;i<n;i++) CMSBoost.BoostBack(p[i]);
}

void Jet_Finder::BoostInFrame(std::vector<AMATOOLS::vec4d> p)
{
  for (int i=0;i<p.size();i++) CMSBoost.Boost(p[i]);
}

void Jet_Finder::BoostBack(std::vector<AMATOOLS::vec4d> p)
{
  for (int i=0;i<p.size();i++) CMSBoost.BoostBack(p[i]);
}

void Jet_Finder::BoostInFrame(AMATOOLS::vec4d & p)
{
  CMSBoost.Boost(p);
}

void Jet_Finder::BoostBack(AMATOOLS::vec4d & p)
{
  CMSBoost.BoostBack(p);
}
