#include "Jet_Finder.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Exception.H"
#include "Blob.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ATOOLS;


/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Jet_Finder::Jet_Finder(const double _ycut,const int _type,const bool _pt_def) : 
  m_ycut(_ycut), m_delta_r(1.), m_type(_type) , m_pt_def(_pt_def),
  m_value(0.), p_frame(NULL)
{
  m_name       = std::string("Jetfinder");
  m_ene        = rpa.gen.Ecms()/2.;
  m_sprime     = m_s = sqr(2.*m_ene); 
  m_smin       = m_ycut * m_s;
  m_smax       = m_s;
  m_shower_pt2 = m_ycut * m_s;

  m_sel_log    = new Selector_Log(m_name);
}

/*---------------------------------------------------------------------

  Special form - flavours etc are known, will operate on momenta only.

  --------------------------------------------------------------------- */


Jet_Finder::Jet_Finder(const int _n,Flavour * _fl,const double _ycut,const int _type,const bool _pt_def) : 
  m_ycut(_ycut), m_delta_r(1.), m_type(_type), m_pt_def(_pt_def),
  m_value(0.), p_frame(NULL)
{
  m_name = std::string("Jetfinder");
  m_fl   = _fl;
  m_n    = _n;
  if (m_type==0) { m_nin = 1; m_nout = _n-1; }
            else { m_nin = 2; m_nout = _n-2; }

  
  p_frame = new Vec4D[m_nin];
  if (m_nin==1) {
    m_ene       = m_fl[0].Mass();
    m_sprime    = m_s = sqr(m_ene); 
    p_frame[0]  = Vec4D(m_ene,0.,0.,0.);
    m_cms_boost = Poincare(p_frame[0]);
  }
  else if (m_nin==2) {
    if((m_type>=3) || (m_type==1)) {
      m_ene      = rpa.gen.Ecms()/2.;
      m_sprime   = m_s = sqr(2.*m_ene); 
      p_frame[0] = Vec4D(m_ene,0.,0., sqrt(sqr(m_ene)-sqr(m_fl[0].Mass())));
      p_frame[1] = Vec4D(m_ene,0.,0.,-sqrt(sqr(m_ene)-sqr(m_fl[1].Mass())));
      if (m_type==3) m_cms_boost = Poincare(p_frame[0]+p_frame[1]);
    }    
    else if (m_type==2) {
      m_ene      = rpa.gen.Ecms()/2.;
      m_sprime   = m_s = sqr(2.*m_ene);
    }
  }
  
  m_smin    = m_ycut * m_s;
  m_smax    = m_s;
  m_sel_log = new Selector_Log(m_name);
}

Jet_Finder::~Jet_Finder() {
  if (p_frame)   delete [] p_frame;
}


/*----------------------------------------------------------------------------------

  Constructing jets, mainly for phase space cuts.

  ----------------------------------------------------------------------------------*/


void Jet_Finder::Init(const Vec4D * p)
{
  PROFILE_HERE;
  if (m_nin==2) {
    switch (m_type) {
    case 4 : return;
    case 3 : {
      msg.Error()<<"Jet_Finder::Init : process-type "<<m_type
		 <<" not implemented yet !"<<std::endl;
      return;
    }
    case 2 : {
      //Initialize the Breit frame
      int lepton(0);
      if (m_fl[0].Strong()) lepton=1;
      int hadron(lepton==1?0:1);
      Vec4D q(p[lepton]);
            
      for (int i=m_nin;i<m_nin+m_nout;i++)
	if (m_fl[i]==m_fl[lepton]) q-=p[i];
      Vec4D store(q);

      double x(-q.Abs2()/(2.*p[hadron]*q)); 
      Vec4D pp(2.*x*p[hadron]+q);
      double gamma(pp[0]/pp.Abs());
      Vec3D eta(Vec3D(pp)/pp.Abs());
      
      m_cms_boost = Poincare(Vec4D(gamma,eta));
      m_cms_boost.Boost(q);
      m_zrot      = Poincare(-1.*q,Vec4D::ZVEC);
      m_zrot.Rotate(q);
      BoostBack(q);
      
      //checks
      if (dabs(q*pp)>1.e-10) 
	msg.Error()<<" ERROR: Jet_Finder::Init could not initialize Breit frame correctly (1) : "
		   <<dabs(q*pp)<<std::endl;
      
      bool check(true);
      for (int i=0;i<3;i++) 
	if (dabs((q[i]-store[i]))>1.e10) check = false; 
      if (!check) msg.Error()<<" ERROR: Jet_Finder::Init could not initialize Breit frame correctly (2) : "
			     <<q-store<<std::endl;
      
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

bool Jet_Finder::ConstructJets(Particle_List * pl, double y_res,int number, bool final_only) 
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
    
    // remove everything not to cluster and create momentum list
    for (Particle_List::iterator it=pl->begin(); it!=pl->end();) {
      if (!(*it)->Flav().IsLepton()) {
	momsout.push_back((*it)->Momentum());
	++it;
      }
      else {
	it=pl->erase(it);
      }
    }

    Vec4D * moms = new Vec4D[m_nin+momsout.size()];
    for (int i=0;i<m_nin;i++) 
      moms[i]=momsin[i];
    for (unsigned int i=m_nin;i<m_nin+momsout.size();i++) 
      moms[i]=momsout[i-m_nin];
    
    Init(moms);
    BoostInFrame(momsin);
    BoostInFrame(momsout);
    
    // delete first two
    pl->erase(pl->begin());
    pl->erase(pl->begin());
  }
  else {
    if (rpa.gen.Beam1().Strong()   || rpa.gen.Beam2().Strong() || 
	rpa.gen.Beam1().IsHadron() || rpa.gen.Beam2().IsHadron()) 
      m_type=4;
    flavsin[0]=rpa.gen.Beam1();
    flavsin[1]=rpa.gen.Beam1();
    
    if (m_type!=2) {
    momsin[0]=Vec4D::ZVEC;
    momsin[0]*=(rpa.gen.Ecms()*0.5);
    momsin[1]=Vec4D(momsin[0][0],-1.*Vec3D(momsin[0]));
    }
    // remove everything not to cluster and create momentum list
    for (Particle_List::iterator it=pl->begin(); it!=pl->end();) {
      if (!(*it)->Flav().IsLepton()) {
	momsout.push_back((*it)->Momentum());
	++it;
      }
      else {
	it=pl->erase(it);
      }
    }
    //DIS case
    if (m_type==2) {
      Vec4D * moms = new Vec4D[m_nin+momsout.size()];
      for (int i=0;i<m_nin;i++) 
	moms[i]=momsin[i];
      for (unsigned int i=m_nin;i<m_nin+momsout.size();i++) 
	moms[i]=momsout[i-m_nin];
      
      Init(moms);
      BoostInFrame(momsin);
      BoostInFrame(momsout);
    }
  }

  // Cluster vectors untill y_res or number reached!
  for (;;) {
    int j,k;
    double yij=YminKt(momsin,flavsin,momsout,j,k);
    if (yij>y_res) break;
    if ((int)momsout.size()<=number) break;
    
    if (j<0) {
      //      momsin[j+2] += momsout[k]; // *AS*   ??!!!
    }
    else {
      momsout[j] += momsout[k];
    }
    for (size_t i=k;i<momsout.size()-1;i++) momsout[i] = momsout[i+1];
    momsout.pop_back();
    for (size_t i=k;i<pl->size()-1;i++) (*pl)[i] = (*pl)[i+1];
    pl->pop_back();
  }
  
  
  // create "complete new particle list"
  int j=0;
  for (Particle_List::iterator it=pl->begin(); it!=pl->end();++it,++j) {
    (*it)= new Particle(**it);
    (*it)->SetFlav(Flavour(kf::jet));
    (*it)->SetMomentum(momsout[j]);
  }

  return true;
}

bool Jet_Finder::ConstructJets(const Particle_List * parts,
			       const std::vector<int> & jets,std::vector<double> & lastys,bool final_only) 
{
  PROFILE_HERE;

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

    for (size_t i=2;i<parts->size();i++) {
      if (!(*parts)[i]->Flav().IsLepton()) {
	momsout.push_back((*parts)[i]->Momentum());
      }
    }
    
    Vec4D * moms = new Vec4D[m_nin+momsout.size()];
    for (int i=0;i<m_nin;i++) 
      moms[i]=momsin[i];
    for (unsigned int i=m_nin;i<m_nin+momsout.size();i++) 
      moms[i]=momsout[i-m_nin];
    
    Init(moms);
    BoostInFrame(momsin);
    BoostInFrame(momsout);
  }
  else {
    if (rpa.gen.Beam1().Strong() || rpa.gen.Beam2().Strong() || 
	rpa.gen.Beam1().IsHadron() || rpa.gen.Beam2().IsHadron()) 
      m_type=4;
    flavsin[0]=rpa.gen.Beam1();
    flavsin[1]=rpa.gen.Beam1();

    if (m_type!=2) {
      momsin[0]=Vec4D::ZVEC;
      momsin[0]*=(rpa.gen.Ecms()*0.5);
      momsin[1]=Vec4D(momsin[0][0],-1.*Vec3D(momsin[0]));
    }
    for (size_t i=0;i<parts->size();i++) {
      if (!(*parts)[i]->Flav().IsLepton()) {
	momsout.push_back((*parts)[i]->Momentum());
      }
    }
    //DIS case
    if (m_type==2) {
      Vec4D * moms = new Vec4D[m_nin+momsout.size()];
      for (int i=0;i<m_nin;i++) 
	moms[i]=momsin[i];
      for (unsigned int i=m_nin;i<m_nin+momsout.size();i++) 
	moms[i]=momsout[i-m_nin];
      
      Init(moms);
      BoostInFrame(momsin);
      BoostInFrame(momsout);
    }
  }

  bool ordered = 1;
  while (((int)momsout.size()<=jets[lastys.size()]) && (lastys.size()<jets.size())) {
    lastys.push_back(-1.);
  }
  while ((int)momsout.size()>jets.back()) {
    if (!ConstructJetSystem(momsin,flavsin,momsout,jets,lastys)) ordered = 0;
  }
  return ordered;
}

bool Jet_Finder::ConstructJetSystem(Vec4D * momsin,Flavour * flavsin,std::vector<Vec4D> & momsout,
				    std::vector<int> jets,std::vector<double> & lastys) 
{
  PROFILE_HERE;
  int j,k;
  bool ordered = 1;
  // Calculate ymin and store for comparison
  lastys.push_back(YminKt(momsin,flavsin,momsout,j,k));
  if (lastys.size()>1) {
    if (lastys.back() > lastys[lastys.size()-2]) ordered = 0; 
  }
  // Erase previous y if not already jets.
  if (((int)momsout.size() > jets[0]) && (lastys.size()>1)) {
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
  for (size_t i=k;i<momsout.size()-1;i++) momsout[i] = momsout[i+1];
  momsout.pop_back();

  return ordered;
}

double Jet_Finder::YminKt(Vec4D * momsin,Flavour * flavsin,std::vector<Vec4D> momsout,int & j1,int & k1)
{
  PROFILE_HERE;
  double ymin = 2.;
  j1=-3; k1=-3;
  double pt2jk,pt2j,pt2k;
  for (size_t j=0;j<momsout.size();j++) {
    if (m_type>=3) {
      pt2j = (sqr(momsout[j][1]) + sqr(momsout[j][2]));
      if (pt2j < ymin*m_s) {
	ymin = pt2j/m_s;
	k1   = j;
	if (momsout[j][3]*momsin[0][3] > 0.) j1 = -2;
	                                else j1 = -1;
      }
      for (size_t k=j+1;k<momsout.size();k++) {
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
      if (m_type==2) {
	int hadron=m_fl[0].Strong()?0:1;
	pt2j = 2.*sqr(momsout[j][0])*(1.-DCos12(momsout[j],momsin[hadron]));
	if (pt2j < ymin*m_sprime) {
	  ymin = pt2j/m_sprime;
	  k1   = j;
	  //cluster to beam 
	  if (hadron) j1 = -1;
	         else j1 = -2;
	}
      }
      for (size_t k=j+1;k<momsout.size();k++) {
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

std::vector<int> ID(size_t id)
{
  std::vector<int> ids;
  for (size_t n(0);id>0;++n) {
    if (id&(1<<n)) {
      ids.push_back(n);
      id-=1<<n;
    }
  }
  return ids;
}

size_t ID(const std::vector<int> &ids)
{
  size_t id(0);
  for (size_t i(0);i<ids.size();++i) 
    if (ids[i]!=0) id+=1<<i;
  return id;
}

Flavour Jet_Finder::GetFlavour(std::string fl)
{
  bool bar(false);
  if (fl=="j") return Flavour(kf::jet);
  if (fl=="G") return Flavour(kf::gluon);
  if (fl.length()>1) {
    if (fl[fl.length()-1]=='b') {
      fl.erase(fl.length()-1,1);
      bar=true;
    }
    else if (fl[fl.length()-1]=='+') {
      fl[fl.length()-1]='-';
      bar=true;
    }
  }
  Flavour flav(kf_table.FromString(fl));
  if (flav.Kfcode()==kf::none) 
    THROW(critical_error,"No flavour for '"+fl+"'.");
  if (bar) flav=flav.Bar();
  return flav;
}

size_t Jet_Finder::FillCombinations(const std::string &name,size_t &cp,
				    const int fl)
{
  size_t sum(0);
  std::vector<int> pos;
  for (size_t i(0);i<name.length();++i) {
    if (name[i]=='[') {
      int open(1);
      for (size_t j(i+1);j<name.length();++j) {
	if (name[j]=='[') ++open;
	if (name[j]==']') --open;
	if (open==0) {
	  pos.push_back(FillCombinations(name.substr(i+1,j-i-1),cp,fl-1));
	  size_t sp(name.rfind('_',i));
	  m_flavs[pos.back()]=GetFlavour(name.substr(sp+1,i-sp-1));
	  sum=sum|pos.back();
	  i=j+2;
	  break;
	}
      }
    }
    else if (name[i]=='_' && name[i-1]!='_') {
      pos.push_back(1<<cp++);
      size_t sp(name.rfind('_',i-1));
      if (sp==std::string::npos) sp=0;
      else ++sp;
      m_flavs[pos.back()]=GetFlavour(name.substr(sp,i-sp));
      sum=sum|pos.back();
    }
  }
  if (name[name.length()-1]!=']') {
    pos.push_back(1<<cp++);
    size_t sp(name.rfind('_',name.length()-1));
    m_flavs[pos.back()]=GetFlavour(name.substr(sp+1,name.length()-sp-1));
    sum=sum|pos.back();
  }
  for (size_t i(0);i<pos.size();++i) {
    for (size_t j(i+1);j<pos.size();++j) {
      if (pos[i]>2 || pos[j]>2) {
	m_combs[pos[i]][pos[j]]=1;
	m_fills.push_back(std::pair<size_t,size_t>(pos[i],pos[j]));
      }
    }
//     if (fl<m_nin+m_nout) {
//       m_combs[pos[i]][sum]=1;
//       m_fills.push_back(std::pair<size_t,size_t>(pos[i],sum));
//     }
  }
  m_mcomb.push_back(pos);
  m_mcomb.back().push_back(sum);
  return sum;
}

void Jet_Finder::FillCombinations()
{
  if (m_combs.empty()) {
    if (m_procname=="") THROW(fatal_error,"Process name not set.");
    size_t size(1<<(m_nin+m_nout));
    m_moms.resize(size);
    m_flavs.resize(size);
    m_combs.resize(size,std::vector<int>(size,0));
    std::string name(m_procname.substr(m_procname.find('_')+1));
    name=name.substr(name.find('_')+1);
    size_t i(0);
    FillCombinations(name,i,m_nin+m_nout);
    if (msg.LevelIsDebugging()) {
      msg.Out()<<METHOD<<"(): Combinations for '"<<m_procname<<"' {\n";
      for (size_t i(0);i<m_combs.size();++i) {
	if (ID(m_combs[i])!=0) {
	  msg.Out()<<"  "<<ID(i)<<"["<<m_flavs[i]<<","
		   <<m_flavs[i].Strong()<<"] & {";
	  for (size_t j(0);j<m_combs[i].size();++j)
	    if (m_combs[i][j]!=0) msg.Out()<<" "<<ID(j)<<"["<<m_flavs[j]<<","
					   <<m_flavs[j].Strong()<<"]";
	  msg.Out()<<" }\n";
	}
      }
      msg.Out()<<"}\n";
      msg.Out()<<METHOD<<"(): Identified clusterings {\n";
      for (size_t j(0);j<m_fills.size();++j)
	msg.Out()<<" ["<<ID(m_fills[j].first)<<","
		     <<ID(m_fills[j].second)<<"]\n";
      msg.Out()<<"}\n";
      msg.Out()<<METHOD<<"(): Momentum combination {\n";
      for (size_t i(0);i<m_mcomb.size();++i) {
	msg.Out()<<"  "<<ID(m_mcomb[i].back())<<" -> {";
	for (size_t j(0);j<m_mcomb[i].size()-1;++j) 
	  msg.Out()<<" "<<ID(m_mcomb[i][j]);
	msg.Out()<<" }\n";
      }
      msg.Out()<<"}\n";
    }
  }
}

void Jet_Finder::PrepareMomList()
{
  for (size_t n(0);n<m_mcomb.size()-1;++n) {
    m_moms[m_mcomb[n].back()]=m_moms[m_mcomb[n].front()];
    for (size_t i(1);i<m_mcomb[n].size()-1;++i)
      m_moms[m_mcomb[n].back()]+=m_moms[m_mcomb[n][i]];
//     msg_Debugging()<<"p["<<ID(m_mcomb[n].back())<<"] = "
// 		   <<m_moms[m_mcomb[n].back()]
// 		   <<" ["<<m_flavs[m_mcomb[n].back()]<<"]\n";
  }
}

bool Jet_Finder::Trigger(const Vec4D * p)
{
  FillCombinations();
  PROFILE_HERE;
  // create copy
  for (int i=0;i<m_nin+m_nout;i++) {
    m_moms[1<<i]=p[i];
//     msg_Debugging()<<"p["<<i<<"] = "<<m_moms[1<<i]<<" ("<<m_flavs[1<<i]<<")\n";
  }

  Init(&m_moms.front());
  BoostInFrame(&m_moms.front());
  PrepareMomList();

  int    j,k;
  bool   trigger(true);
  double ymin=0.;
  ymin = YminKt(&m_moms.front(),j,k); 
  if (ymin < m_ycut) trigger = false;
  
  m_value = ymin;
  return (1-m_sel_log->Hit(1-trigger));
}

void Jet_Finder::BuildCuts(Cut_Data * cuts) 
{
  FillCombinations();
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = m_fl[i].SelMass();
    if (m_fl[i].Strong()) {                
      /* 
	 minimal energies : 
	 either   :  E^2 > kt^2 > y_cut s      
	             (hadron-hadron collisions)
	 or       :  4 min{E_i^2,E_j^2} > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	             ycut s' > ycut s_min   
	             (lepton-lepton collisions)
      */
      if (m_type>=2 && (m_combs[1<<i][1<<0] || m_combs[1<<i][1<<1])) {
	cuts->energymin[i] = Max(sqrt(1. * m_ycut * m_s),cuts->energymin[i]);
	if (m_type==4) {
	  cuts->cosmax[0][i] = cuts->cosmax[1][i] = cuts->cosmax[i][0] = cuts->cosmax[i][1] =  
	    Min(cuts->cosmax[0][i],sqrt(1.-4.*m_ycut));
	  cuts->etmin[i] = Max(sqrt(m_ycut * m_s),cuts->etmin[i]);
	}
	if (m_type==2) {
	  int hadron=m_fl[0].Strong()?0:1;
	  if (m_combs[1<<i][1<<hadron]) {
	    cuts->cosmax[hadron][i] = cuts->cosmax[i][hadron] = 
	      Min(cuts->cosmax[hadron][i],sqrt(1.-4.*m_ycut));
	    cuts->cosmin[hadron][i] = cuts->cosmin[i][hadron] = 
	      Max(cuts->cosmin[hadron][i],-sqrt(1.-4.*m_ycut));
	    cuts->etmin[i] = Max(sqrt(m_ycut * m_s),cuts->etmin[i]);
	  }
	}
      }
      else cuts->energymin[i] = Max(sqrt(m_ycut * m_smin/4.),cuts->energymin[i]);
      
      for (int j=i+1; j<m_nin+m_nout; ++j) {
	if (m_fl[j].Strong() && m_combs[1<<i][1<<j]) {
	  /* 
	     minimal scut :
	     either   :  s_ij = 2 E_i E_j (1-cos(ij)) > 2 min{E_i^2,E_j^2} (1-cos(ij)) > 
	     m_ycut s' > ycut s_min   
	     (lepton-lepton collisions)
	     or       :  similarly .... have to think ...
               	         (hadron-hadron collisions)
 
	  */
	  if (m_type>=2) cuts->scut[j][i] = cuts->scut[i][j] 
			   = Max(cuts->scut[i][j],sqr(m_delta_r)*m_ycut*m_s);
	  else cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],m_ycut*m_smin);
	}
      }
    }
  }
}

void   Jet_Finder::UpdateCuts(double sprime,double y,Cut_Data * cuts) 
{
  if (m_type>1) return;
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = Max(sqrt(m_ycut * sprime/4.),cuts->energymin[i]);
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
	cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],m_ycut*sprime);
      }
    }
  }
}

double Jet_Finder::YminKt(Vec4D * p,int & j1,int & k1)
{
  PROFILE_HERE;
  double ymin = 2.;
  double pt2jk,pt2j,pt2k;
  for (size_t ps(0);ps<m_fills.size();++ps) {
    int j(m_fills[ps].first), k(m_fills[ps].second);
    Vec4D pj(p[j]), pk(p[k]);
//     msg_Debugging()<<"test "<<ID(j)<<"["<<m_flavs[j]<<"] & "
// 		   <<ID(k)<<"["<<m_flavs[k]<<"]\n";
    if (j&k) {
      if (j>k) pj=-pj;
      else pk=-pk;
    }
    if (m_flavs[k].Strong()) {
      if (m_type>=3) {
	pt2k=pk.PPerp2();
	if (m_pt_def) pt2k+=pk.Abs2();
	if (j<3) {
	  if (pt2k<ymin*m_s) {
	    ymin=pt2k/m_s;
	    j1=j;
	    k1=k;
	  } 
	}
	else {
	  if (m_flavs[j].Strong()) {
	    pt2j=pj.PPerp2();
	    if (m_pt_def) pt2j+=pj.Abs2();
	    pt2jk=2.*Min(pt2j,pt2k)*
	      (Coshyp(DEta12(pj,pk))-CosDPhi12(pj,pk))/sqr(m_delta_r);
	    if (pt2jk<ymin*m_s) {
	      ymin=pt2jk/m_s;
	      j1=j;
	      k1=k;
	    }
	  }
	}
      }
      else {
	if (m_type==2) {
	  int hadron=m_fl[0].Strong()?0:1;
	  if (j==1<<hadron) {
	    pt2k=2.*sqr(pk[0])*(1.-DCos12(pk,p[hadron]));
	    if (pt2k<ymin*m_sprime) {
	      ymin=pt2k/m_sprime;
	      j1=j;
	      k1=k;
	    }
	  }
	}
	if (m_flavs[j].Strong()) {
	  pt2jk=2.*sqr(Min(pj[0],pk[0]))*(1.-DCos12(pj,pk));
	  if (pt2jk<ymin*m_sprime) {
	    ymin=pt2jk/m_sprime;
	    j1=j;
	    k1=k;
	  }
	}
      }
    }
  }
  return ymin;
}

/*----------------------------------------------------------------------------------

  Jet measure, mainly for showering and merging

  ----------------------------------------------------------------------------------*/

double Jet_Finder::MTij2(Vec4D p1,Vec4D p2)
{
  double mt12_2(0.);
  //check for DIS situation
  if (m_type>=2) {
    double pt1_2(p1.PPerp2()), mt1_2(pt1_2);
    double pt2_2(p2.PPerp2()), mt2_2(pt2_2);
    if (!m_pt_def) { 
      mt1_2 = pt1_2+dabs(p1.Abs2());
      mt2_2 = pt2_2+dabs(p2.Abs2());
    }
    if (IsZero(pt1_2/(pt1_2+pt2_2)) || IsZero(pt2_2/(pt1_2+pt2_2))) return mt1_2+mt2_2;
    else {
      if (m_type==2) mt12_2 = 2.*Min(mt1_2,mt2_2) * (1.-DCos12(p1,p2))/sqr(m_delta_r);
                else mt12_2 = 2.*Min(mt1_2,mt2_2) * 
		       (Coshyp(DEta12(p1,p2)) - CosDPhi12(p1,p2))/sqr(m_delta_r);
    }
  }
  else mt12_2               = 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2));
  return mt12_2;
}

bool Jet_Finder::TwoJets(const Vec4D & p1,const bool fix) 
{
  double crit = fix ? m_smin : m_shower_pt2;
  if (m_type>=2) {
    if (m_pt_def) { if(p1.PPerp2() < crit ) return false; }
             else { if(p1.MPerp2() < crit ) return false; }
  }
  else {
    msg.Out()<<"WARNING in Jet_Finder::TwoJets(Vec4D &) "<<std::endl
	     <<"    Still not implemented for mode "<<m_type<<std::endl;
  }
  return true;
}


/*----------------------------------------------------------------------------------

  Used in showers.

  ----------------------------------------------------------------------------------*/



bool Jet_Finder::TwoJets(const Vec4D & _p1,const Vec4D & _p2,const bool fix)
{
  Vec4D p1=_p1;
  Vec4D p2=_p2;

  BoostInFrame(p1);
  BoostInFrame(p2);

  double pt1_2(0.), pt2_2(0.), pt12_2(0.);
  double crit = fix ? m_smin : m_shower_pt2;

  if (m_type>=2) {
    if (m_pt_def) { pt1_2  = p1.PPerp2(); pt2_2  = p2.PPerp2(); }
             else { pt1_2  = p1.MPerp2(); pt2_2  = p2.MPerp2(); }
    if (pt1_2  < crit || pt2_2  < crit ) return 0;
    
    if (m_type==2) pt12_2 = 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2))/sqr(m_delta_r);
              else pt12_2 = 2.*Min(pt1_2,pt2_2) * (Coshyp(DEta12(p1,p2)) - CosDPhi12(p1,p2))/sqr(m_delta_r);
  }
  else {
    pt12_2 = 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2));

  //   std::cout<<" --- Jetfinder ee mode : "<<pt12_2<<" "<<m_shower_pt2<<" --- ";
  }
  if (pt12_2 < crit ) return 0;
  return 1;
}


bool Jet_Finder::TwoJets(double & E2,double & z,double & costheta,bool mode)
{
  double pt12_2(0.);
  if (mode == 1) {
    pt12_2 = -1000.;
    msg.Out()<<"WARNING in Jet_Finder::TwoJets(Vec4D &) "<<std::endl
	     <<"    Still not implemented for mode "<<m_type<<std::endl;
  }
  else pt12_2 = 2.*E2*sqr(Min(z,1.- z))*(1.-costheta)/sqr(m_delta_r);

  if (pt12_2 < m_shower_pt2 ) return 0;
  return 1;
}

/*----------------------------------------------------------------------------------

  Utilities

  ----------------------------------------------------------------------------------*/
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
  for (int i=0;i<m_n;i++) {
    m_cms_boost.Boost(p[i]);
    m_zrot.Rotate(p[i]);
  }
}

void Jet_Finder::BoostBack(Vec4D * p)
{
  for (int i=0;i<m_n;i++) {
    m_zrot.RotateBack(p[i]);
    m_cms_boost.BoostBack(p[i]);
  }
}

void Jet_Finder::BoostInFrame(std::vector<Vec4D> p)
{
  for (size_t i=0;i<p.size();i++) {
    m_cms_boost.Boost(p[i]);
    m_zrot.Rotate(p[i]);
  }
}

void Jet_Finder::BoostBack(std::vector<Vec4D> p)
{
  for (size_t i=0;i<p.size();i++) {
    m_zrot.RotateBack(p[i]);
    m_cms_boost.BoostBack(p[i]);
  }
}

void Jet_Finder::BoostInFrame(Vec4D & p)
{
  m_cms_boost.Boost(p);
  m_zrot.Rotate(p);
}

void Jet_Finder::BoostBack(Vec4D & p)
{
  m_zrot.RotateBack(p);
  m_cms_boost.BoostBack(p);
}

void Jet_Finder::SetDeltaR(double dr) 
{ 
  if (dr<=1.e-6) {
    msg.Error()<<"ERROR in Jet_Finder::SetDeltaR("<<dr<<") : "<<std::endl
	       <<"   delta_R to small, will ignore it and leave it at delta_R ="<<m_delta_r<<"."<<std::endl;
    return;
  }
  m_delta_r = dr; 
}

double Jet_Finder::ActualValue() const 
{
  return m_value; 
}
