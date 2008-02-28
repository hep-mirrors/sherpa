#include "Jet_Finder.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "Blob.H"
#include "Data_Reader.H"
#include <cassert>

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

Jet_Finder::Jet_Finder(const std::string &_ycut,const int _type) : 
  m_mass_scheme(1), m_value(0.), p_frame(NULL)
{
  m_ycut       = 2.0;
  m_delta_r    = 1.;
  m_type       = _type;
  m_cuttag     = _ycut;
  m_name       = std::string("Jetfinder");
  m_ene        = rpa.gen.Ecms()/2.;
  m_sprime     = m_s = sqr(2.*m_ene); 
  m_smax       = m_s;

  m_sel_log    = new Selector_Log(m_name);

  int help;
  Data_Reader reader(" ",";","!","=");
  if (reader.ReadFromFile(help,"JF_MASS_SCHEME")) {
    m_mass_scheme = help;
    msg_Tracking()<<METHOD<<": Set mass scheme to "<<m_mass_scheme<<"."<<std::endl;
  }
  //PRINT_INFO("  Jet finder mass scheme set to: "<<m_mass_scheme);
}

/*---------------------------------------------------------------------

  Special form - flavours etc are known, will operate on momenta only.

  --------------------------------------------------------------------- */


Jet_Finder::Jet_Finder(const int _n,Flavour * _fl,
		       const std::string &_ycut,const int _type) : 
  m_mass_scheme(1), m_value(0.), p_frame(NULL)
{
  m_ycut    = 2.0;
  m_delta_r = 1.;
  m_type    = _type;
  m_cuttag  = _ycut;

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
  
  m_smax    = m_s;
  m_sel_log = new Selector_Log(m_name);

  int help;
  Data_Reader reader(" ",";","!","=");
  if (reader.ReadFromFile(help,"JF_MASS_SCHEME")) {
    m_mass_scheme = help;
    msg_Tracking()<<METHOD<<": Set mass scheme to "<<m_mass_scheme<<"."<<std::endl;
  }
  //PRINT_INFO("  Jet finder mass scheme set to: "<<m_mass_scheme);
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
      msg_Error()<<"Jet_Finder::Init : process-type "<<m_type
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
	msg_Error()<<" ERROR: Jet_Finder::Init could not initialize Breit frame correctly (1) : "
		   <<dabs(q*pp)<<std::endl;
      
      bool check(true);
      for (int i=0;i<3;i++) 
	if (dabs((q[i]-store[i]))>1.e10) check = false; 
      if (!check) msg_Error()<<" ERROR: Jet_Finder::Init could not initialize Breit frame correctly (2) : "
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
      msg_Error()<<"This process-type is unknown!"<<std::endl;
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
	pt2jk  = 2.*Min(momsout[j].PSpat2(),momsout[k].PSpat2())*(1.-DCos12(momsout[j],momsout[k]));
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

Flavour Jet_Finder::GetFlavour(std::string fl)
{
  bool bar(false);
  if (fl=="j") return Flavour(kf::jet);
  if (fl=="Q") return Flavour(kf::quark);
  if (fl=="G") return Flavour(kf::gluon);
  if (fl=="P") return Flavour(kf::photon);
  if (fl.length()>1) {
    if (fl[fl.length()-1]=='b') {
      fl.erase(fl.length()-1,1);
      bar=true;
    }
    else if ((fl[0]=='W' || fl[0]=='H')) {
      if (fl[fl.length()-1]=='-') {
	fl[fl.length()-1]='+';
	bar=true;
      }
    }
    else if (fl[fl.length()-1]=='+') {
      fl[fl.length()-1]='-';
      bar=true;
    }
  }
  Flavour flav(kf_table.FromString(fl));
  if (fl=="Q") return Flavour(kf::quark);
  if (flav.Kfcode()==kf::none) 
    THROW(critical_error,"No flavour for '"+fl+"'.");
  if (bar) flav=flav.Bar();
  return flav;
}

size_t Jet_Finder::FillCombinations(const std::string &name,
				    const std::string &ycut,
				    const std::string &gycut,
				    size_t &cp,const int fl)
{
  bool ex(false);
  size_t sum(0), sp(0);
  std::vector<int> pos;
  std::string cut(ycut), ccut(cut), ncut(cut);
  std::string gcut(gycut), cgcut(gcut), ngcut(gcut);
  for (size_t i(0);i<name.length();++i) {
    if (name[i]=='[') {
      int open(1);
      for (size_t j(i+1);j<name.length();++j) {
	if (name[j]=='[') ++open;
	if (name[j]==']') --open;
	if (open==0) {
	  for (size_t ci(0);ci<cut.length();++ci) {
	    if (cut[ci]=='[') {
	      int copen(1);
	      for (size_t cj(ci+1);cj<cut.length();++cj) {
		if (cut[cj]=='[') ++copen;
		if (cut[cj]==']') --copen;
		if (copen==0) {
		  if (ccut==ycut) ccut=cut.substr(0,ci);
		  ncut=cut.substr(ci+1,cj-ci-1);
		  cut=cut.substr(cj+1);
		}
	      }
	    }
	  }
	  for (size_t ci(0);ci<gcut.length();++ci) {
	    if (gcut[ci]=='[') {
	      int copen(1);
	      for (size_t cj(ci+1);cj<gcut.length();++cj) {
		if (gcut[cj]=='[') ++copen;
		if (gcut[cj]==']') --copen;
		if (copen==0) {
		  if (cgcut==gycut) cgcut=gcut.substr(0,ci);
		  ngcut=gcut.substr(ci+1,cj-ci-1);
		  gcut=gcut.substr(cj+1);
		}
	      }
	    }
	  }
	  pos.push_back(FillCombinations
			(name.substr(i+1,j-i-1),ncut,ngcut,cp,fl-1));
	  m_flavs[pos.back()]=GetFlavour(name.substr(sp,i-sp));
	  sum=sum|pos.back();
	  sp=i=j+2;
	  break;
	}
      }
    }
    else if (name[i]=='_') {
      if (name[i-1]!='_' && (i<2 || name.rfind("nu",i)!=i-2) && 
	  name[i+1]!='R' && name[i+1]!='L' &&
	  (name[i+1]<48 || name[i+1]>57)) {
	pos.push_back(1<<cp++);
	m_flavs[pos.back()]=GetFlavour(name.substr(sp,i-sp));
	sum=sum|pos.back();
	ex=true;
      }
      if (ex && name[i+1]!='_') {
	sp=i+1;
	ex=false;
      }
    }
  }
  if (name[name.length()-1]!=']') {
    pos.push_back(1<<cp++);
    m_flavs[pos.back()]=GetFlavour(name.substr(sp,name.length()-sp));
    sum=sum|pos.back();
  }
  Algebra_Interpreter interpreter;
  interpreter.AddTag("E_CMS",ToString(rpa.gen.Ecms()));
  m_cycut=ToType<double>(interpreter.Interprete(ccut));
  m_gcycut=ToType<double>(interpreter.Interprete(cgcut));
  for (size_t i(0);i<pos.size();++i) {
    //for(std::map<size_t,size_t>::const_iterator it=m_map.begin(); it!=m_map.end(); ++it) {
    //  if(pos[i]<0) THROW(fatal_error,"This code is negative!");
    //  if(it->first==(size_t)pos[i]) THROW(fatal_error,"This code has been already kept!");
    //}
    //size_t cde(m_map.size());
    //m_map[pos[i]]=cde;
    //PRINT_INFO("pos,"<<i<<": "<<pos[i]<<"  "<<m_map[pos[i]]<<"  "<<m_flavs[pos[i]]);
    bool isi((pos[i]&(1<<0)) || (pos[i]&(1<<1)));
    m_ycuts[pos[i]][pos[i]]=isi?m_cycut:m_cycut*sqr(m_delta_r);
    m_gycuts[pos[i]][pos[i]]=isi?m_gcycut:m_gcycut*sqr(m_delta_r);
    for (size_t j(i+1);j<pos.size();++j) {
      if (pos[i]>2 || pos[j]>2) {
	m_ycuts[pos[i]][pos[j]]=isi?m_cycut:m_cycut*sqr(m_delta_r);
	m_gycuts[pos[i]][pos[j]]=isi?m_gcycut:m_gcycut*sqr(m_delta_r);
	m_ycut=Min(m_ycut,m_cycut);
	m_gycut=Min(m_gycut,m_gcycut);
	m_fills[fl].push_back(std::pair<size_t,size_t>(pos[i],pos[j]));
      }
    }
  }
  m_mcomb.push_back(pos);
  m_mcomb.back().push_back(sum);
  return sum;
}

void Jet_Finder::FillCombinations()
{
  if (m_ycuts.empty()) {
    if (m_procname=="") THROW(fatal_error,"Process name not set.");
    //size_t size(1<<(m_nin+m_nout)); PRINT_INFO(m_nin<<","<<m_nout<<","<<size<<","<<(1<<0));
    m_moms.clear();
    m_flavs.clear();
    m_ycuts.clear();
    m_gycuts.clear();
    m_fills.resize(m_nin+m_nout+1);
    std::string name(m_procname.substr(m_procname.find('_')+1));
    name=name.substr(name.find('_')+1);
    size_t i(0);
    FillCombinations(name,m_cuttag,rpa.gen.Variable("Y_CUT"),i,m_nin+m_nout);
/*
for(std::map<size_t,Flavour>::const_iterator it=m_flavs.begin(); it!=m_flavs.end(); ++it) {
PRINT_INFO(it->first<<":"<<it->second);}
for(std::map<size_t,std::map<size_t,double> >::const_iterator it=m_ycuts.begin(); it!=m_ycuts.end(); ++it) {
for(std::map<size_t,double>::const_iterator jt=m_ycuts[it->first].begin(); jt!=m_ycuts[it->first].end(); ++jt) {
PRINT_INFO(it->first<<","<<jt->first<<":"<<jt->second<<"="<<m_ycuts[it->first][jt->first]);
}}
*/
    if (msg_LevelIsDebugging()) {
      msg_Out()<<METHOD<<"(): Combinations for '"<<m_procname<<"' {\n";
      double s(sqr(rpa.gen.Ecms()));
      for(std::map<size_t,std::map<size_t,double> >::const_iterator it=m_ycuts.begin(); it!=m_ycuts.end(); ++it) {
        if(it->second.size()>1) {
	  msg_Out()<<"  "<<it->first<<" : "<<ID(it->first)<<"["<<m_flavs[it->first]<<","<<m_flavs[it->first].Strong()<<"] & {";
          for(std::map<size_t,double>::const_iterator jt=(it->second).begin(); jt!=(it->second).end(); ++jt)
	    if(it->first!=jt->first)
	      msg_Out()<<" "<<ID(jt->first)<<"["<<m_flavs[jt->first]<<","
	               <<m_flavs[jt->first].Strong()<<",("<<sqrt(m_ycuts[it->first][jt->first]*s)
	               <<","<<sqrt(m_gycuts[it->first][jt->first]*s)<<")]";
	  msg_Out()<<" }\n";
        }
      }
      msg_Out()<<"}\n";
      msg_Out()<<METHOD<<"(): Identified clusterings {\n";
      for (size_t i(0);i<m_fills.size();++i)
	for (size_t j(0);j<m_fills[i].size();++j)
	  msg_Out()<<"  ["<<ID(m_fills[i][j].first)<<","<<ID(m_fills[i][j].second)<<"] ("<<i<<")\n";
      msg_Out()<<"}\n";
      msg_Out()<<METHOD<<"(): Momentum combination {\n";
      for (size_t i(0);i<m_mcomb.size();++i) {
	msg_Out()<<"  "<<ID(m_mcomb[i].back())<<" -> {";
	for (size_t j(0);j<m_mcomb[i].size()-1;++j) msg_Out()<<" "<<ID(m_mcomb[i][j]);
	msg_Out()<<" }\n";
      }
      msg_Out()<<"}\n";
    }
  }
}

void Jet_Finder::PrepareMomList(const Vec4D* vec) {
  for (int i(m_nin+m_nout-1);i>=0;--i) m_moms[1<<i]=vec[i];
  for (size_t n(0);n<m_mcomb.size()-1;++n) {
    m_moms[m_mcomb[n].back()]=m_moms[m_mcomb[n].front()];
    for (size_t i(1);i<m_mcomb[n].size()-1;++i)
      m_moms[m_mcomb[n].back()]+=m_moms[m_mcomb[n][i]];
#ifdef BOOST_Decays
    Poincare cms(m_moms[m_mcomb[n].back()]);
    for (size_t i(0);i<m_mcomb[n].size()-1;++i) {
      cms.Boost(m_moms[m_mcomb[n][i]]);
      cc+=m_moms[m_mcomb[n][i]];
    }
    static double accu(sqrt(Accu()));
    Vec4D::SetAccu(accu);
    if (!(Vec3D(cc)==Vec3D()) || 
	!IsEqual(cc.Abs2(),m_moms[m_mcomb[n].back()].Abs2())) 
      msg_Error()<<METHOD<<"(): CMS boost failure. sum = "
		 <<cc<<" "<<cc.Abs2()<<" vs. "
		 <<m_moms[m_mcomb[n].back()].Abs2()<<"\n";
    Vec4D::ResetAccu();
#endif
//     msg_Debugging()<<"p["<<ID(m_mcomb[n].back())<<"] = "
//  		   <<m_moms[m_mcomb[n].back()]
//  		   <<" ["<<m_flavs[m_mcomb[n].back()]<<"] "<<cc<<"\n";
  }
}

bool Jet_Finder::Trigger(const Vec4D * p)
{
  FillCombinations();
  PROFILE_HERE;
  // create copy
  std::vector<Vec4D> vec(m_nin+m_nout);
  for (int i(0);i<m_nin+m_nout;++i) {
    vec[i]=p[i];
    //msg_Debugging()<<"p["<<i<<"] = "<<vec[i]<<" ("<<m_flavs[1<<i]<<")\n";
  }

  Init(&vec.front());
  BoostInFrame(&vec.front());
  PrepareMomList(&vec.front());

  int    j,k;
  double ymin(2.0);
  msg_Debugging()<<METHOD<<"() {\n";
  for (short unsigned int cl(0);cl<m_fills.size();++cl) {
    ymin=Min(ymin,YminKt(j,k,cl));
    if (ymin<0.0) return 1-m_sel_log->Hit(true);
  }
  msg_Debugging()<<"} -> q_min = "<<sqrt(ymin*m_s)<<"\n";
  m_value = ymin;
  return 1-m_sel_log->Hit(false);
}

void Jet_Finder::BuildCuts(Cut_Data * cuts) 
{
  FillCombinations();
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = m_fl[i].SelMass();
    if (m_fl[i].Strong()) {                
      if (m_type==1) {
	cuts->energymin[i] = Max(sqrt(GetScaledYcut(1<<i,1<<i) * m_s/4.),
				 cuts->energymin[i]);
      }
      else {
	cuts->energymin[i] = Max(sqrt(GetScaledYcut(1<<i,1<<i) * m_s),
				 cuts->energymin[i]);
	double cut(GetYcut(1<<0,1<<i));
	if (GetYcut(1<<1,1<<i)>0.0) { 
	  if (cut>0.0) cut=Min(cut,GetYcut(1<<1,1<<i));
	  else cut=GetYcut(1<<1,1<<i);
	}
	if (m_type==4 && cut>0.0 ) { 
	  if (!m_fl[i].IsMassive()) {
	    cuts->cosmax[0][i] = cuts->cosmax[1][i] = 
	      cuts->cosmax[i][0] = cuts->cosmax[i][1] =  
	      Min(cuts->cosmax[0][i],sqrt(1.-4.*cut));
	  }
	  switch (m_mass_scheme) {
	  case 1:
	  case 10:
	  case 11:
	    cuts->etmin[i] = Max(sqrt(cut * m_s-sqr(m_fl[i].SelMass())),cuts->etmin[i]);
	    break;
	  default:
	    cuts->etmin[i] = Max(sqrt(cut * m_s),cuts->etmin[i]);	    
	  }
	}
	if (m_type==2) {
	  int hadron=m_fl[0].Strong()?0:1;
	  double cut(GetYcut(1<<hadron,1<<i));
	  if (cut>0.0) {
	    cuts->cosmax[hadron][i] = cuts->cosmax[i][hadron] = 
	      Min(cuts->cosmax[hadron][i],sqrt(1.-4.*cut));
	    cuts->cosmin[hadron][i] = cuts->cosmin[i][hadron] = 
	      Max(cuts->cosmin[hadron][i],-sqrt(1.-4.*cut));
	    cuts->etmin[i] = Max(sqrt(cut * m_s),cuts->etmin[i]);
	  }
	}
      }
      
      for (int j=i+1; j<m_nin+m_nout; ++j) {
	if (m_fl[j].Strong() && GetYcut(1<<i,1<<j)>0.0) {
          if (m_type>=2 && (m_fl[i].IsMassive() || m_fl[j].IsMassive())) {
	    double scut(0.);
	    switch (m_mass_scheme) {
	    case 1:
	      scut = Max(sqr(m_fl[i].SelMass()+m_fl[j].SelMass()),
			 GetYcut(1<<i,1<<j)*m_s-sqr(m_fl[i].SelMass())-sqr(m_fl[j].SelMass()));
	      break;
	    case 0:
	    case 3:
	    case 10:
	    case 20:
	      scut = sqr(m_fl[i].SelMass()+m_fl[j].SelMass())+GetYcut(1<<i,1<<j)*m_s;
	      break;
	    case 11:
	    case 21:
	      scut = sqr(m_fl[i].SelMass())+sqr(m_fl[j].SelMass())+GetYcut(1<<i,1<<j)*m_s;
	    }
	    cuts->scut[i][j] = cuts->scut[j][i] = Max(cuts->scut[i][j],scut);
	  }
          else {
	    cuts->scut[i][j] = cuts->scut[j][i]
	      = Max(cuts->scut[i][j],GetYcut(1<<i,1<<j)*m_s);
	  }
        }
      }

    }
  }
}

void   Jet_Finder::UpdateCuts(double sprime,double y,Cut_Data * cuts) 
{
  if (m_type>1) return;
  for (int i=m_nin; i<m_nin+m_nout; ++i) {
    if (m_fl[i].Strong()) {                
      if (m_type==1)
	cuts->energymin[i] = Max(sqrt(GetYcut(1<<i,1<<i) * m_s/4.),
				 cuts->energymin[i]);
      else
	cuts->energymin[i] = Max(sqrt(GetYcut(1<<i,1<<i) * m_s),
				 cuts->energymin[i]);
      for (int j=i+1; j<m_nin+m_nout; ++j) {
	if (m_fl[j].Strong() && GetYcut(1<<i,1<<j)>0.0) {
	  if (m_type>=2) cuts->scut[j][i] = cuts->scut[i][j] 
			   = Max(cuts->scut[i][j],GetYcut(1<<i,1<<j)*m_s);
	  else cuts->scut[i][j] = cuts->scut[j][i] = 
		 Max(cuts->scut[i][j],GetYcut(1<<i,1<<j)*m_s);
	}
      }
    }
  }
}

double Jet_Finder::GetYcut(const size_t& i,const size_t& j) const {
  std::map<size_t,std::map<size_t,double> >::const_iterator it=m_ycuts.find(i);
  if(it==m_ycuts.end()) return -1.0;
  std::map<size_t,double>::const_iterator jt=(it->second).find(j);
  if(jt==(it->second).end()) return -1.0;
  return jt->second;
}
double Jet_Finder::GetScaledYcut(const size_t &i,const size_t &j) const {
  double ycut(GetYcut(i,j));
  return i<3||j<3?ycut:ycut/sqr(m_delta_r);
}
double Jet_Finder::GetGlobalYcut(const size_t &i,const size_t &j) const {
  std::map<size_t,std::map<size_t,double> >::const_iterator it=m_gycuts.find(i);
  if(it==m_gycuts.end()) return -1.0;
  std::map<size_t,double>::const_iterator jt=(it->second).find(j);
  if(jt==(it->second).end()) return -1.0;
  return jt->second;
}
double Jet_Finder::GetScaledGlobalYcut(const size_t &i,const size_t &j) const {
  double gycut(GetGlobalYcut(i,j));
  return i<3||j<3?gycut:gycut/sqr(m_delta_r);
}

double Jet_Finder::YminKt(int & j1,int & k1,int cl)
{
  PROFILE_HERE;
  double ymin = 2.;
  double pt2jk,pt2j,pt2k;
  for (size_t ps(0);ps<m_fills[cl].size();++ps) {
    int j(m_fills[cl][ps].first), k(m_fills[cl][ps].second);
    Vec4D pj(m_moms[j]), pk(m_moms[k]);
    double ycut(m_ycuts[j][k]), mj(m_flavs[j].Mass()), mk(m_flavs[k].Mass());
    msg_Debugging()<<"  "<<ID(j)<<"["<<m_flavs[j]<<","<<mj<<"] & "
   		   <<ID(k)<<"["<<m_flavs[k]<<","<<mk<<"], qcut = "
		   <<sqrt(ycut*m_s)<<"/"<<sqrt(m_gycuts[j][k]*m_s);
    if (m_flavs[k].Strong()) {
      double add(0.0);
      if (m_mass_scheme==1) add+=dabs(pk.Abs2());
      else if (m_mass_scheme==3) add+=dabs(pk.Abs2()-mk*mk);
      if (m_type>=3) {
	pt2k=CPerp2(pk);
	if (j<3) {
	  msg_Debugging()<<", is -> ptk = "<<sqrt(pt2k)<<" ("
			 <<(pt2k>=ycut*m_s)<<(pt2k<ycut*m_s?")\n":")");
	  if (add+pt2k<ycut*m_s) return -1.0;
	  if (add+pt2k<ymin*m_s) {
	    ymin=(add+pt2k)/m_s;
	    j1=j;
	    k1=k;
	  } 
	}
	else if (m_flavs[j].Strong()) {
	  pt2j=CPerp2(pj);
	  if (m_mass_scheme==1) add+=dabs(pj.Abs2());
	  else if (m_mass_scheme==3) add+=dabs(pj.Abs2()-mj*mj);
	  // delta r is taken into account in m_ycuts !
	  pt2jk=2.*Min(pt2j,pt2k)*CDij(pj,pk);
	  msg_Debugging()<<", fs -> ptjk = "<<sqrt(pt2jk)<<" ("
			 <<(pt2jk>=ycut*m_s)<<(pt2jk<ycut*m_s?")\n":")");
	  if (add+pt2jk<ycut*m_s) return -1.0;
	  if (add+pt2jk<ymin*sqr(m_delta_r)*m_s) {
	    ymin=(add+pt2jk)/sqr(m_delta_r)/m_s;
	    j1=j;
	    k1=k;
	  }
	}
      }
      else {
	if (m_type==2) {
	  int hadron=m_fl[0].Strong()?0:1;
	  if (j==1<<hadron) {
	    pt2k=2.*sqr(pk[0])*(1.-DCos12(pk,pj));
	    if (pt2k<ycut*m_sprime) return -1.0;
	    if (pt2k<ymin*m_sprime) {
	      ymin=pt2k/m_sprime;
	      j1=j;
	      k1=k;
	    }
	  }
	}
	if (j>2 && m_flavs[j].Strong()) {
	  pt2jk=2.0*Min(pj.PSpat2(),pk.PSpat2())*(1.-DCos12(pj,pk));
	  msg_Debugging()<<", fs -> ptjk = "<<sqrt(pt2jk)<<" s' = "
			 <<m_sprime<<" ("<<(pt2jk>=ycut*m_sprime)
			 <<(pt2jk<ycut*m_sprime?")\n":")");
	  if (pt2jk<ycut*m_sprime) return -1.0;
	  if (pt2jk/sqr(m_delta_r)<ymin*m_sprime) {
	    ymin=pt2jk/sqr(m_delta_r)/m_sprime;
	    j1=j;
	    k1=k;
	  }
	}
      }
    }
    msg_Debugging()<<"\n";
  }
  return ymin;
}

/*----------------------------------------------------------------------------------

  Jet measure, mainly for showering and merging

  ----------------------------------------------------------------------------------*/

double Jet_Finder::MTij2(Vec4D p1,Vec4D p2,double m1,double m2)
{
  double mt12_2(0.);
  //check for DIS situation
  if (m_type>=2) {
    double pt1_2(CPerp2(p1)), pt2_2(CPerp2(p1)), add(0.0);
    if (m_mass_scheme==1) add+=dabs(p1.Abs2())+dabs(p2.Abs2());
    else if (m_mass_scheme==3) add+=dabs(p1.Abs2()-m1*m1)+dabs(p2.Abs2()-m2*m2);
    if (IsZero(pt1_2) && IsZero(pt2_2)) return Max((p1+p2).Abs2(),add);
    if (IsZero(pt1_2/(pt1_2+pt2_2)) || 
	IsZero(pt2_2/(pt1_2+pt2_2))) return add+pt1_2+pt2_2;
    else {
      if (m_type==2) mt12_2 = add+2.*Min(pt1_2,pt2_2)*
	(1.-DCos12(p1,p2))/sqr(m_delta_r);
      else mt12_2 = (add+2.*Min(pt1_2,pt2_2)*CDij(p1,p2))/sqr(m_delta_r);
    }
  }
  else {
    double add(0.0);
    if (m_mass_scheme==1) add+=dabs(p1.Abs2())+dabs(p2.Abs2());
    else if (m_mass_scheme==3)
      add+=dabs(p1.Abs2()-m1*m1)+dabs(p2.Abs2()-m2*m2);
    mt12_2 = 2.0*Min(p1.PSpat2(),p2.PSpat2())*(1.0-DCos12(p1,p2));
  }
  return mt12_2;
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

double Jet_Finder::CPerp2(Vec4D& p) 
{
  switch (m_mass_scheme) {
  case 10:
  case 11:
    return p.MPerp2();
  case 20:
  case 21:
    return p.EPerp2();
  default:
    return p.PPerp2();
  }
}

double Jet_Finder::CDij(Vec4D& p1,Vec4D& p2) 
{
  switch (m_mass_scheme) {
  case 10:
    return (Coshyp(DEta12(p1,p2))-CosDPhi12(p1,p2))*sqrt(p1.PPerp2()/p1.MPerp2()*p2.PPerp2()/p2.MPerp2());
  case 11:
    return p1*p2/sqrt(p1.MPerp2()*p2.MPerp2());
  case 20:
    return (Coshyp(DEta12(p1,p2))-CosDPhi12(p1,p2))*sqrt(p1.PPerp2()/p1.EPerp2()*p2.PPerp2()/p2.EPerp2());
  case 21:
    return p1*p2/sqrt(p1.EPerp2()*p2.EPerp2());
  default:
    return Coshyp(DEta12(p1,p2))-CosDPhi12(p1,p2);   
  }
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

double Jet_Finder::ActualValue() const 
{
  return m_value; 
}
