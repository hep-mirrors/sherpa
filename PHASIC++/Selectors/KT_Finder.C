#include "PHASIC++/Selectors/KT_Finder.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

KT_Finder::KT_Finder(const std::string &_ycut,const int _type) : 
  Selector_Base("KTFinder"), m_value(0.), p_frame(NULL)
{
  m_ycut       = 2.0;
  m_delta_r    = 1.;
  m_type       = _type;
  m_cuttag     = _ycut;
  m_ene        = rpa->gen.Ecms()/2.;
  m_sprime     = m_s = sqr(2.*m_ene); 
  m_smax       = m_s;

  m_sel_log    = new Selector_Log(m_name);
}

/*---------------------------------------------------------------------

  Special form - flavours etc are known, will operate on momenta only.

  --------------------------------------------------------------------- */


KT_Finder::KT_Finder(const int _n,Flavour * _fl,
		       const std::string &_ycut,const int _type) : 
  Selector_Base("KTFinder"), m_value(0.), p_frame(NULL)
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
      m_ene      = rpa->gen.Ecms()/2.;
      m_sprime   = m_s = sqr(2.*m_ene); 
      p_frame[0] = Vec4D(m_ene,0.,0., sqrt(sqr(m_ene)-sqr(m_fl[0].Mass())));
      p_frame[1] = Vec4D(m_ene,0.,0.,-sqrt(sqr(m_ene)-sqr(m_fl[1].Mass())));
      if (m_type==3) m_cms_boost = Poincare(p_frame[0]+p_frame[1]);
    }    
    else if (m_type==2) {
      m_ene      = rpa->gen.Ecms()/2.;
      m_sprime   = m_s = sqr(2.*m_ene);
    }
  }
  
  m_smax    = m_s;
  m_sel_log = new Selector_Log(m_name);
}

KT_Finder::~KT_Finder() {
  if (p_frame)   delete [] p_frame;
}


/*----------------------------------------------------------------------------------

  Constructing jets, mainly for phase space cuts.

  ----------------------------------------------------------------------------------*/


void KT_Finder::Init(const Vec4D * p)
{
  if (m_nin==2) {
    switch (m_type) {
    case 4 : return;
    case 3 : {
      msg_Error()<<"KT_Finder::Init : process-type "<<m_type
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
	msg_Error()<<" ERROR: KT_Finder::Init could not initialize Breit frame correctly (1) : "
		   <<dabs(q*pp)<<std::endl;
      
      bool check(true);
      for (int i=0;i<3;i++) 
	if (dabs((q[i]-store[i]))>1.e10) check = false; 
      if (!check) msg_Error()<<" ERROR: KT_Finder::Init could not initialize Breit frame correctly (2) : "
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

bool KT_Finder::ConstructJets(Particle_List * pl, double y_res,int number, bool final_only) 
{
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
    if (rpa->gen.Beam1().Strong()   || rpa->gen.Beam2().Strong() || 
	rpa->gen.Beam1().IsHadron() || rpa->gen.Beam2().IsHadron()) 
      m_type=4;
    flavsin[0]=rpa->gen.Beam1();
    flavsin[1]=rpa->gen.Beam1();
    
    if (m_type!=2) {
    momsin[0]=Vec4D::ZVEC;
    momsin[0]*=(rpa->gen.Ecms()*0.5);
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
    (*it)->SetFlav(Flavour(kf_jet));
    (*it)->SetMomentum(momsout[j]);
  }

  return true;
}

bool KT_Finder::ConstructJets(const Particle_List * parts,
			       const std::vector<int> & jets,std::vector<double> & lastys,bool final_only) 
{
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
    if (rpa->gen.Beam1().Strong() || rpa->gen.Beam2().Strong() || 
	rpa->gen.Beam1().IsHadron() || rpa->gen.Beam2().IsHadron()) 
      m_type=4;
    flavsin[0]=rpa->gen.Beam1();
    flavsin[1]=rpa->gen.Beam1();

    if (m_type!=2) {
      momsin[0]=Vec4D::ZVEC;
      momsin[0]*=(rpa->gen.Ecms()*0.5);
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

bool KT_Finder::ConstructJetSystem(Vec4D * momsin,Flavour * flavsin,std::vector<Vec4D> & momsout,
				    std::vector<int> jets,std::vector<double> & lastys) 
{
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

double KT_Finder::YminKt(Vec4D * momsin,Flavour * flavsin,std::vector<Vec4D> momsout,int & j1,int & k1)
{
  double ymin = 2.;
  j1=-3; k1=-3;
  for (size_t j=0;j<momsout.size();j++) {
    if (m_type>=3) {
      double pt2j=Min(MTij2(momsout[j],momsin[0]),
		      MTij2(momsout[j],momsin[1]));
      if (pt2j < ymin*m_s) {
	ymin = pt2j/m_s;
	k1   = j;
	if (momsout[j][3]*momsin[0][3] > 0.) j1 = -2;
	                                else j1 = -1;
      }
      for (size_t k=j+1;k<momsout.size();k++) {
	double pt2jk = MTij2(momsout[j],momsout[k]);
	if (pt2jk<ymin*m_s) {
	  ymin = pt2jk/m_s;
	  j1 = j;k1 = k;
	}
      }
    }
    else {
      if (m_type==2) {
	int hadron=m_fl[0].Strong()?0:1;
	double pt2j = 2.*sqr(momsout[j][0])*
	  (1.-DCos12(momsout[j],momsin[hadron]));
	if (pt2j < ymin*m_sprime) {
	  ymin = pt2j/m_sprime;
	  k1   = j;
	  //cluster to beam 
	  if (hadron) j1 = -1;
	         else j1 = -2;
	}
      }
      for (size_t k=j+1;k<momsout.size();k++) {
	double pt2jk  = 2.*Min(momsout[j].PSpat2(),momsout[k].PSpat2())*
	  (1.-DCos12(momsout[j],momsout[k]));
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

Flavour KT_Finder::GetFlavour(std::string fl)
{
  bool bar(false);
  if (fl=="j") return Flavour(kf_jet);
  if (fl=="Q") return Flavour(kf_quark);
  if (fl=="G") return Flavour(kf_gluon);
  if (fl=="P") return Flavour(kf_photon);
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
  if (fl=="Q") return Flavour(kf_quark); // why again?
  Flavour flav(s_kftable.KFFromIDName(fl));
  if (flav.Kfcode()==kf_none) 
    THROW(critical_error,"No flavour for '"+fl+"'.");
  if (bar) flav=flav.Bar();
  return flav;
}

size_t KT_Finder::FillCombinations(const std::string &name,
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
	  sp=i=j+3;
	  break;
	}
      }
    }
    else if (name[i]=='_') {
      if (name[i-1]!='_' && name[i+1]=='_') {
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
  interpreter.AddTag("E_CMS",ToString(rpa->gen.Ecms()));
  m_cycut=ToType<double>(interpreter.Interprete(ccut));
  m_gcycut=ToType<double>(interpreter.Interprete(cgcut));
  for (size_t i(0);i<pos.size();++i) {
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

void KT_Finder::FillCombinations()
{
  if (m_ycuts.empty()) {
    if (p_proc==NULL) THROW(fatal_error,"Process not set.");
    m_procname=p_proc->Process()->Name();
    if (m_procname.find("__QCD")!=std::string::npos)
      m_procname=m_procname.substr(0,m_procname.find("__QCD"));
    if (m_procname.find("__EW")!=std::string::npos)
      m_procname=m_procname.substr(0,m_procname.find("__EW"));
    m_moms.clear();
    m_flavs.clear();
    m_ycuts.clear();
    m_gycuts.clear();
    m_fills.resize(m_nin+m_nout+1);
    std::string name(m_procname.substr(m_procname.find('_')+1));
    name=name.substr(name.find("__")+2);
    size_t i(0);
    FillCombinations(name,m_cuttag,rpa->gen.Variable("Y_CUT"),i,m_nin+m_nout);
    if (msg_LevelIsDebugging()) {
      msg_Out()<<METHOD<<"(): Combinations for '"<<m_procname<<"' {\n";
      double s(sqr(rpa->gen.Ecms()));
      for (std::map<size_t,std::map<size_t,double> >::const_iterator
	     iit(m_ycuts.begin());iit!=m_ycuts.end();++iit) {
	size_t i(iit->first);
	if (iit->second.size()>1)
	  msg_Out()<<"  "<<ID(i)<<"["<<m_flavs[i]<<","
		   <<m_flavs[i].Strong()<<"] & {";
	for (std::map<size_t,double>::const_iterator
	       jit(iit->second.begin());jit!=iit->second.end();++jit) {
	  size_t j(jit->first);
	  if (i!=j) 
	    msg_Out()<<" "<<ID(j)<<"["<<m_flavs[j]<<","
		     <<m_flavs[j].Strong()<<",("<<sqrt(m_ycuts[i][j]*s)
		     <<","<<sqrt(m_gycuts[i][j]*s)<<")]";
	}
	if (iit->second.size()>1) msg_Out()<<" }\n";
      }
      msg_Out()<<"}\n";
      msg_Out()<<METHOD<<"(): Identified clusterings {\n";
      for (size_t i(0);i<m_fills.size();++i)
	for (size_t j(0);j<m_fills[i].size();++j)
	  msg_Out()<<"  ["<<ID(m_fills[i][j].first)<<","
		   <<ID(m_fills[i][j].second)<<"] ("<<i<<")\n";
      msg_Out()<<"}\n";
      msg_Out()<<METHOD<<"(): Momentum combination {\n";
      for (size_t i(0);i<m_mcomb.size();++i) {
	msg_Out()<<"  "<<ID(m_mcomb[i].back())<<" -> {";
	for (size_t j(0);j<m_mcomb[i].size()-1;++j) 
	  msg_Out()<<" "<<ID(m_mcomb[i][j]);
	msg_Out()<<" }\n";
      }
      msg_Out()<<"}\n";
    }
  }
}

void KT_Finder::PrepareMomList(const Vec4D* vec) {
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

bool KT_Finder::Trigger(const Vec4D_Vector &p)
{
  FillCombinations();
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

void KT_Finder::BuildCuts(Cut_Data * cuts) 
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
	  cuts->etmin[i] = Max(sqrt(cut * m_s-sqr(m_fl[i].SelMass())),cuts->etmin[i]);
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
	    double scut=sqr(m_fl[i].SelMass())+sqr(m_fl[j].SelMass())+GetYcut(1<<i,1<<j)*m_s;
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

double KT_Finder::GetYcut(const size_t& i,const size_t& j) const {
  std::map<size_t,std::map<size_t,double> >::const_iterator it=m_ycuts.find(i);
  if(it==m_ycuts.end()) return -1.0;
  std::map<size_t,double>::const_iterator jt=(it->second).find(j);
  if(jt==(it->second).end()) return -1.0;
  return jt->second;
}
double KT_Finder::GetScaledYcut(const size_t &i,const size_t &j) const {
  double ycut(GetYcut(i,j));
  return i<3||j<3?ycut:ycut/sqr(m_delta_r);
}
double KT_Finder::GetGlobalYcut(const size_t &i,const size_t &j) const {
  std::map<size_t,std::map<size_t,double> >::const_iterator it=m_gycuts.find(i);
  if(it==m_gycuts.end()) return -1.0;
  std::map<size_t,double>::const_iterator jt=(it->second).find(j);
  if(jt==(it->second).end()) return -1.0;
  return jt->second;
}
double KT_Finder::GetScaledGlobalYcut(const size_t &i,const size_t &j) const {
  double gycut(GetGlobalYcut(i,j));
  return i<3||j<3?gycut:gycut/sqr(m_delta_r);
}

double KT_Finder::YminKt(int & j1,int & k1,int cl)
{
  double ymin = 2.;
  Vec4D p0(m_moms[1]), p1(m_moms[2]);
  double m0(m_flavs[1].SelMass()), m1(m_flavs[2].SelMass());
  for (size_t ps(0);ps<m_fills[cl].size();++ps) {
    int j(m_fills[cl][ps].first), k(m_fills[cl][ps].second);
    Vec4D pj(m_moms[j]), pk(m_moms[k]);
    double ycut(m_ycuts[j][k]);
    double mj(m_flavs[j].SelMass()), mk(m_flavs[k].SelMass());
    msg_Debugging()<<"  "<<ID(j)<<"["<<m_flavs[j]<<","<<mj<<"] & "
		   <<ID(k)<<"["<<m_flavs[k]<<","<<mk<<"], qcut = "
		   <<sqrt(ycut*m_s)<<"/"<<sqrt(m_gycuts[j][k]*m_s);
    if (m_flavs[k].Strong()) {
      if (m_type>=3) {
	if (j<3) {
	  double pt2k=Min(MTij2(p0,pk,m0,mk),MTij2(p1,pk,m1,mk));
	  msg_Debugging()<<", is -> ptk = "<<sqrt(pt2k)<<" ("
			 <<(pt2k>=ycut*m_s)<<(pt2k<ycut*m_s?")\n":")");
	  if (pt2k<ycut*m_s) return -1.0;
	  if (pt2k<ymin*m_s) {
	    ymin=pt2k/m_s;
	    j1=j;
	    k1=k;
	  } 
	}
	else if (m_flavs[j].Strong()) {
	  // delta r is taken into account in m_ycuts !
	  double pt2jk=MTij2(pj,pk,mj,mk);
	  msg_Debugging()<<", fs -> ptjk = "<<sqrt(pt2jk)<<" ("
			 <<(pt2jk>=ycut*m_s)<<(pt2jk<ycut*m_s?")\n":")");
	  if (pt2jk<ycut*m_s) return -1.0;
	  if (pt2jk<ymin*sqr(m_delta_r)*m_s) {
	    ymin=pt2jk/sqr(m_delta_r)/m_s;
	    j1=j;
	    k1=k;
	  }
	}
      }
      else {
	if (m_type==2) {
	  int hadron=m_fl[0].Strong()?0:1;
	  if (j==1<<hadron) {
	    double pt2k=2.*sqr(pk[0])*(1.-DCos12(pk,pj));
	    if (pt2k<ycut*m_sprime) return -1.0;
	    if (pt2k<ymin*m_sprime) {
	      ymin=pt2k/m_sprime;
	      j1=j;
	      k1=k;
	    }
	  }
	}
	if (j>2 && m_flavs[j].Strong()) {
	  double pt2jk=MTij2(pj,pk,mj,mk);
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

double KT_Finder::MTij2(Vec4D p1,Vec4D p2,double m1,double m2)
{
  if (m_type>=2) {
    double pt1_2(p1.PPerp2()), pt2_2(p2.PPerp2());
    if (IsZero(pt1_2) && IsZero(pt2_2)) return 2.0*p1*p2;
    if (IsZero(pt1_2)) return CPerp2(p2);
    if (IsZero(pt2_2)) return CPerp2(p1);
    return 2.0*Min(CPerp2(p1),CPerp2(p2))*CDij(p1,p2)/sqr(m_delta_r);
  }
  return 2.0*Min(p1.PSpat2(),p2.PSpat2())*(1.0-DCos12(p1,p2));
}

double KT_Finder::DCos12(Vec4D & p1,Vec4D & p2)
{
  return Vec3D(p1)*Vec3D(p2)/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

double KT_Finder::CPerp2(Vec4D& p) 
{
  return dabs(p.Abs2())+p.PPerp2();
}

double KT_Finder::CDij(Vec4D& p1,Vec4D& p2) 
{
  return p1*p2/sqrt(CPerp2(p1)*CPerp2(p2));
}


void KT_Finder::BoostInFrame(Vec4D * p)
{
  for (int i=0;i<m_n;i++) {
    m_cms_boost.Boost(p[i]);
    m_zrot.Rotate(p[i]);
  }
}

void KT_Finder::BoostBack(Vec4D * p)
{
  for (int i=0;i<m_n;i++) {
    m_zrot.RotateBack(p[i]);
    m_cms_boost.BoostBack(p[i]);
  }
}

void KT_Finder::BoostInFrame(std::vector<Vec4D> p)
{
  for (size_t i=0;i<p.size();i++) {
    m_cms_boost.Boost(p[i]);
    m_zrot.Rotate(p[i]);
  }
}

void KT_Finder::BoostBack(std::vector<Vec4D> p)
{
  for (size_t i=0;i<p.size();i++) {
    m_zrot.RotateBack(p[i]);
    m_cms_boost.BoostBack(p[i]);
  }
}

void KT_Finder::BoostInFrame(Vec4D & p)
{
  m_cms_boost.Boost(p);
  m_zrot.Rotate(p);
}

void KT_Finder::BoostBack(Vec4D & p)
{
  m_zrot.RotateBack(p);
  m_cms_boost.BoostBack(p);
}

void KT_Finder::SetDeltaR(double dr) 
{ 
  if (dr<=1.e-6) {
    msg_Error()<<METHOD<<"(): \\delta_R to small, ignore and set to "
               <<m_delta_r<<"."<<std::endl;
    return;
  }
  m_delta_r=dr; 
}

namespace PHASIC{

DECLARE_ND_GETTER(KT_Finder_Getter,"JetFinder",Selector_Base,Selector_Key,true);

Selector_Base *KT_Finder_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<2) THROW(critical_error,"Invalid syntax");
  int type(0);
  if (key.p_proc->NIn()==2) {
    int instrong(0);
    for (size_t j=0;j<key.p_proc->NIn();j++) { 
      if (key.p_proc->Process()->
	  Flavours()[j].Strong()) instrong++; 
    }
    if (instrong==0) type = 1;
    if (instrong==1) type = 2;
    if (instrong==2) type = 4;
  }
  KT_Finder *jf(new KT_Finder(key.p_proc->NIn()+key.p_proc->NOut(),
			      (Flavour*)&key.p_proc->Process()->
			      Flavours().front(),key[0][0],type));
  jf->SetDeltaR(ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1])));
  jf->SetProcess(key.p_proc);
  return jf;
}

void KT_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"kt jet finder"; 
}

}
