#include "Combine_Table.H"
#include "Run_Parameter.H"
#include "Random.H"
#include <iomanip>

using namespace SHERPA;
using namespace AMEGIC;
using namespace ATOOLS;

int Combine_Table::s_all=0;

using std::endl;

// ============================================================
//    class Combine_Data
// ============================================================

Combine_Key::Combine_Key():i(0),j(0)
{
}

Combine_Key::Combine_Key(int _i, int _j, const Flavour & _flav) :
  i(_i),j(_j),flav(_flav)
{
}

bool SHERPA::operator<(const Combine_Key & a, const Combine_Key & b) {
  if (a.i < b.i) return true;
  if (a.i > b.i) return false;
  if (a.j < b.j) return true;
  if (a.j > b.j) return false;
  //  if (a.flav.Kfcode()==kf::none && b.flav.Kfcode()!=kf::none) return true;
  if (a.flav.Kfcode() > b.flav.Kfcode()) return true;
  return false;
}

// ============================================================

Combine_Data::Combine_Data():pt2ij(0),down(0)
{
}

Combine_Data::Combine_Data(double _pt2ij, int _ngraph):
  pt2ij(_pt2ij),sij(0.),prop(0.),coupling(0.),weight(0.),down(0) 
{
  graphs.push_back(_ngraph);
}


Combine_Data::Combine_Data(const Combine_Data & cd):
  pt2ij(cd.pt2ij),sij(cd.sij),prop(cd.prop),coupling(cd.coupling),
  weight(cd.weight),strong(cd.strong),graphs(cd.graphs),down(0) 
{
  if (cd.down) {
    msg.Out()<<"Combine_Data::Combine_Data(const Combine_Data &)\n"
	     <<" you are not supposed to use this constructor "<<std::endl;
  }
}

Combine_Data::~Combine_Data() 
{
  if (down!=0) delete down;
}

std::ostream& SHERPA::operator<< (std::ostream& s,const Combine_Key & ck)
{
  s<<" "<<ck.i<<"&"<<ck.j;
  if (ck.flav.Kfcode()!=kf::none) s<<"["<<std::setw(6)<<ck.flav<<"]";
  else s<<std::string(8,' ');
  return s;
}

std::ostream& SHERPA::operator<< (std::ostream& s,const Combine_Data & cd)
{
  s<<" "<<std::setw(8)<<cd.pt2ij<<" "<<std::setw(8)<<cd.sij<<" "
   <<std::setw(11)<<cd.weight<<" ";
  if (cd.prop!=0.) s<<" ("<<std::setw(11)<<cd.prop<<","<<std::setw(11)<<cd.coupling<<") ";
  else s<<std::string(27,' ');
  for (size_t k=0; k<cd.graphs.size(); ++k)
    s<<cd.graphs[k]<<","<<std::flush;
  s<<"     ";
  if (cd.down)
    s<<" "<<cd.down->m_no<<std::endl;
  else 
    s<<" #"<<std::endl;
  return s;
}

Combine_Table::Combine_Table(Jet_Finder * _jf,Vec4D * _moms, Combine_Table * _up,
			     int isrmode, int isrshoweron):
  p_up(_up),p_legs(0),m_gwin(0),m_isr1on(isrmode&1),m_isr2on((isrmode&2)/2),
  m_isrshoweron(isrshoweron),p_jf(_jf),p_moms(_moms)
{
  m_mode=0;
  m_no=s_all++;
}

bool Combine_Table::Combinable(const Leg & a , const Leg & b, 
			       Vertex_Info & vinfo) const
// 			       Flavour & fl, Complex *& cpl, 
// 			       Color_Function *& color) const 
{
  // 1.) check if both points have common mother
  if ((a->prev == b->prev) && (a->prev != 0)) {
    vinfo.fl=a->prev->fl;
    vinfo.cpl = a->prev->cpl;
    vinfo.color = a->prev->Color;
    vinfo.mode  = 1;
    //    std::cout<<"fls="<<vinfo.fl<<","<<a->fl<<","<<b->fl<<"\n";
    return 1;
  }

  // 2.) check if "a" is daughter of "b"
  if (a->prev == &b)   {
    if (&a==b->left) {
      vinfo.fl  = b->right->fl;
      vinfo.cpl = b->cpl;
      vinfo.color = b->Color;
      vinfo.mode  = 2;
    }
    else if (&a==b->right) {
      vinfo.fl  = b->left->fl;
      vinfo.cpl = b->left->cpl;
      vinfo.color = b->left->Color;
      vinfo.mode  = 3;
    }
    //    std::cout<<"fls="<<vinfo.fl<<","<<a->fl<<","<<b->fl<<"\n";
    return 1;
  }

  // 3.) check if "b" is daughter of "a"
  if (b->prev == &a)  {
    if (&b==a->left) {
      vinfo.fl  = a->right->fl;
      vinfo.cpl = a->cpl;
      vinfo.color = a->Color;
      vinfo.mode  = 4;
    }
    else if (&b==a->right) {
      vinfo.fl  = a->left->fl;
      vinfo.cpl = a->cpl;
      vinfo.color = a->Color;
      vinfo.mode  = 5;
    }
    //    std::cout<<"fls="<<vinfo.fl<<","<<a->fl<<","<<b->fl<<"\n";
    return 1;
  }
  // else legs not combinable
  return 0;
}

Leg Combine_Table::CombinedLeg(Leg * legs, int i, int j)
{
  Leg & a = legs[i];
  Leg & b = legs[j];

  Leg mo;

  if ( (a->prev == b->prev) && (a->prev != 0) ) {
    // combinable-type (1.)
    mo.SetPoint(a->prev);
  } 
  else if (&a == b->left) {
    // combinable-type (2.a)
    mo.SetPoint(b->right);
  } 
  else if (&a == b->right) {
    // combinable-type (2.b)
    mo.SetPoint(b->left);
  } 
  else  if (&b == a->left) {
    // combinable-type (3.a)
    mo.SetPoint(a->right);
  } 
  else  if (&b == a->right) {
    // combinable-type (3.b)
    mo.SetPoint(a->left);
  } 
  else {
    msg.Error()<<" ERROR: cannot combine legs!"<<std::endl;
  }

  // fix charge incase initial state has wrong
  int icharge;
  if (i<2)  icharge = a.ExtraAnti()*a->fl.IntCharge() - b.ExtraAnti()*b->fl.IntCharge();
  else      icharge = a->fl.IntCharge() + b->fl.IntCharge();

  if (icharge!=mo->fl.IntCharge()) {
    mo.SetAnti(-1);
  }    
  
  return mo;
}
  
Leg * Combine_Table::CombineLegs(Leg *legs, int i, int j, int _nlegs) 
{
  Leg * alegs = new Leg[_nlegs];
  // assume i < j 

  for (int l=0; l<j; ++l) {
    if (l==i) alegs[i] = CombinedLeg(legs,i,j);
    else      alegs[l] = legs[l];
  }
  for (int l=j+1; l<=_nlegs; ++l) alegs[l-1] = legs[l];
  return alegs;
}


void Combine_Table::CombineMoms(Vec4D* _moms , int i, int j, int maxl) 
{
  // assume i < j
  for (int l=0; l<j; ++l) {
    if (l==i) { 
      if (i<2) p_moms[i] = _moms[i] - _moms[j];      
      else     p_moms[i] = _moms[i] + _moms[j];
    }
    else       p_moms[l] = _moms[l];
  }
  for (int l=j+1; l<=maxl; ++l) p_moms[l-1]=_moms[l];
}

void Combine_Table::CombineMoms(Vec4D * _moms ,int i,int j,int maxl,Vec4D *& omoms) 
{
  omoms = new Vec4D[maxl];
  // assume i < j
  for (int l=0; l<j; ++l) {
    if (l==i) {
      if (i<2) omoms[i] = _moms[i]-_moms[j];      
      else     omoms[i] = _moms[i]+_moms[j];
    }
    else       omoms[l] = _moms[l];
  }
  for (int l=j+1; l<=maxl; ++l) omoms[l-1]=_moms[l];
}

void Combine_Table::FillTable(Leg **_legs,int _nlegs, int _nampl)
{
  // store information
  p_legs  = _legs;
  m_nlegs = _nlegs;
  m_nampl = _nampl;

  // determine possible combinations and corresponding y_ij  if nlegs>4
  if (m_nlegs>4) {
    if (msg.LevelIsDebugging()) {
      for (int k=0;k<m_nampl;++k) {
	msg.Out()<<"Combine_Table for Graph "<<k<<std::endl
		 <<"=============================="<<endl
		 <<(&p_legs[k][0]);
      }
      msg.Out()<<"=============================="<<endl;
    }
    int start=0;
    // cluster initial state only if isrshower and isr_x is on. 
    if (!m_isrshoweron) start=2;
    for (int i=start; i<m_nlegs; ++i) {  
      if (!m_isr1on && i==0) i=1;
      if (!m_isr2on && i==1) i=2;
      for (int j=i+1; j<m_nlegs; ++j) {
	// never combine "0&1" !
	if (j==1) j=2;
	// check if leg i is combinable with leg j in any graph
	for (int k=0;k<m_nampl;++k) {
	  Vertex_Info vinfo;
	  if (Combinable(p_legs[k][i],p_legs[k][j],vinfo)) {  
	    AddPossibility(i,j,k,vinfo); // insert graph k with combination i&j in table
	  } 
	}	 
      }
    }
  }
}

double Combine_Table::ColorFactor(int i, int j, AMEGIC::Color_Function * const color, 
				  unsigned int mode)
{
  if (!((mode==1) || (mode==4 && i<2))) {
    msg.Out()<<" Combine_Table::ColorFactor("<<i<<","<<j<<") "
	     <<mode<<"-mode not covered yet\n";
    msg.Out()<<color->String()<<std::endl;
    //    return 1.;
  }

  // assume combinable type 1.) and final state clustering !!!
  if (color->Type()==cf::None)
    return 1.;
  if (color->Type()==cf::F) {
    if (i<2) return 2.*3.;
    return 3.;
  }
  if (color->Type()==cf::D) {
    if (color->StringArg(0)=='0' || color->StringArg(1)=='0')
      return 1.;
    return 3.;
  }
  if (color->Type()==cf::G) {
    if (color->StringArg(0)=='0' || color->StringArg(1)=='0')
      return 1.;
    return 8.;
  }
  if (color->Type()==cf::T) {
    if (color->StringArg(0)=='0')
      return 1./2.;
    return 4./3.;
  }
  return 1.;
}


void Combine_Table::AddPossibility(int i, int j, int ngraph, const Vertex_Info & vinfo) 
{
  CD_Iterator cit=m_combinations.find(Combine_Key(i,j));
  if (cit!=m_combinations.end()) {
    // add graph only ("i&j" row exists already)
    cit->second.graphs.push_back(ngraph);
    if (vinfo.fl.Strong()) cit->second.strong=true; 
    //    we should probably only change to strong status if 
    //    strong graph is also to be used in shower initialization !!!
    if (m_mode==1) {
      cit=m_combinations.find(Combine_Key(i,j,vinfo.fl));
      if (cit!=m_combinations.end()) {
	cit->second.graphs.push_back(ngraph);
      }
      else {
	double cfac = ColorFactor(i,j,vinfo.color,vinfo.mode);
	//	std::cout<<vinfo.color->String()<<"="<<cfac<<std::endl;
	Combine_Data cd(0.,ngraph);
	cd.strong   = vinfo.fl.Strong();
	cd.coupling = cfac*(norm(vinfo.cpl[0])+norm(vinfo.cpl[1]));
	m_combinations[Combine_Key(i,j,vinfo.fl)]=cd;
      }
    }
  }
  else {
    // add new "i&j" combination 
    //    std::cout<<vinfo.color->String()<<"="<<cfac<<std::endl;
    Combine_Data cd(0.,ngraph);
    cd.strong=vinfo.fl.Strong();
    m_combinations[Combine_Key(i,j)]=cd;
    if (m_mode==1) {
      double cfac = ColorFactor(i,j,vinfo.color,vinfo.mode);
      cd.coupling = cfac*(norm(vinfo.cpl[0])+norm(vinfo.cpl[1]));
      m_combinations[Combine_Key(i,j,vinfo.fl)]=cd;
    }
  }
}


Combine_Table * Combine_Table::CalcJet(int nl,double _x1,double _x2, ATOOLS::Vec4D * _moms) 
{
  int prefer_ew_clustering = 0;

  if (p_up==0) {
    m_x1 = _x1;
    m_x2 = _x2;
  }
  Combine_Table * ct=0;
  CD_List & cl=m_combinations;
  if (cl.size()==0) return this;

  // change momenta to actual values    
  if (_moms!=0) {
    for (int l=0;l<nl;++l)
      p_moms[l]=_moms[l];
  }
  
  // boost in CMS frame and rotate to z-axis (store old moms)
  Vec4D * save_moms=0;
  save_moms = new Vec4D[nl];
  for (int i=0;i<nl;++i) save_moms[i]=p_moms[i];
  
  Poincare cms,zaxis;
  bool did_boost=0;
  // boost needed?
  if (!(Vec3D(p_moms[0])==Vec3D(-1.*p_moms[1]))) {
    cms   = Poincare(p_moms[0]+p_moms[1]);
    for (int i=0;i<nl;++i) cms.Boost(p_moms[i]);
    zaxis = Poincare(p_moms[0],Vec4D::ZVEC);
    for (int i=0;i<nl;++i) zaxis.Rotate(p_moms[i]);
    did_boost=1;
  }
  
  // calculate pt2ij and determine "best" combination
  double pt2max   = sqr(rpa.gen.Ecms());
  double pt2min   = pt2max;
  double pt2min2  = 4.*pt2max;
  double prop_max = 0.;
  double ewpt2min = pt2max;

  CD_Iterator ewcit=cl.end();
  CD_Iterator pwcit=cl.end();
  CD_Iterator cwin2=cl.end();
  m_cwin=cl.end();
  for (CD_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    CD_Iterator tit = CalcPropagator(cit);
    // if flav==none calc  s, pt2ij
    // else determine coupling, propagator
    double pt2ij = cit->second.pt2ij;
    double prop  = tit->second.weight;

    if (prop>prop_max) {
      pwcit=tit;
      prop_max=prop;
    }

    if (tit==cit) {
      if (cit->first.i < 2  && pt2ij < pt2min ) {
	// check if is combination has right direction:
	double d = p_moms[cit->first.i][3] * p_moms[cit->first.j][3];
	if (d<0. ) {
	  pt2ij*=1.001;
	  cit->second.pt2ij = pt2ij;
	} 
	if (pt2ij<pt2min2) {
	  pt2min2 = pt2ij;
	  cwin2   = cit;
	}
      
	// check if combined momenta is physical even in LAB frame
	double e = save_moms[cit->first.i][0] - save_moms[cit->first.j][0];
	if (e<0. ) {
	  cit->second.pt2ij = pt2ij = pt2max;
	}
	else {
	  // check if combined momenta are physical in their CMS frame
	  Vec4D s1 = save_moms[cit->first.i] - save_moms[cit->first.j];
	  Vec4D s2 = save_moms[1-cit->first.i];
	  Poincare test(s1+s2);
	  test.Boost(s1);
	  test.Boost(s2);
	  if (s1[0]<0. || s2[0]<0.) {
	    cit->second.pt2ij = pt2ij = pt2max;
	  }
	}
      }
    }
    
    if (cit->second.strong==0 && pt2ij>pt2max && ewpt2min==pt2max) {
      ewpt2min=pt2ij;
      ewcit=cit;
    }
    // make sure min search starts with the highest scale
    if (pt2ij>pt2max && pt2min==pt2max) {
      pt2max=pt2min;
      pt2min=pt2ij;
      m_cwin = cit;
    }
    
    if (pt2ij<pt2min) {
      pt2min = pt2ij;
      m_cwin = cit;
    }
    if (pt2ij<pt2min2) {
      pt2min2 = pt2ij;
      cwin2   = cit;
    }

    if (cit->second.strong==0 && pt2ij<ewpt2min) {
      ewpt2min = pt2ij;
      ewcit=cit;
    }
  }
  if (m_cwin==cl.end()) 
    m_cwin=cwin2;


  if (m_mode==1) {
    double disc = ran.Get()*pwcit->second.weight;
    double sum = 0.;
    do {
      ++pwcit;
      sum+=pwcit->second.weight;
    } while (disc>sum);
    m_cwin=pwcit;
  }


    
  if (ewcit!=cl.end() && prefer_ew_clustering && m_cwin!=ewcit ) { 
    // More stuff to come here for more sophisticated treatment of ew stuff. 
    m_cwin=ewcit;
  }
  
  // check if boosted  (restore saved moms)
  if (did_boost) {
    for (int i=0;i<nl;++i) p_moms[i]=save_moms[i];
  }
  delete [] save_moms;
  
  --nl;
  
  // if number of legs is still greater 4 Cluster once more
  // if number of legs equals 4, determine end situation
  if (nl>=4) {
    if (!m_cwin->second.down) {
      Leg ** alegs = new Leg*[m_cwin->second.graphs.size()];
      for (size_t k=0;k<m_cwin->second.graphs.size();++k) {
	alegs[k]   = CombineLegs(p_legs[m_cwin->second.graphs[k]],m_cwin->first.i,m_cwin->first.j,nl);
      }
      Vec4D * amoms;
      CombineMoms(p_moms,m_cwin->first.i,m_cwin->first.j,nl,amoms); // generate new momenta
      m_cwin->second.down=new Combine_Table(p_jf,amoms,this,m_isr1on+2*m_isr2on,m_isrshoweron);
      m_cwin->second.down->FillTable(alegs,nl,m_cwin->second.graphs.size());   // initialise Combine_Table
    } 
    else {
      m_cwin->second.down->CombineMoms(p_moms,m_cwin->first.i,m_cwin->first.j,nl);
    }
    // update x1 (before calling CalcJet again
    if (m_cwin->first.i<2) {
      double z=m_cwin->second.down->Sprime()/Sprime();
      if (m_cwin->first.i==0) {
	m_cwin->second.down->m_x1=m_x1*z;
	m_cwin->second.down->m_x2=m_x2;
      }
      else {
	m_cwin->second.down->m_x1=m_x1;
	m_cwin->second.down->m_x2=m_x2*z;
      }
    }
    else {
      m_cwin->second.down->m_x1=m_x1;
      m_cwin->second.down->m_x2=m_x2;
    }
    
    ct   = m_cwin->second.down->CalcJet(nl,_x1,_x2);
    m_gwin = m_cwin->second.down->m_gwin;
    m_gwin = m_cwin->second.graphs[m_gwin];   // translate back
    
  } 
  else {
    msg.Error()<<"ERROR in Combine_Table::CalcJet :  nlegs < 4. Abort."<<std::endl;
    abort();
  }
  return ct;
}


CD_Iterator  Combine_Table::CalcPropagator(CD_Iterator & cit)
{
  if (cit->first.flav.Kfcode()==kf::none) {
    // determine general features of the i&j combination first 
    if (cit->first.i<2)
      cit->second.sij = (p_moms[cit->first.i]-p_moms[cit->first.j]).Abs2();
    else
      cit->second.sij = (p_moms[cit->first.i]+p_moms[cit->first.j]).Abs2();
    cit->second.pt2ij = p_jf->MTij2(p_moms[cit->first.i], p_moms[cit->first.j]);
    return cit;
  }
  else {
    // retrieve information from i&j master
    CD_Iterator father = m_combinations.find(Combine_Key(cit->first.i,cit->first.j));
    if (father!=m_combinations.end()) {
      // calculate propagator and sum them
      cit->second.pt2ij = father->second.pt2ij;
      double sij    = cit->second.sij = father->second.sij;
      double mass2  = sqr(cit->first.flav.Mass());
      double width2 = sqr(cit->first.flav.Width());
      if (mass2==0.) {
	cit->second.prop = 1./sqr(sij);
      }
      else {
	cit->second.prop = 1./(sqr(sij-mass2) + mass2*width2);
      }
      cit->second.weight = cit->second.prop*cit->second.coupling;
      father->second.weight += cit->second.weight;
      return father;
    }
    else {
      msg.Out()<<"WARNING: in Combine_Table: "<<cit->first<<std::endl;
    }
  }
  return cit;
}


Combine_Table::~Combine_Table()
{
  delete [] p_moms;
 
  for (int k=0;k<m_nampl;++k) 
    delete [] p_legs[k];
  delete [] p_legs;
}

std::ostream& SHERPA::operator<< (std::ostream& s ,const Combine_Table & ct) 
{
  if (&ct) {
    s<<std::endl<<" Combine_Table "<<ct.m_no<<" (up=";
    if (ct.p_up) s<<ct.p_up->m_no<<")"<<std::endl; else s<<"#)"<<std::endl;
    s<<" x1="<<ct.m_x1<<" x2="<<ct.m_x2<<endl;
    s<<" ==============="<<std::endl;
    s<<"moms="<<ct.p_moms<<std::endl;
    for (int l=0; l<ct.m_nlegs; ++l) 
      s<<" "<<l<<":"<<ct.p_moms[l]<<std::endl;
    s<<" ---------------"<<std::endl;
    const CD_List & cl=ct.m_combinations;
    if (cl.size()>0) {
      s<<"     with "<<cl.size()<<" combinations"<<std::endl;
      for (CD_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
 	s<<cit->first;
	s<<cit->second; 
      }
      for (CD_Const_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
	if (cit->second.down) {
	  s<<*cit->second.down<<std::endl;
	}
      }
      // test output
      for (int k=0; k<ct.m_nampl; ++k) {
	for (int l=0; l<ct.m_nlegs; ++l) 
	  s<<" "<<ct.p_legs[k][l]->fl<<"("<<ct.p_legs[k][l].ExtraAnti()<<")  ";
	s<<std::endl;
      }

    }
    else {
      for (int k=0; k<ct.m_nampl; ++k) {
	for (int l=0; l<ct.m_nlegs; ++l) 
	  s<<" "<<ct.p_legs[k][l]->fl<<"("<<ct.p_legs[k][l].ExtraAnti()<<")  ";
	s<<std::endl;
      }
    }
  } 
  else
    s<<"***empty Combine_Table***"<<std::endl;
  return s;
}

double Combine_Table::Sprime() const
{
  if (!p_moms) {
    return 0;
  }
  return (p_moms[0]+p_moms[1]).Abs2();
}
