#include "Combine_Table.H"
#include "Run_Parameter.H"

using namespace SHERPA;
using namespace AMEGIC;
using namespace ATOOLS;

int Combine_Table::all=0;

using std::endl;

// ============================================================
//    class Combine_Data
// ============================================================

Combine_Data::Combine_Data():i(0),j(0),pt2ij(0),down(0)
{
}

Combine_Data::Combine_Data(int _i, int _j, double _pt2ij, int _ngraph):
  i(_i),j(_j),pt2ij(_pt2ij),down(0) 
{
  graphs.push_back(_ngraph);
}


Combine_Data::Combine_Data(const Combine_Data & cd):
  i(cd.i),j(cd.j),pt2ij(cd.pt2ij),graphs(cd.graphs),down(cd.down) 
{
}

Combine_Data::~Combine_Data() 
{
  if (down!=0) delete down;
  down=0;
}

std::ostream& SHERPA::operator<< (std::ostream & s ,Combine_Data & cd) 
{
    s<<" "<<cd.i<<"&"<<cd.j<<"   "<<cd.pt2ij<<"    "<<cd.strong<<"    "<<std::flush;
    for (size_t k=0; k<cd.graphs.size(); ++k)
    s<<cd.graphs[k]<<","<<std::flush;
    s<<"     ";
    if (cd.down)
    s<<" "<<cd.down->no<<std::endl;
    else 
    s<<" #"<<std::endl;
    return s;
}

Combine_Table::Combine_Table(Jet_Finder * _jf,Vec4D * _moms, Combine_Table * _up,
			     int isrmode, int isrshoweron):
  up(_up),legs(0),gwin(0),m_isr1on(isrmode&1),m_isr2on((isrmode&2)/2),
  m_isrshoweron(isrshoweron),jf(_jf),moms(_moms)
{
  no=all++;
}

bool Combine_Table::Combinable(const Leg & a , const Leg & b, int & strong) const 
{
  //std::cout<<"Combinable("<<a->fl<<","<<b->fl<<") : ";
  strong = 0;
  // 1.) check if both points have common mother
  if ((a->prev == b->prev) && (a->prev != 0)) {
    strong=a->prev->fl.Strong();
    //std::cout<<" Yes! "<<std::endl;
    return 1;
  }

  // 2.) check if "a" is daughter of "b"
  if (a->prev == &b)   {
    if (&a==b->left) strong=b->right->fl.Strong();
    else if (&a==b->right) strong=b->left->fl.Strong();
    //std::cout<<" Yes! "<<std::endl;
    return 1;
  }
  // 3.) check if "b" is daughter of "a"
  if (b->prev == &a)  {
    if (&b==a->left) strong=a->right->fl.Strong();
    else if (&b==a->right) strong=a->left->fl.Strong();
    //std::cout<<" Yes! "<<std::endl;
    return 1;
  }

  //std::cout<<" No! "<<std::endl;
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
      if (i<2) moms[i] = _moms[i] - _moms[j];      
      else     moms[i] = _moms[i] + _moms[j];
    }
    else       moms[l] = _moms[l];
  }
  for (int l=j+1; l<=maxl; ++l) moms[l-1]=_moms[l];
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
  legs  = _legs;
  nlegs = _nlegs;
  nampl = _nampl;

  // determine possible combinations and corresponding y_ij  if nlegs>4
  if (nlegs>4) {
    for (int k=0;k<nampl;++k) {
      //std::cout<<" Graph "<<k<<std::endl;
      //std::cout<<"=============================="<<endl;
      //std::cout<<(&legs[k][0]);
    }
    //std::cout<<"=============================="<<endl;

    int start=0;
    // cluster initial state only if isrshower and isr_x is on. 
    if (!m_isrshoweron) start=2;
    //std::cout<<" isr : "<<m_isrshoweron<<","<<m_isr1on<<","<<m_isr2on<<" ."<<std::endl;
    for (int i=start; i<nlegs; ++i) {  
      if (!m_isr1on && i==0) i=1;
      if (!m_isr2on && i==1) i=2;
      for (int j=i+1; j<nlegs; ++j) {
	// never combine "0&1" !
	if (j==1) j=2;
	// check if leg i is combinable with leg j in any graph
	for (int k=0;k<nampl;++k) {
	  int strong = 0;
	  //std::cout<<" gr : "<<k<<"  (i,j) : ("<<i<<","<<j<<")"<<endl;
	  if (Combinable(legs[k][i],legs[k][j],strong)) {  
	    AddPossibility(i,j,k,strong); // insert graph k with combination i&j in table
	  } 
	}	 
      }
    }
  }
}


void Combine_Table::AddPossibility(int i, int j, int ngraph, int strong) 
{
  if (combinations.size()==0) {
    // add new "i&j" combination to empty table
    combinations.push_back(Combine_Data(i,j,0.,ngraph));
    combinations.back().strong=strong;
  }
  else if ((combinations.back().i==i)&&(combinations.back().j==j)) {
    // add graph only ("i&j" row exists already)
    combinations.back().graphs.push_back(ngraph);
    if (strong) combinations.back().strong=strong; 
    //    we should probably only change to strong status if 
    //    strong graph is also to be used in shower initialization !!!
  } 
  else {
    // add new "i&j" combination 
    combinations.push_back(Combine_Data(i,j,0.,ngraph));
    combinations.back().strong=strong;
  }
}


Combine_Table * Combine_Table::CalcJet(int nl,double _x1,double _x2, ATOOLS::Vec4D * _moms) 
{
  int prefer_ew_clustering = 0;

  if (up==0) {
    x1 = _x1;
    x2 = _x2;
  }
  Combine_Table * ct=0;
  CD_List & cl=combinations;
  if (cl.size()==0) return this;

  // change momenta to actual values    
  if (_moms!=0) {
    for (int l=0;l<nl;++l)
      moms[l]=_moms[l];
  }
  
  // boost in CMS frame and rotate to z-axis (store old moms)
  Vec4D * save_moms=0;
  save_moms = new Vec4D[nl];
  for (int i=0;i<nl;++i) save_moms[i]=moms[i];
  
  Poincare cms,zaxis;
  bool did_boost=0;
  // boost needed?
  if (!(Vec3D(moms[0])==Vec3D(-1.*moms[1]))) {
    cms   = Poincare(moms[0]+moms[1]);
    for (int i=0;i<nl;++i) cms.Boost(moms[i]);
    zaxis = Poincare(moms[0],Vec4D::ZVEC);
    for (int i=0;i<nl;++i) zaxis.Rotate(moms[i]);
    did_boost=1;
    //      for (int i=0;i<nl;++i) std::cout<<i<<" : "<<moms[i]<<" ("<<moms[i].Abs2()<<std::endl;
  }
  
  // calculate pt2ij and determine "best" combination
  double pt2max   = sqr(rpa.gen.Ecms());
  double pt2min   = pt2max;
  double ewpt2min = pt2max;
  CD_Iterator ewcit=cl.end();
  for (CD_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
    double pt2ij = cit->pt2ij = jf->MTij2(moms[cit->i], moms[cit->j]);
    if (cit->i < 2  && pt2ij < pt2min ) {
      // check if is combination has right direction:
      double d = moms[cit->i][3] * moms[cit->j][3];
      if (d<0. ) {
	pt2ij*=1.001;
	cit->pt2ij = pt2ij;
      } 
      
      // check if combined momenta is physical even in LAB frame
      double e = save_moms[cit->i][0] - save_moms[cit->j][0];
      if (e<0. ) {
	//	  std::cout<<"lv ("<<cit->i<<","<<cit->j<<") "<<pt2ij<<endl;
	cit->pt2ij = pt2ij = pt2max;
      }
      else {
	// check if combined momenta are physical in their CMS frame
	Vec4D s1 = save_moms[cit->i] - save_moms[cit->j];
	Vec4D s2 = save_moms[1-cit->i];
	Poincare test(s1+s2);
	test.Boost(s1);
	test.Boost(s2);
	if (s1[0]<0. || s2[0]<0.) {
	  //	    std::cout<<"cv ("<<cit->i<<","<<cit->j<<") "<<pt2ij<<endl; 
	  cit->pt2ij = pt2ij = pt2max;
	}
      }
    }
    
    if (cit->strong==0 && pt2ij>pt2max && ewpt2min==pt2max) {
      ewpt2min=pt2ij;
      ewcit=cit;
    }
    // make sure min search starts with the highest scale
    if (pt2ij>pt2max && pt2min==pt2max) {
      pt2max=pt2min;
      pt2min=pt2ij;
      cwin = cit;
    }
    
    if (pt2ij<pt2min) {
      //	std::cout<<"or: ("<<cit->i<<","<<cit->j<<") "<<pt2ij<<" vs. "<<pt2min<<endl;       
      pt2min = pt2ij;
      cwin = cit;
    }
    if (cit->strong==0 && pt2ij<ewpt2min) {
      ewpt2min = pt2ij;
      ewcit=cit;
    }
  }
    
  if (ewcit!=cl.end() && prefer_ew_clustering && cwin!=ewcit ) { 
    // More stuff to come here for more sophisitcated treatment of ew stuff. 
    cwin=ewcit;
  }
  
  // check if boosted  (restore saved moms)
  if (did_boost) {
    for (int i=0;i<nl;++i) moms[i]=save_moms[i];
  }
  delete [] save_moms;
  
  --nl;
  
  // if number of legs is still greater 4 Cluster once more
  // if number of legs equals 4, determine end situation
  if (nl>=4) {
    if (!cwin->down) {
      Leg ** alegs = new Leg*[cwin->graphs.size()];
      for (size_t k=0;k<cwin->graphs.size();++k) {
	alegs[k]   = CombineLegs(legs[cwin->graphs[k]],cwin->i,cwin->j,nl);
      }
      Vec4D * amoms;
      CombineMoms(moms,cwin->i,cwin->j,nl,amoms); // generate new momenta
      cwin->down=new Combine_Table(jf,amoms,this,m_isr1on+2*m_isr2on,m_isrshoweron);
      cwin->down->FillTable(alegs,nl,cwin->graphs.size());   // initialise Combine_Table
    } 
    else {
      cwin->down->CombineMoms(moms,cwin->i,cwin->j,nl);
    }
    // update x1 (before calling CalcJet again
    if (cwin->i<2) {
      double z=cwin->down->Sprime()/Sprime();
      if (cwin->i==0) {
	cwin->down->x1=x1*z;
	cwin->down->x2=x2;
      }
      else {
	cwin->down->x1=x1;
	cwin->down->x2=x2*z;
      }
    }
    else {
      cwin->down->x1=x1;
      cwin->down->x2=x2;
    }
    
    ct   = cwin->down->CalcJet(nl,_x1,_x2);
    gwin = cwin->down->gwin;
    gwin = cwin->graphs[gwin];   // translate back
    
  } 
  else {
    msg.Error()<<" ERROR:  nlegs < 4 !!!!!!!!!"<<std::endl;
    abort();
  }
  return ct;
}

Combine_Table::~Combine_Table()
{
  delete [] moms;
 
  for (int k=0;k<nampl;++k) 
    delete [] legs[k];
  delete [] legs;
}

std::ostream& SHERPA::operator<< (std::ostream& s ,Combine_Table * ct) 
{
  if (ct) {
    s<<std::endl<<" Combine_Table "<<ct->no<<" (up=";
    if (ct->up) s<<ct->up->no<<")"<<std::endl; else s<<"#)"<<std::endl;
    s<<" x1="<<ct->x1<<" x2="<<ct->x2<<endl;
    s<<" ==============="<<std::endl;
    s<<"moms="<<ct->moms<<std::endl;
    for (int l=0; l<ct->nlegs; ++l) 
      s<<" "<<l<<":"<<ct->moms[l]<<std::endl;
    s<<" ---------------"<<std::endl;
    CD_List & cl=ct->combinations;
    if (cl.size()>0) {
      s<<"     with "<<cl.size()<<" combinations"<<std::endl;
      for (CD_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
	s<<(*cit); 
      }
      
      for (CD_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
	if (cit->down) {
	  s<<cit->down<<std::endl;
	}
      }

      // test output
      for (int k=0; k<ct->nampl; ++k) {
	for (int l=0; l<ct->nlegs; ++l) 
	  s<<" "<<ct->legs[k][l]->fl<<"("<<ct->legs[k][l].ExtraAnti()<<")  ";
	s<<std::endl;
      }

    }
    else {
      for (int k=0; k<ct->nampl; ++k) {
	for (int l=0; l<ct->nlegs; ++l) 
	  s<<" "<<ct->legs[k][l]->fl<<"("<<ct->legs[k][l].ExtraAnti()<<")  ";
	s<<std::endl;
      }
    }
  } 
  else
    s<<"***empty Combine_Table***"<<std::endl;
  return s;
}

double Combine_Table::Sprime() 
{
  if (!moms) {
    return 0;
  }
  return (moms[0]+moms[1]).Abs2();
}
