#include "Combine_Table.H"
#include "Run_Parameter.H"

using namespace MOCAIC;
using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;

int Combine_Table::all=0;


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

Combine_Data::~Combine_Data() {
  if (down==0) delete down;
}

std::ostream& MOCAIC::operator<< (std::ostream & s ,Combine_Data & cd) {
    s<<" "<<cd.i<<"&"<<cd.j<<"   "<<cd.pt2ij<<"    "<<std::flush;
    for (int k=0; k<cd.graphs.size(); ++k)
    s<<cd.graphs[k]<<","<<std::flush;
    s<<"     ";
    if (cd.down)
    s<<" "<<cd.down->no<<std::endl;
    else 
    s<<" #"<<std::endl;
}

Combine_Table::Combine_Table(Jet_Finder * _jf,vec4d * _moms, Combine_Table * _up):
  jf(_jf),moms(_moms),legs(0),gwin(0),up(_up)
{
  msg.Debugging()<<"creating new Combine_Table::Combine_Table() "<<std::endl;
  no=all++;
}

inline bool Combine_Table::Combinable(const Leg & a , const Leg & b) const {
  // 1.) check if both points have common mother
  if ((a->prev == b->prev) && (a->prev != 0))   return 1;

  // 2.) check if "a" is daughter of "b"
  if (a->prev == &b)                            return 1;

  // 3.) check if "b" is daughter of "a"
  if (b->prev == &a)                            return 1;

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
  if (i<2)  icharge = a.ExtraAnti()*a->fl.icharge() - b.ExtraAnti()*b->fl.icharge();
  else      icharge = a->fl.icharge() + b->fl.icharge();

  if (icharge!=mo->fl.icharge()) {
    std::cout<<" changing "<<mo->fl<<" to "<<Flavour(mo->fl).bar()<<", "<<std::endl;
    mo.SetAnti(-1);
  }    
  
  return mo;
}
  
Leg * Combine_Table::CombineLegs(Leg *legs, int i, int j, int _nlegs) {
  Leg * alegs = new Leg[_nlegs];
  // assume i < j 

  for (int l=0; l<j; ++l) {
    if (l==i) alegs[i] = CombinedLeg(legs,i,j);
    else      alegs[l] = legs[l];
  }
  for (int l=j+1; l<=_nlegs; ++l) alegs[l-1] = legs[l];
  return alegs;
}


void Combine_Table::CombineMoms(vec4d* _moms , int i, int j, int maxl) {
  
  msg.Debugging()<<"CombineMoms(i="<<i<<",  j="<<j<<",  maxl="<<maxl<<")"<<std::endl;
  // assume i < j
  for (int l=0; l<j; ++l) {
    if (l==i) { 
      if (i<2) moms[i] = _moms[i] - _moms[j];      
      else     moms[i] = _moms[i] + _moms[j];
    }
    else       moms[l] = _moms[l];
  }
  for (int l=j+1; l<=maxl; ++l) 
    moms[l-1]=_moms[l];
  if (rpa.gen.Debugging()) {
    for (int l=0; l<maxl;++l) 
      msg.Out()<<" moms["<<l<<"]="<<moms[l]<<std::endl;
  }
}

void Combine_Table::CombineMoms(vec4d * _moms ,int i,int j,int maxl,vec4d *& omoms) {
  omoms = new vec4d[maxl];
  msg.Debugging()<<"CombineMoms(i="<<i<<",  j="<<j<<",  maxl="<<maxl<<",  omoms="<<omoms<<")"<<std::endl;
  
  // assume i < j
  for (int l=0; l<j; ++l) {
    if (l==i) {
      if (i<2) omoms[i] = _moms[i]-_moms[j];      
      else     omoms[i] = _moms[i]+_moms[j];
    }
    else       omoms[l] = _moms[l];
  }
  for (int l=j+1; l<=maxl; ++l) 
    omoms[l-1]=_moms[l];

  if (rpa.gen.Debugging()) {
    for (int l=0; l<maxl;++l) 
      msg.Out()<<" omoms["<<l<<"]="<<omoms[l]<<std::endl;
  }
}

void Combine_Table::FillTable(Leg **_legs,int _nlegs, int _nampl)
{
  msg.Tracking()<<"in Combine_Table::FillTable"<<std::endl;
  // store information
  legs  = _legs;
  nlegs = _nlegs;
  nampl = _nampl;

  // determine possible combinations and corresponding y_ij  if nlegs>4
  if (nlegs>4) {
    for (int i=2; i<nlegs; ++i) {   // *AS* 2 for outgoing only
      for (int j=i+1; j<nlegs; ++j) {
	// never combine "0&1" !
	if (j==1) j=2;
	// check if leg i is combinable with leg j in any graph
	for (int k=0;k<nampl;++k) {
	  if (Combinable(legs[k][i],legs[k][j])) {  
	    AddPossibility(i,j,k); // insert graph k with combination i&j in table
	  } 
	}	 
      }
    }
  }
  msg.Tracking()<<"out Combine_Table::FillTable"<<std::endl;
}


void Combine_Table::AddPossibility(int i, int j, int ngraph) {
  if (combinations.size()==0) {
    // add new "i&j" combination to empty table
    combinations.push_back(Combine_Data(i,j,0.,ngraph));
  }
  else if ((combinations.back().i==i)&&(combinations.back().j==j)) {
    // add graph only ("i&j" row exists already)
    combinations.back().graphs.push_back(ngraph);
  } 
  else {
    // add new "i&j" combination 
    combinations.push_back(Combine_Data(i,j,0.,ngraph));
  }
}


Combine_Table * Combine_Table::CalcJet(int nl, AMATOOLS::vec4d * _moms) {
  Combine_Table * ct=0;
  CD_List & cl=combinations;
  msg.Tracking()<<"in Combine_Table::CalcJet "<<std::endl;
  if (cl.size()==0) {
    return this;
  } 
  else {
    // change momenta to actual values    
    if (_moms!=0) {
      for (int l=0;l<nl;++l)
      moms[l]=_moms[l];
    }
     
    // calculate pt2ij and determine "best" combination
    double pt2min = sqr(rpa.gen.Ecms());
    for (CD_Iterator cit=cl.begin(); cit!=cl.end(); ++cit) {
      double pt2ij = cit->pt2ij = jf->PTij(moms[cit->i], moms[cit->j]);
      if (pt2ij<pt2min) {
	pt2min = pt2ij;
	cwin = cit;
      }
    } 
    msg.Debugging()<<" Winner:"<<(*cwin)<<std::endl;

    --nl;
    
    // if number of legs is still greater 4 Cluster once more
    // if number of legs equals 4, determine end situation
    if (nl>=4) {
      if (!cwin->down) {
	Leg ** alegs = new Leg*[cwin->graphs.size()];
	for (int k=0;k<cwin->graphs.size();++k) {
	  alegs[k]   = CombineLegs(legs[cwin->graphs[k]],cwin->i,cwin->j,nl);
	}
	vec4d * amoms;
	CombineMoms(moms,cwin->i,cwin->j,nl,amoms); // generate new momenta
	cwin->down=new Combine_Table(jf,amoms,this);
	cwin->down->FillTable(alegs,nl,cwin->graphs.size());   // initialise Combine_Table
      } 
      else {
	cwin->down->CombineMoms(moms,cwin->i,cwin->j,nl);
      }
      ct   = cwin->down->CalcJet(nl);
      gwin = cwin->down->gwin;
      gwin = cwin->graphs[gwin];   // translate back
    } 
    else {
      msg.Error()<<" ERROR:  nlegs < 4 !!!!!!!!!"<<std::endl;
      abort();
    }
  }
  msg.Debugging()<<"out Combine_Table::CalcJet  gwin="<<gwin<<"   table no.:"<<ct->no<<std::endl;
  return ct;
}

Combine_Table::~Combine_Table()
{
  //  std::cout<<this;
  msg.Tracking()<<" removing momenta"<<std::endl;
  delete [] moms;
 
  //  msg.Out()<<"WE STILL HAVE to delete SOME legs inside  combi!!!"<<std::endl; // !!!!!!!!!!!!!
  msg.Tracking()<<" removing legs"<<std::endl;

  //  std::cout<<" nampl="<<nampl<<std::endl;
  for (int k=0;k<nampl;++k) 
    delete [] legs[k];
  delete [] legs;
  msg.Tracking()<<" done "<<std::endl;
}

std::ostream& MOCAIC::operator<< (std::ostream& s ,Combine_Table * ct) {
  if (ct) {
    s<<std::endl<<" Combine_Table "<<ct->no<<" (up=";
    if (ct->up) s<<ct->up->no<<")"<<std::endl; else s<<"#)"<<std::endl;
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
