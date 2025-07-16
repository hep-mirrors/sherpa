#include "AddOns/Analysis/Triggers/Midpoint_Cone.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/MathTools.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;


Midpoint_Cone::Midpoint_Cone(ATOOLS::Particle_Qualifier_Base * const qualifier,
		  int reco, double f) : 
  Jet_Algorithm_Base(qualifier), _M_reco(reco), _M_midpoint(0), m_f(f)
{
}

Midpoint_Cone::~Midpoint_Cone() 
{
}


void Midpoint_Cone::AddToKtlist(double kt2) {
  if (p_kts) {
    p_kts->push_back(kt2);
  }
}

void Midpoint_Cone::AddToJetlist(const Vec4D & mom, bool bf) {
  if (p_jets) {
    if(!bf) p_jets->push_back(new Particle(p_jets->size(),Flavour(kf_jet),mom));
    else    p_jets->push_back(new Particle(p_jets->size(),Flavour(kf_bjet),mom));
  }
}

bool Midpoint_Cone::ConstructJets(const Particle_List *pl, Particle_List * jets,
				     std::vector<double> * kt2, double R)
{
  // assume empty containers
  p_jets = jets;
  p_kts  = kt2;

  // create vector list from particle list
  int n=0;
  _M_pp.clear();
  _M_pp.push_back(Vec4D(0.,0.,0.,0.));
  bool * bflag = new bool[pl->size()];
  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    if ((*p_qualifier)(*it)) {
      //    if (!(*it)->Flav().IsLepton()) {
      _Vector p((*it)->Momentum());
      _M_pp.push_back(p); 
      bflag[n] = (((*it)->Flav()).Kfcode()==kf_b)&& !m_bflag;
      ++n;
    }
  }

  // cluster
  _M_do_clustering(R, m_f);

  _M_pp.clear();
  delete [] bflag;

  // finalize (sort and release used containers)
  SortPT(p_jets);

  p_jets=0;
  p_kts =0;

  return true;
}


void Midpoint_Cone::_M_iterate_cone(double R, const _Vector& c)
  {
    //---- calculate the centroid ----
    _Proto __prev, __curr;
    _Proto *prev = &__prev, *curr = &__curr;
    
    //---- check the cone ----
    double dy, dphi, R2 = R*R;
    double y = c.eta, phi = c.phi;
    unsigned int nt = _M_pp.size()-1;

    for(unsigned int i = 1; i <= nt; i++) {
      dy = y - _M_pp[i].eta;
      dphi = _S_dphi(phi - _M_pp[i].phi);
      if(dy*dy+dphi*dphi <= R2) {
	curr->pars.push_back(i);
	this -> _M_assadd(curr->psum, _M_pp[i]);
      }
    }
    
    //----- this cone contains a single parton -----
    if(curr->pars.size() == 1) {
       this -> _M_was_it_already_found(curr);
      return;
    }
    int ccc=0;
    //----- iterate the cone -----
    do {
      ccc++;
      _Proto *tmp = prev; 
      prev = curr; curr = tmp;
      
      y = prev->psum.eta, phi = prev->psum.phi;
      curr->psum = _Vector();
      curr->pars.clear();
      
      for(unsigned int i = 1; i <= nt; i++) {
	dy = y - _M_pp[i].eta;
	dphi = _S_dphi(phi - _M_pp[i].phi);
	if(dy*dy+dphi*dphi <= R2) {
	  curr -> pars.push_back(i);
	  this -> _M_assadd(curr->psum, _M_pp[i]);
	}
      }
      if(curr->pars.size() == 0) return;
       if(ccc==100){
	 cout<<"Midpoint_Cone::_M_iterate_cone does not converge! "<<_M_reco<<endl;
// 	 abort();
       }
//         if(ccc>20){
//         cout<<" content of preproto: "<<endl;
//          std::list<unsigned int>::const_iterator b=curr->pars.begin();
//          for(;b!=curr->pars.end();b++) cout<<*b<<"("<<_M_pp[*b].phi<<","<<_M_pp[*b].eta<<") "<<_M_pp[*b].p<<endl;;
//         cout<<ccc<<"  centroid:           "<<curr->psum.p<<"  "<<curr->psum.phi<<","<<curr->psum.eta<<endl;
// 	}
    } while(curr->pars != prev->pars && ccc<100);
//        cout<<" content of proto: "<<endl;
//        std::list<unsigned int>::const_iterator b=curr->pars.begin();
//        for(;b!=curr->pars.end();b++) cout<<*b<<"("<<_M_pp[*b].phi<<","<<_M_pp[*b].eta<<") "<<_M_pp[*b].p<<endl;;
//        cout<<"  centroid:           "<<curr->psum.p<<"  "<<curr->psum.phi<<","<<curr->psum.eta<<endl;
    
    this -> _M_was_it_already_found(curr);
    return;
  }
 
  
  void Midpoint_Cone::_M_was_it_already_found(const _Proto *curr)
  {
    std::list<_Proto>::iterator b = _M_proto.begin();
    
    while(b != _M_proto.end()) {
      if(b->pars == curr->pars) return;
      b++;
    }
    
    // else add to the list of proto-jets
    _M_proto.push_back(*curr);
  }

  void Midpoint_Cone::_M_do_clustering(double R, double f)
  {
    //  number of particles
    unsigned int nt = _M_pp.size()-1;
    
    //  initialize
    _M_proto.clear();
    
    if(_M_midpoint == 0) {
      //   find the primary proto-jets
      for(unsigned int i = 1; i <= nt; i++) 
	{
	this -> _M_iterate_cone(R, _M_pp[i]);
	}
      //   add the midpoints
      _Vector pij;
      unsigned int np = _M_proto.size();
      std::list<_Proto>::iterator pj, pi = _M_proto.begin();
      
      for(unsigned int i = 0; i < np; i++, pi++) {
	pj = pi; pj++;
	for(unsigned int j = i+1; j < np; j++, pj++) {
	  pij = pi->psum; this->_M_assadd(pij, pj->psum);
	  this -> _M_iterate_cone(R, pij);
	}
      }      
    } else if(_M_midpoint == 1) {
      cout<<"getting seeds for "<<nt<<" partons"<<endl;
      //  find the all possible centroid
      for(unsigned int c = 1; c <= nt; c++) {
	unsigned int i, *p = new unsigned int[nt];
	
	for(i = 0; i < c; i++) p[i] = i+1;
	this->_M_check_trial_cone(c, p, R);
	
	while(p[0] != nt-c+1){
	  bool flag = false;
	  
	  i = 0; ++p[c-1];
	  while(i < c && p[c-i-1] > nt-i) {
	    ++p[c-i-2]; flag = true; i++;
	  }
	  
	  if(flag)
	    for(unsigned int j = c-i; j < c; j++)
	      p[j] = p[j-1] + 1;
	  
	  this->_M_check_trial_cone(c, p, R);
	}
	delete[] p;
      }
    } else throw "undefined cone variant";
    
    //---- analizing the protojet ----
    unsigned int npj, ns = 0;
    std::list<_Proto>::iterator i, b;
    
    
    while((npj = _M_proto.size()) > 0) {
      //---- sort the proto-jets by transverse energy ----
      if(npj > 1) _M_proto.sort(pT_sort());
      
      //---- checking the shared partons ----
      ns = 0; i = b = _M_proto.begin(); i++;
      for(; i != _M_proto.end(); i++)
	if((ns = this->_M_check_shared_towers(*b, *i)) > 0) break;
//       cout<<"_M_do_clustering: "<<_M_reco<<" ";
//       if (i != _M_proto.end())cout<<"i: "<<i->psum.p;
//       if (b != _M_proto.end())cout<<"b: "<<b->psum.p;
//       cout<<endl;
      if(ns == 0) {
	//---- no shared parton ----
	switch(_M_reco) {
	case 0: 
	  AddToJetlist(b->psum.p,false); 
	  AddToKtlist(b->psum.p.PPerp2());
	  break;
	case 1:
	  {
	    double pT = b->psum.pt, ei = b->psum.eta, fi = b->psum.phi;
	    Vec4D p(pT*Vec4D(std::cosh(ei),std::cos(fi),std::sin(fi),std::sinh(ei))); 
	    AddToJetlist(p,false); 
	    AddToKtlist(pT*pT);	    
	  }
	  break;
	}
	
	_M_proto.pop_front();
      } else this -> _M_split_merge(f, b, i);
    }
  }
  

  void Midpoint_Cone::
  _M_check_trial_cone(unsigned int n, unsigned int *idx, double R)
  {
    //---- calculate the centroid ----
    _Proto cone;
    
    for(unsigned int i = 0; i < n; i++)
      this -> _M_assadd(cone.psum, _M_pp[idx[i]]);
    
    //---- check the cone ----
    double dy, dphi, R2 = R*R;
    double y = cone.psum.eta, phi = cone.psum.phi;
    
    unsigned int nt = _M_pp.size()-1;
    for(unsigned int i = 1; i <= nt; i++) {
      dy = y - _M_pp[i].eta;
      dphi = _S_dphi(phi - _M_pp[i].phi);
      if(dy*dy+dphi*dphi <= R2) 
	cone.pars.push_back(i);
    }
    
    //---- compare the parton lists ----
    if(n != cone.pars.size()) return;
    else {
      std::list<unsigned int>::const_iterator b = cone.pars.begin();
      for(unsigned int i = 0; i < n; i++)
	if(idx[i] != *(b++)) return;
    }
    
    //---- this is a stable, add the centroid to list of the protojets ----
    _M_proto.push_back(cone);
  }
  

  unsigned int Midpoint_Cone::
  _M_check_shared_towers(const _Proto& x, const _Proto& y)
  {
    typedef std::list<unsigned int>::const_iterator const_iterator;
    
    if(x.pars.size() == 1 && y.pars.size() == 1) 
      return 0;
    
    _M_overlap.clear();
    for(const_iterator i = x.pars.begin(); i != x.pars.end(); i++)
      for(const_iterator j = y.pars.begin(); j != y.pars.end(); j++)
	if(*i == *j) _M_overlap.push_back(*i);
    
    return _M_overlap.size();
  }


  void Midpoint_Cone::
  _M_split_merge(double f, std::list<_Proto>::iterator b, std::list<_Proto>::iterator i)
  {
    //---- start the merge/split algorithm -----
    _Vector ps;
//     cout<<"_M_split_merge: "<<f<<" "<<b->psum.p<<" "<<i->psum.p<<endl;
    std::list<unsigned int>::const_iterator io;
    for(io = _M_overlap.begin(); io != _M_overlap.end(); io++)
      this -> _M_assadd(ps, _M_pp[*io]);
//     cout<<ps.p<<endl;    
    double Es = ps.pt, En = i->psum.pt;
    if(Es > f*En) {
      //---- merge ----
      for(io = _M_overlap.begin(); io != _M_overlap.end(); io++)
	i->pars.remove(*io);
      
      b->pars.merge(i->pars);
      _M_proto.erase(i);
      
      //---- recalculate the centroid ----
      io = b->pars.begin(); b->psum = _M_pp[*io]; io++;
      for(; io != b->pars.end(); io++)
	this->_M_assadd(b->psum, _M_pp[*io]);
    } else {
      //---- split ----
      bool b_flag = false, i_flag = false;
      double b_y = b->psum.eta, b_phi = b->psum.phi;
      double i_y = i->psum.eta, i_phi = i->psum.phi;
      double s_y, b_dy, b_dphi, s_phi, i_dy, i_dphi;
      
      for(io = _M_overlap.begin(); io != _M_overlap.end(); io++) {
	s_y = _M_pp[*io].eta; s_phi = _M_pp[*io].phi;
	
	b_dy = b_y-s_y; b_dphi = _S_dphi(b_phi-s_phi);
	i_dy = i_y-s_y; i_dphi = _S_dphi(i_phi-s_phi);
	
	if(b_dy*b_dy+b_dphi*b_dphi < i_dy*i_dy+i_dphi*i_dphi) { 
	  i_flag = true;
	  i->pars.remove(*io);
	} else {
	  b_flag = true;
	  b->pars.remove(*io);
	}
      }
      
      //---- recalculate the centroids ----
      if(i_flag) {
	bool inflag = false;
	std::list<_Proto>::iterator ip;
	
	//---- Is this proto-jet already in the list? ---- 
	for(ip = _M_proto.begin(); ip != _M_proto.end(); ip++)
	  if(ip != i && ip->pars == i->pars) inflag = true;
	
	//---- Delete if it is in. Otherwise recalculate the centroid. ----
	if(inflag) _M_proto.erase(i);
	else {
	  io = i->pars.begin(); i->psum = _M_pp[*io]; io++;
	  for(; io != i->pars.end(); io++) 
	    this->_M_assadd(i->psum, _M_pp[*io]);
	}
      }
      
      if(b_flag) {
	bool inflag = false;
	std::list<_Proto>::iterator ip;

	//---- Is this proto-jet already in the list? ---- 
	for(ip = _M_proto.begin(); ip != _M_proto.end(); ip++)
	  if(ip != b && ip->pars == b->pars) inflag = true;
	
	//---- Delete if it is in. Otherwise recalculate the centroid. ----
	if(inflag) _M_proto.erase(b);
	else {
	  io = b->pars.begin(); b->psum = _M_pp[*io]; io++;
	  for(; io != b->pars.end(); io++)
	    this->_M_assadd(b->psum, _M_pp[*io]);
	}
      }
    }
  }


  double Midpoint_Cone::_S_dphi(double phi)
  {
    static const double twopi = 2.*M_PI;
    
    if(phi >= M_PI) phi = std::fmod(M_PI+phi, twopi) - M_PI;
    else if(phi < -M_PI) phi = -std::fmod(M_PI-phi, twopi) + M_PI;
    
    return phi;
  }
  
  void Midpoint_Cone::_M_assadd(_Vector& p1, const _Vector& p2) const
  {
    switch(_M_reco) {
    case 0: 
      p1.p += p2.p;
      p1.pt = p1.p.PPerp();
      p1.eta = p1.p.Y();
      p1.phi = p1.p.Phi();
      break;
    case 1:
      { //--- pT weighted schemes ---
        double pi = p1.pt,  pj = p2.pt;
        double ei = p1.eta, ej = p2.eta; 
        double fi = p1.phi, fj = p2.phi;
        
        //--- weighted sum ---
        p1.pt = pi + pj;
	p1.eta = (pi*ei + pj*ej)/(pi+pj);
	p1.phi = _S_dphi(fi + pj*_S_dphi(fj-fi)/(pi+pj));
      }
      break;
    }
  }
