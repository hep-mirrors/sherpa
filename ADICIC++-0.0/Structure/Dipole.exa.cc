//bof
//Version: 1 ADICIC++-0.0/2004/06/05

//Implementation of the Dipole_Particle structure of Dipole.H.



//#include ""





//using;    //is already done in Dipole.C





//=============================================================================



int Dipole_Particle::s_count=0;    //so far not faked
int Dipole_Particle::s_maxcount=0;

const int& Dipole_Particle::InStore=Dipole_Particle::s_count;



//=============================================================================



Dipole_Particle::Dipole_Particle(const Dipole_Particle& dipa)
  : m_num(++s_maxcount), Name(m_num),
    m_typ(dipa.m_typ), m_tag(dipa.m_tag),
    m_pac(dipa.m_pac), Parton(m_pac),
    m_tow(list<Dipole*>()) {

  //It cannot belong to any dipole since it is just created, right now.
  //Thus, we use an empty tower.

  ++s_count; this->SetPacNum(); this->cc_info();
}





Dipole_Particle::Dipole_Particle(const Dipole_Quark_Base& Q, const Vec4D& P)
  : m_num(++s_maxcount), Name(m_num),
    m_typ(Positive), m_tag(Q.Tag()),
    m_pac( Particle(m_num,Q(),P) ), Parton(m_pac),
    m_tow(std::list<Dipole*>()) {

  ++s_count; this->bc_info();
}





Dipole_Particle::Dipole_Particle(const Dipole_Antiquark_Base& A,
				 const Vec4D& P)
  : m_num(++s_maxcount), Name(m_num),
    m_typ(Negative), m_tag(A.Tag()),
    m_pac( Particle(m_num,A(),P) ), Parton(m_pac),
    m_tow(std::list<Dipole*>()) {

  ++s_count; this->bc_info();
}





Dipole_Particle::Dipole_Particle(const Dipole_Gluon_Base& G, const Vec4D& P)
  : m_num(++s_maxcount), Name(m_num),
    m_typ(Nil), m_tag(G.Tag()),
    m_pac( Particle(m_num,G(),P) ), Parton(m_pac),
    m_tow(std::list<Dipole*>()) {

  ++s_count; this->bc_info();
}





Dipole_Particle::~Dipole_Particle() {

  if(m_tow.empty()) {

    --s_count; this->sd_info();
    if(m_num && m_num==s_maxcount) --s_maxcount;
    if(!s_count) s_maxcount=0;

  } else {

    this->Test();    //Pointer handling test.
    Dipole_Particle* pcopy=this->Copy();
#ifdef DIPOLE_OUTPUT
    cout<<"    1.Test: Old tower: Empty? ..."<<m_tow.empty();
    cout<<"    New tower: Empty? ..."<<pcopy->m_tow.empty()<<endl;
#endif
    list<Dipole*>::iterator iter=pcopy->m_tow.begin();
    pcopy->m_tow.splice(iter,this->m_tow);
    iter=pcopy->m_tow.begin();
#ifdef DIPOLE_OUTPUT
    cout<<"    2.Test: Old tower: Empty? ..."<<m_tow.empty();
    cout<<"    New tower: Empty? ..."<<pcopy->m_tow.empty()<<endl;
#endif
    Dipole::Gate gate;
    gate(*iter,this,pcopy,true);
    ++iter;
    for(; iter!=pcopy->m_tow.end(); ++iter) gate(*iter,this,pcopy,false);

  }

}





Dipole_Particle& Dipole_Particle::operator=(const Dipole_Particle& dipa) {

  if(this==&dipa) return *this;

  m_tag=dipa.m_tag;
  m_pac=dipa.m_pac;
  this->SetPacNum();
  //Towering! All dipoles carrying this particle have to be updated.
  for(list<Dipole*>::iterator it=m_tow.begin(); it!=m_tow.end(); ++it) {
    (*it)->UpdateType();
    (*it)->UpdateMass();
  }

  return *this;

}





const bool Dipole_Particle::operator==(const Dipole_Particle& dp) const {
  cerr<<"\nSorry :o( Method has not been implemented yet.\n";
  assert(0); exit(1);
}





void Dipole_Particle::ShowDipoles() const {
  cout<<endl<<"[Tower]"; this->nm_info(); cout<<endl; 
  cout<<string(14,'=')<<endl<<"dip  ph  type"<<endl;
  for(list<Dipole*>::const_iterator it=m_tow.begin(); it!=m_tow.end(); ++it)
    cout<<setiosflags(ios::left)
	<<setw(5)<<(*it)->Name<<setw(4)<<(*it)->PointerHandling()
	<<setw(5)<<(*it)->IsType()<<endl
	<<resetiosflags(ios::left);
  cout<<string(14,'=')<<endl;
}





//Available:
//(g,Nil)[2] (g,Pos) (g,Neg) (q,Pos) (q,Nil) (qbar,Neg) (qbar,Nil)
//Not available:
//(q,Neg) (qbar,Pos)





const Particle& Dipole_Particle::Quarkize(const Dipole_Quark_Base& Q) {
  //(g,Nil)->(q,Nil)
  //(g,Pos)->(q,Pos)
  if( m_tag==0 && (m_typ==Nil || m_typ==Positive) ) {    //SM gluons!
    if(Q.Tag()==1) {    //SM quarks!
      m_tag=Q.Tag();
      m_pac.SetFlav(Q());
      //Towering!
      for(list<Dipole*>::iterator it=m_tow.begin(); it!=m_tow.end(); ++it)
	(*it)->UpdateType();    //updates and checks the corresponding dipoles.
    }
  }
  return m_pac;
}





const Particle& Dipole_Particle::Antiquarkize(const Dipole_Antiquark_Base& A) {
  //(g,Nil)->(qbar,Nil)
  //(g,Neg)->(qbar,Neg)
  if( m_tag==0 && (m_typ==Nil || m_typ==Positive) ) {    //SM gluons!
    if(A.Tag()==-1) {    //SM antiquarks!
      m_tag=A.Tag();
      m_pac.SetFlav(A());
      //Towering!
      for(list<Dipole*>::iterator it=m_tow.begin(); it!=m_tow.end(); ++it)
	(*it)->UpdateType();
    }
  }
  return m_pac;
}





const Particle& Dipole_Particle::Gluonize() {
  //(q,Pos/Nil)->(g,Pos/Nil)
  //(qbar,Neg/Nil)->(qbar,Neg/Nil)
  if( m_tag== 1 && (m_typ==Positive || m_typ==Nil) ||
      m_tag==-1 && (m_typ==Negative || m_typ==Nil) ) {    //SM quarks!
    m_tag=info.gluon.g.Tag();
    m_pac.SetFlav(info.gluon.g());
    //Towering!
    for(list<Dipole*>::iterator it=m_tow.begin(); it!=m_tow.end(); ++it)
      (*it)->UpdateType();
  }
  return m_pac;
}



//=============================================================================



Dipole_Particle* Dipole_Particle::Copy() const {

  //Note: Dipole_Particle s have always an empty tower.
  //Thus: Copy is used only for derived types, see ~Dipole_Particle().

  this->cp_info();

  Dipole_Particle* pcopy;
  switch(m_typ) {
  case Positive: pcopy=new Dipole_Branch(*(Dipole_Branch*)this); break;
  case Negative: pcopy=new Dipole_Antibranch(*(Dipole_Antibranch*)this); break;
  case Nil:      pcopy=new Dipole_Glubranch(*(Dipole_Glubranch*)this); break;
  default: cerr<<"\nBug :o( in:\n"; assert(0);
  }

  assert(pcopy);

  pcopy->sd_info();

  --s_count;
  --s_maxcount;
  pcopy->m_num=this->m_num;
  pcopy->SetPacNum();

  return pcopy;

}



//=============================================================================



//The interface structures.
//=========================



Dipole_Branch::Dipole_Branch(const Particle& par)
  : Dipole_Particle(Positive,par) {
  assert(par.Flav().IsQuark() && !par.Flav().IsAnti());
  //later: this->SetTag()=...;
  //test: this->SetTag()=9; -> Dipole::Type assertion failed
}
Dipole_Antibranch::Dipole_Antibranch(const Particle& par)
  : Dipole_Particle(Negative,par) {
  assert(par.Flav().IsQuark() && par.Flav().IsAnti());
  //later: this->SetTag()=...;
}
Dipole_Glubranch::Dipole_Glubranch(const Particle& par)
  : Dipole_Particle(Nil,par) {
  assert(par.Flav().IsGluon());
  //later: this->SetTag()=...;
}



//=============================================================================





//eof
