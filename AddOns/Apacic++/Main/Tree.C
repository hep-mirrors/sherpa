#include "AddOns/Apacic++/Main/Tree.H"
#include "ATOOLS/Org/Message.H"

using namespace APACIC;
using namespace ATOOLS;
using std::flush;
using std::endl;

namespace APACIC {
void StreamTree(std::ostream & s, Knot * mo)
{ 
  s<<(*mo)<<"\n";
  if (mo->left) StreamTree(s, mo->left);
  if (mo->right) StreamTree(s, mo->right);
}

std::ostream &operator<<(std::ostream & s,const Tree &tree)
{
  s<<"-------------------"<<endl;
  Knot * mo = tree.GetRoot();
  if (mo) while (mo->prev) mo=mo->prev;
  if (mo) StreamTree(s, mo);
  s<<"-------------------"<<endl;
  return s;
}
}

// ------------------ static ----------------
Knot_List * Tree::s_knots=NULL;

//-----------------------------------------------------------------------
//--------------------------- Constructors ------------------------------
//----------------------------------------------------------------------- 


Tree::Tree() {
  if (!s_knots) s_knots = new Knot_List;
  p_root      = 0;
  p_save_root = 0;
}

Tree::Tree(Tree * tree) {
  if (!s_knots) s_knots = new Knot_List;
  p_root  = NewKnot(tree->p_root);
  Links(p_root,tree->p_root);
  p_save_root = 0;
}

void Tree::ResetKnots() {
  if (!s_knots) return;
  for (Knot_List::iterator kit=s_knots->begin(); kit!=s_knots->end(); ++kit) {
    delete (*kit); 
  }
  s_knots->erase(s_knots->begin(),s_knots->end());
  p_root = 0;
}


Tree::~Tree() {
  Reset();
  if (s_knots) {
    delete s_knots;
    s_knots=NULL;
  }

  Knot * help;
  help = p_save_root;
  if (help) while (help->prev) help = help->prev;
  DeleteKnot(help);
}

//-----------------------------------------------------------------------
//--------------------------- Initialize New knots ----------------------
//----------------------------------------------------------------------- 

void Tree::Links(Knot * act,Knot * ink) {
  if (ink->left) {
    act->left = NewKnot(ink->left);
    Links(act->left,ink->left);
  }
  if (ink->right) {
    act->right = NewKnot(ink->right);
    Links(act->right,ink->right);
  }
  if (ink->prev) {
    act->prev = NewKnot(ink->prev);
    Links(act->prev,ink->prev);
  }
}

Knot * Tree::NewKnot(ATOOLS::Flavour fl, ATOOLS::Vec4D p, double t, double x1) {
  Knot * newk = new Knot;
  s_knots->push_back(newk);
  newk->kn_no     = s_knots->size();
  Particle * newp = new Particle(newk->kn_no,fl,p);
  newk->part      = newp;
  newk->t         = t;
  newk->x         = x1;
  newk->maxpt2    = 0.;
  newk->left      = 0;
  newk->right     = 0;
  newk->prev      = 0;
  if (!p_root) {
    p_root = newk;
  }
  return newk;
}

Knot * Tree::NewKnot(Knot * ink) {
  Knot * newk = new Knot;
  s_knots->push_back(newk);
  newk->kn_no  = s_knots->size();
  newk->t      = ink->t;
  newk->tout   = ink->tout;
  newk->maxpt2 = ink->maxpt2;
  newk->x      = ink->x;
  newk->E2     = ink->E2;
  newk->costh  = ink->costh;
  newk->phi    = ink->phi;
  newk->thcrit = ink->thcrit;

  newk->left      = 0;
  newk->right     = 0;
  newk->prev      = 0;
  if (!p_root) {
    p_root = newk;
  }
  newk->part      = new Particle(*ink->part);
  return newk;
}

Knot * Tree::NewKnot(Particle * _inpart)
{
  Knot * newk     = new Knot;
  s_knots->push_back(newk);
  newk->kn_no     = s_knots->size();
  if (_inpart==NULL) {
    newk->part    = new Particle(newk->kn_no);
  }
  else {
    newk->part    = new Particle(*_inpart);
  }
  newk->left      = 0;
  newk->right     = 0;
  newk->prev      = 0;
  if (!p_root) {
    p_root = newk;
  }
  return newk;
}

//-----------------------------------------------------------------------
//--------------------------- Resetting the tree ------------------------
//----------------------------------------------------------------------- 
/*
void Tree::ResetDaughters(Knot * in) {
  if (!(in)) return;
  
  if (in->left)  ResetDaughters(in->left);
  in->left=0;
  if (in->right) ResetDaughters(in->right);
  in->right=0;

  for (Knot_Iterator kit=s_knots->begin(); kit!=s_knots->end(); ++kit) {
    if ((*kit) == in) {
      delete (in);
      s_knots->erase(kit);
      return;
    }
  }
}
*/


Knot * Tree::GetInitiator() const {
  Knot * help;
  help = p_root;
  if (help)  while (help->prev) help = help->prev;
  return help;
}

void Tree::UpdateDaughters(Knot * mo)
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"):\n";
  msg_Indent();
  if (mo->part->Momentum()==Vec4D(0.,0.,0.,0.)) return;
  mo->E2 = sqr(mo->part->Momentum()[0]);
  if (mo->left==NULL) return;
  Knot *d1(mo->left), *d2(mo->right);
  if (d1->part->Momentum()==Vec4D(0.,0.,0.,0.)) return;
  double z=d1->part->Momentum()[0]/mo->part->Momentum()[0];
  if (z>0.0 && z<1.0) mo->z=z;
  UpdateDaughters(d1);
  UpdateDaughters(d2);
}


//-----------------------------------------------------------------------
//--------------------------- New frames for the tree -------------------
//----------------------------------------------------------------------- 

void Tree::BoRo(ATOOLS::Poincare & lorenz) 
{
  Knot * mo(GetRoot());
  while (mo->prev) mo = mo->prev;
  mo->part->SetMomentum(lorenz*mo->part->Momentum());
  BoRoDaughters(lorenz,mo);
}

void Tree::BoRoDaughters(ATOOLS::Poincare & lorenz, Knot * mo) 
{
  if (mo->left) {
    mo->left->part->SetMomentum(lorenz*mo->left->part->Momentum());
    BoRoDaughters(lorenz,mo->left);
    mo->right->part->SetMomentum(lorenz*mo->right->part->Momentum());
    BoRoDaughters(lorenz,mo->right);      
  }
}

void Tree::BoRo(ATOOLS::Poincare & lorenz, Knot * mo) 
{
  mo->part->SetMomentum(lorenz*mo->part->Momentum());
  BoRoDaughters(lorenz,mo);
}

bool Tree::Restore(Knot *const k,Knot *const r) const
{
  if (r==NULL) return false;
  if (k->kn_no==r->kn_no) {
    k->CopyData(r);
    return true;
  }
  if (Restore(k,r->left)) return true;
  if (Restore(k,r->right)) return true;
  return false;
}

bool Tree::Restore(Knot *const k) const
{
  if (p_save_root==NULL) {
    msg_Error()<<METHOD<<"(): No Storage."<<std::endl;
    return false;
  }
  Knot *ini(p_save_root);
  while (ini->prev) ini=ini->prev;
  if (!Restore(k,ini)) {
    msg_Error()<<METHOD<<"(): Knot "<<k->kn_no<<" not found."<<std::endl;
    return true;
  }
  return false;
}

Knot * Tree::CopyKnot(Knot * a, Knot * prev) 
{
  if (!a) return 0;
  Knot * nk = new Knot(a);
  nk->prev  = prev;
  nk->left  = CopyKnot(a->left,nk);
  nk->right = CopyKnot(a->right,nk);
  return nk;
}

void Tree::CopyBackKnot(Knot * a, Knot * b) 
{
  if (!b) return;

  if (!a || a==b) {
    msg_Error()<<" Error in  Tree::CopyBackKnot "<<a<<" "<<b<<std::endl;
    Knot * b = p_save_root;
    if (b) {
      while (b->prev) {
	b = b->prev;
      }
    }
    StreamTree(std::cout,b);
    abort();
  }

  if (a->decay!=NULL) {
    Knot *l(a->left), *d(a->decay);
    while (d==l->decay) {
      d=l->decay;
      l=l->left;
    }
    if (l->prev!=a) {
      DeleteKnot(a->right);
      a->right=l->prev->right;
      l->prev->right=l->prev->left=NULL;
      DeleteKnot(a->left);
      a->left=l;
    }
  }
  a->CopyData(b);
  if (b->left) CopyBackKnot(a->left,b->left);
  else {
    DeleteKnot(a->left);
    a->left=0;
  }
  if (b->right) CopyBackKnot(a->right,b->right);
  else {
    DeleteKnot(a->right);
    a->right=0;
  }
  if (!(b->prev)) a->prev=0;
}

void Tree::DeleteKnot(Knot * b) {
  if (!b) return;
  DeleteKnot(b->left);
  DeleteKnot(b->right);
  for (Knot_List::iterator kit(s_knots->begin());
       kit!=s_knots->end();++kit) {
    if (*kit==b) {
      s_knots->erase(kit);
      kit=s_knots->begin();
    }
  }
  delete b;
}

void Tree::ClearStore()
{
  Knot * help;
  help = p_save_root;
  if (help)  while (help->prev)  help = help->prev;
  DeleteKnot(help); 
  p_save_root = 0;
}

void Tree::Store()
{
  Knot * help;
  help = p_save_root;
  if (help)  while (help->prev)  help = help->prev;
  DeleteKnot(help);  

  p_save_root = CopyKnot(GetInitiator(),0);

  if (p_root!=GetInitiator())
    while (p_save_root->right) p_save_root = p_save_root->right;
}

void Tree::Restore()
{
  Knot * a = p_root;
  Knot * b = p_save_root;
  if (b) {
    while (b->prev) {
      b = b->prev;
      a = a->prev;
    }
    if (a->prev) {
      Knot *ini(a->prev);
      if (ini->left==a) ini->left=NULL;
      else ini->right=NULL;
      DeleteKnot(ini);
      a->prev=NULL;
    }
  }
  CopyBackKnot(a,b);
}

bool Tree::SingleCheckStructure(Knot *mo, Knot*gr, bool fixit)
{
  if (!mo) return true;
  if (!(mo->prev==gr)) {
    if (!fixit) {
      msg_Error()<<"ERROR in Tree::SingleCheckStructure :"<<std::endl
		 <<"   Relation of "<<(*gr)<<" and "<<(*mo)<<" badly defined."<<std::endl
		 <<"   Return false."<<std::endl;
      return false;
    }
    else {
      msg_Error()<<"WARNING in Tree::SingleCheckStructure :"<<std::endl
		 <<"   Relation of "<<(*gr)<<" and "<<(*mo)<<" badly defined."<<std::endl
		 <<"   Try to repair it and hope for the best."<<std::endl;
      mo->prev=gr;
    }
  }
  return (SingleCheckStructure(mo->left,mo,fixit) && SingleCheckStructure(mo->right,mo,fixit));
}

bool Tree::CheckStructure(bool fixit)
{
  Knot * mo(GetInitiator());
  return (SingleCheckStructure(mo->left,mo,fixit) && SingleCheckStructure(mo->right,mo,fixit));
}

bool Tree::CheckMomentumConservation() const 
{
  bool success(CheckMomentumConservation(GetInitiator()));
  return success;
}

bool Tree::CheckMomentumConservation(Knot *const knot) const 
{
  if (knot==NULL) return true;
  msg_Indent();
  bool success(true);
  static double accu(sqrt(rpa->gen.Accu()));
  if (knot->left!=NULL && !IsEqual(knot->part->Momentum(),Vec4D(),accu)) {
    msg_Debugging()<<"fmc check "<<knot->kn_no<<"\n"; 
    Vec4D p(knot->part->Momentum());
    Vec4D p1(knot->left->part->Momentum()), p2(knot->right->part->Momentum());
    if (!IsEqual(p,p1+p2,accu)) {
      msg_Error()<<METHOD<<"(): Four momentum not conserved in knot "
		 <<knot->kn_no<<"\n   p      = "<<p
		 <<"\n   p_miss = "<<(p-p1-p2)
		 <<"\n   p1     = "<<p1<<"\n   p2     = "<<p2<<std::endl;
      success=false;
    }
  }
  if (!CheckMomentumConservation(knot->left)) success=false;
  if (!CheckMomentumConservation(knot->right)) success=false;
  return success;
}
