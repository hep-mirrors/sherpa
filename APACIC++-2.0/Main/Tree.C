#include "Tree.H"

using namespace APACIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using std::flush;
using std::endl;

namespace APACIC {
void StreamTree(std::ostream & s, Knot * mo)
{ 
  s<<(*mo);
  if (mo->left) StreamTree(s, mo->left);
  if (mo->right) StreamTree(s, mo->right);
}

std::ostream & operator<<(std::ostream & s, Tree * tree)
{
  s<<"-------------------"<<endl;
  Knot * mo = tree->GetRoot();
  if (mo) while (mo->prev) mo=mo->prev;
  if (mo) StreamTree(s, mo);
  s<<"-------------------"<<endl;
  return s;
}
}

//-----------------------------------------------------------------------
//--------------------------- Constructors ------------------------------
//----------------------------------------------------------------------- 
Tree::Tree() {
  p_knots     = new Knot_List;
  p_root      = 0;
  p_save_root = 0;
}

Tree::Tree(Tree * tree) {
  p_knots = new Knot_List;
  p_root  = NewKnot(tree->p_root);
  Links(p_root,tree->p_root);
  p_save_root = 0;
}

Tree::~Tree() {
  std::cout<<" deleting "<<this<<std::endl;
  Reset();
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

Knot * Tree::NewKnot(APHYTOOLS::Flavour fl, AMATOOLS::Vec4D p, double t, double x1) {
  Knot * newk = new Knot;
  p_knots->push_back(newk);
  newk->kn_no     = p_knots->size();
  Parton * newp   = new Parton(newk->kn_no,fl,p);
  newk->part      = newp;
  newk->t         = t;
  newk->x         = x1;
  newk->maxpt2    = 0.;
  newk->left      = 0;
  newk->right     = 0;
  newk->prev      = 0;
  if (!p_root) {
    p_root = newk;
    // msg.Debugging()<<"Knot::NewKnot : "<<newk<<" serves as new root."<<endl;
  }
  return newk;
};

Knot * Tree::NewKnot(Knot * ink) {
  Knot * newk = new Knot;
  p_knots->push_back(newk);
  newk->kn_no  = p_knots->size();
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
    // msg.Debugging()<<"Knot::NewKnot : "<<newk<<" serves as new root."<<endl;
  }
  newk->part      = new Parton(ink->part);
  return newk;
}

Knot * Tree::NewKnot(Parton * _inpart)
{
  Knot * newk     = new Knot;
  p_knots->push_back(newk);
  newk->kn_no     = p_knots->size();
  if (_inpart==NULL) {
    newk->part    = new Parton(newk->kn_no);
  }
  else {
    newk->part    = _inpart;
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


void Tree::Restore(Knot * start, Knot * stop) {
  if (start==stop) return;
  if (start==0) return;

  if (start->left) start->left->prev=0;
  if (start->right) start->right->prev=0;
  Restore(start->left,stop);
  Restore(start->right,stop);
  for (Knot_Iterator kit=p_knots->begin(); kit!=p_knots->end(); ++kit) {
    if ((*kit) == start) {
      delete (*kit);
      p_knots->erase(kit);
      return;
    }
  }
}


void Tree::Restore(Knot * in) {

  if (!(in->prev)) return;
  if (in->prev->left != in)  {
    ResetDaughters(in->prev->left);
    in->prev->left=0;
  }
  if (in->prev->right != in) {
    ResetDaughters(in->prev->right);
    in->prev->right=0;
  }
  Restore(in->prev);
  for (Knot_Iterator kit=p_knots->begin(); kit!=p_knots->end(); ++kit) {
    if ((*kit) == in->prev) {
      delete (*kit);
      p_knots->erase(kit);
      return;
    }
  }
}

void Tree::ResetDaughters(Knot * in) {
  if (!(in)) return;
  //msg.Debugging()<<"Tree::ResetDaughters."<<endl
  //		 <<"   Delete daughters of/and "<<in<<" "<<in->kn_no<<endl;

  
  if (in->left)  ResetDaughters(in->left);
  in->left=0;
  if (in->right) ResetDaughters(in->right);
  in->right=0;

  for (Knot_Iterator kit=p_knots->begin(); kit!=p_knots->end(); ++kit) {
    if ((*kit) == in) {
      delete (in);
      p_knots->erase(kit);
      return;
    }
  }
}

/*  obsolete new version see below
void Tree::Restore() {
  for (Knot_Iterator kit=p_knots->begin(); kit!=p_knots->end(); ++kit) {
    if ( ( (*kit)->part->Info() != 'G') && ( (*kit)->part->Info() != 'H') ) {
      //msg.Debugging()<<"Tree::Restore : Delete knot "<<(*kit)->kn_no
      //	     <<"  "<<(*kit)->part->Flav()<<"  "<<(*kit)->part->Info()<<endl;
      delete (*kit);
      p_knots->erase(kit);
    }
    else {
      //msg.Debugging()<<"Tree::Restore : Kept knot "<<(*kit)->kn_no
      //	     <<"  "<<(*kit)->part->Flav()<<"  "<<(*kit)->part->Info()<<endl;
      if ((*kit)->prev) {
	if (((*kit)->prev->part->Info() != 'G') && 
	    ((*kit)->prev->part->Info() != 'H')) (*kit)->prev  = 0;
      }
      if ((*kit)->left) {
	if (((*kit)->left->part->Info() != 'G') && 
	    ((*kit)->left->part->Info() != 'H')) (*kit)->left  = 0;
      }
      if ((*kit)->right) {
	if (((*kit)->right->part->Info() != 'G') && 
	    ((*kit)->right->part->Info() != 'H')) (*kit)->right = 0;
      }
    }
  }
}
*/

Knot * Tree::GetInitiator() {
  Knot * help;
  help = p_root;
  if (help)  while (help->prev) help = help->prev;
  return help;
}

Knot * Tree::GetRoot() { return p_root; } ;
  


//-----------------------------------------------------------------------
//--------------------------- New frames for the tree -------------------
//----------------------------------------------------------------------- 

void Tree::BoRo(AMATOOLS::Poincare & lorenz) 
{
  Knot * mo= GetRoot();
  if (mo)  while (mo->prev) mo = mo->prev;
  // msg.Debugging()<<"changing Knot "<<mo->kn_no<<" from "<<mo->part->Momentum()<<endl;
  mo->part->SetMomentum(lorenz*mo->part->Momentum());
  // msg.Debugging()<<"                to "<<mo->part->Momentum()<<endl;
  BoRoDaughters(lorenz,mo);
}

void Tree::BoRoDaughters(AMATOOLS::Poincare & lorenz, Knot * mo) 
{
  if (mo->left) {
    // msg.Debugging()<<"changing Knot "<<mo->left->kn_no<<" from "<<mo->left->part->Momentum()<<endl;
    mo->left->part->SetMomentum(lorenz*mo->left->part->Momentum());
    // msg.Debugging()<<"                to "<<mo->left->part->Momentum()<<endl;
    BoRoDaughters(lorenz,mo->left);
  }
  if (mo->right) {
    // msg.Debugging()<<"changing Knot "<<mo->right->kn_no<<" from "<<mo->right->part->Momentum()<<endl;
    mo->right->part->SetMomentum(lorenz*mo->right->part->Momentum());
    // msg.Debugging()<<"                to "<<mo->right->part->Momentum()<<endl;
    BoRoDaughters(lorenz,mo->right);      
  }
}

void Tree::BoRo(AMATOOLS::Poincare & lorenz, Knot * mo) 
{
  // msg.Debugging()<<"changing Knot "<<mo->kn_no<<" from "<<mo->part->Momentum()<<endl;
  mo->part->SetMomentum(lorenz*mo->part->Momentum());
  // msg.Debugging()<<"                to "<<mo->part->Momentum()<<endl;
  BoRoDaughters(lorenz,mo);
}

Knot * Tree::CopyKnot(Knot * a, Knot * prev) 
{
  if (!a) return 0;
  //cout<<" copy: "<<*a<<endl;
  Knot * nk = new Knot(a);
  nk->prev  = prev;
  nk->left  = CopyKnot(a->left,nk);
  nk->right = CopyKnot(a->right,nk);
  return nk;
}

void Tree::CopyBackKnot(Knot * a, Knot * b) 
{
  if (!b) return;
  //cout<<"a copy back "<<*a<<" from "<<*b;

  if (!a || a==b) {
    std::cout<<" ERROR: Tree::CopyBackKnot : no knot to copy to! "<<endl;
    std::cout<<" from "<<*b<<endl;
    std::cout<<"in"<<endl;
    Knot * b = p_save_root;
    if (b) {
      while (b->prev) {
	b = b->prev;
      }
    }
    StreamTree(std::cout,b);
    std::cout<<"to"<<endl;
    std::cout<<this<<endl;
    abort();
  }

  a->CopyData(b);
  if (b->left) CopyBackKnot(a->left,b->left);
  else a->left=0;
  if (b->right) CopyBackKnot(a->right,b->right);
  else a->right=0;
  if (!(b->prev))
    a->prev=0;
  //std::cout<<" to "<<*a<<endl;

}

void Tree::DeleteKnot(Knot * b) {
  if (!b) return;
  DeleteKnot(b->left);
  DeleteKnot(b->right);

  //std::cout<<" delete"<< *b<<endl;
  delete b;
}

void Tree::Store()
{
  Knot * help;
  help = p_save_root;
  if (help)  while (help->prev)  help = help->prev;
  DeleteKnot(help);  

  p_save_root = CopyKnot(GetInitiator(),0);
  //std::cout<<" Stored:"<<endl;
  //  StreamTree(std::cout,p_save_root);
  

  while (p_save_root->right) p_save_root = p_save_root->right;
}

void Tree::Restore()
{
  Knot * a = p_root;
  Knot * b = p_save_root;
  if (b) {
    while (b->prev) {
      b = b->prev;
      if (!a) {
	std::cout<<"ERROR in Tree::Restore() line 354"<<endl;
      }
      a = a->prev;
    }
  }
  if (!a) {
    std::cout<<"ERROR in Tree::Restore() line 360"<<endl;
  }
  CopyBackKnot(a,b);
}
