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
  while (mo->prev) mo=mo->prev;
  if (mo) StreamTree(s, mo);
  s<<"-------------------"<<endl;
  return s;
}
}

//-----------------------------------------------------------------------
//--------------------------- Constructors ------------------------------
//----------------------------------------------------------------------- 
Tree::Tree() {
  knots   = new Knot_List;
  root    = 0;
}

Tree::Tree(Tree * tree) {
  knots = new Knot_List;
  root  = NewKnot(tree->root);
  Links(root,tree->root);
}

Tree::~Tree() {
  Reset();
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
  knots->push_back(newk);
  newk->kn_no     = knots->size();
  Parton * newp   = new Parton(newk->kn_no,fl,p);
  newk->part      = newp;
  newk->t         = t;
  newk->x         = x1;
  newk->maxpt2    = 0.;
  newk->left      = 0;
  newk->right     = 0;
  newk->prev      = 0;
  if (!root) {
    root = newk;
    // msg.Debugging()<<"Knot::NewKnot : "<<newk<<" serves as new root."<<endl;
  }
  return newk;
};

Knot * Tree::NewKnot(Knot * ink) {
  Knot * newk = new Knot;
  knots->push_back(newk);
  newk->kn_no  = knots->size();
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
  if (!root) {
    root = newk;
    // msg.Debugging()<<"Knot::NewKnot : "<<newk<<" serves as new root."<<endl;
  }
  newk->part      = new Parton(ink->part);
  return newk;
};


Knot * Tree::NewKnot()
{
  Knot * newk     = new Knot;
  knots->push_back(newk);
  newk->kn_no     = knots->size();
  Parton * newp   = new Parton;
  newk->part      = newp;
  newk->part->SetNumber(newk->kn_no);
  newk->left      = 0;
  newk->right     = 0;
  newk->prev      = 0;
  if (!root) {
    root = newk;
    // msg.Debugging()<<"Knot::NewKnot : "<<newk<<" serves as new root."<<endl;
  }
  return newk;
}

//-----------------------------------------------------------------------
//--------------------------- Resetting the tree ------------------------
//----------------------------------------------------------------------- 


void Tree::Restore(Knot * in) {
  //msg.Debugging()<<"Tree::Restore."<<endl
  //		 <<"   Delete all predecessors of "<<in<<" "<<in->kn_no<<endl;

  if (!(in->prev)) return;
  if (in->prev->left != in)  ResetDaughters(in->prev->left);
  if (in->prev->right != in) ResetDaughters(in->prev->right);
  Restore(in->prev);
  for (Knot_Iterator kit=knots->begin(); kit!=knots->end(); ++kit) {
    if ((*kit) == in->prev) {
      delete (*kit);
      knots->erase(kit);
      return;
    }
  }
}

void Tree::ResetDaughters(Knot * in) {
  if (!(in)) return;
  //msg.Debugging()<<"Tree::ResetDaughters."<<endl
  //		 <<"   Delete daughters of/and "<<in<<" "<<in->kn_no<<endl;
  
  if (in->left)  ResetDaughters(in->left);
  if (in->right) ResetDaughters(in->right);
  for (Knot_Iterator kit=knots->begin(); kit!=knots->end(); ++kit) {
    if ((*kit) == in) {
      delete (in);
      knots->erase(kit);
      return;
    }
  }
}

void Tree::Restore() {
  for (Knot_Iterator kit=knots->begin(); kit!=knots->end(); ++kit) {
    if ( ( (*kit)->part->Info() != 'G') && ( (*kit)->part->Info() != 'H') ) {
      //msg.Debugging()<<"Tree::Restore : Delete knot "<<(*kit)->kn_no
      //	     <<"  "<<(*kit)->part->Flav()<<"  "<<(*kit)->part->Info()<<endl;
      delete (*kit);
      knots->erase(kit);
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

Knot * Tree::GetInitiator() {
  Knot * help;
  help = root;
  while (help->prev) help = help->prev;
  return help;
}

Knot * Tree::GetRoot() { return root; } ;
  


//-----------------------------------------------------------------------
//--------------------------- New frames for the tree -------------------
//----------------------------------------------------------------------- 

void Tree::BoRo(AMATOOLS::Poincare & lorenz) {
  Knot * mo= GetRoot();
  while (mo->prev) mo = mo->prev;
  // msg.Debugging()<<"changing Knot "<<mo->kn_no<<" from "<<mo->part->Momentum()<<endl;
  mo->part->SetMomentum(lorenz*mo->part->Momentum());
  // msg.Debugging()<<"                to "<<mo->part->Momentum()<<endl;
  BoRoDaughters(lorenz,mo);
};

void Tree::BoRoDaughters(AMATOOLS::Poincare & lorenz, Knot * mo) {
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
};

void Tree::BoRo(AMATOOLS::Poincare & lorenz, Knot * mo) {
  // msg.Debugging()<<"changing Knot "<<mo->kn_no<<" from "<<mo->part->Momentum()<<endl;
  mo->part->SetMomentum(lorenz*mo->part->Momentum());
  // msg.Debugging()<<"                to "<<mo->part->Momentum()<<endl;
  BoRoDaughters(lorenz,mo);
};



