//bof
//Version: 1 ADICIC++-0.0/2004/03/12

//Flavtest.C - testing the dipole-flavour structure.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <ioextra>
#include <enumextra>
#include "Dipole.H"


//#define _OUTPUT _OUTPUT





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





class Base {
public:
  virtual void Hello() { cout<<"Base says hello."<<endl;}
};
class Real : public Base {
public:
  void Hello() { cout<<"Hello for Real."<<endl;}
};
class Anti : public Base {
public:
  void Hello() { cout<<"Hello for Anti."<<endl;}
};





int main() {

  cout<<endl;
  cout<<"=============================================================="<<endl;

  Real r;
  r.Hello();
  Anti a;
  a.Hello();

  cout<<endl;

  Base* h; h=&r;
  //h=static_cast<Real*>(h);
  cout<<typeid(h).name()<<endl;
  cout<<typeid(*h).name()<<endl;
  cout<<typeid(r).name()<<endl;
  cout<<typeid(a).name()<<endl;
  h->Hello();

  Base* g; g=new Anti;
  g->Hello();

  cout<<endl;

  cout<<"=============================================================="<<endl;
  cout<<endl;

  //abort();

  //assert(1==2);    //it works

  const int* pint=new int(7);
  cout<<*pint<<endl;
  delete pint;    //this is pretty strange
  cout<<*pint<<endl;
  //*pint=5;
  //cout<<*pint<<endl;

  const Dipole* const pdip=new Dipole;
  delete pdip;    //this is pretty strange

  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"======================================="<<endl;
  cout<<" Testing the dipole-flavour structure. "<<endl;
  cout<<"======================================="<<endl;

  {

    //ParticleInit("./data");    //not necessary due to Dipole_Flavour_Init

    Dipole_Quark_Base* base;
    //Dipole_Quark_Base Base;    //indicates an error as wished
    base=new Dipole_Quark_U;
    const Flavour& fl1=(*base)();
    cout<<"Test flavour: "<<fl1<<" from group: "<<base->Tag()<<endl;
    Dipole_Gluon_G glu;
    const Flavour& fl2=glu();
    cout<<"Test flavour: "<<fl2<<" from group: "<<glu.Tag()<<endl;
    delete base;

    base=new Dipole_Quark_D;
    cout<<"Test flavour: "<<(*base)()<<" "<<base->Tag()<<endl; delete base;
    base=new Dipole_Quark_S;
    cout<<"Test flavour: "<<(*base)()<<" "<<base->Tag()<<endl; delete base;
    base=new Dipole_Quark_C;
    cout<<"Test flavour: "<<(*base)()<<" "<<base->Tag()<<endl; delete base;
    base=new Dipole_Quark_B;
    cout<<"Test flavour: "<<(*base)()<<" "<<base->Tag()<<endl; delete base;
    base=new Dipole_Quark_T;
    cout<<"Test flavour: "<<(*base)()<<" "<<base->Tag()<<endl; delete base;

    Dipole_Antiquark_Base* anba;
    Dipole_Gluon_Base* guba;

    anba=new Dipole_Antiquark_D;
    cout<<"Test flavour: "<<(*anba)()<<" "<<anba->Tag()<<endl; delete anba;
    anba=new Dipole_Antiquark_U;
    cout<<"Test flavour: "<<(*anba)()<<" "<<anba->Tag()<<endl; delete anba;
    anba=new Dipole_Antiquark_S;
    cout<<"Test flavour: "<<(*anba)()<<" "<<anba->Tag()<<endl; delete anba;
    anba=new Dipole_Antiquark_C;
    cout<<"Test flavour: "<<(*anba)()<<" "<<anba->Tag()<<endl; delete anba;
    anba=new Dipole_Antiquark_B;
    cout<<"Test flavour: "<<(*anba)()<<" "<<anba->Tag()<<endl; delete anba;
    anba=new Dipole_Antiquark_T;
    cout<<"Test flavour: "<<(*anba)()<<" "<<anba->Tag()<<endl; delete anba;
    guba=new Dipole_Gluon_G;
    cout<<"Test flavour: "<<(*guba)()<<" "<<guba->Tag()<<endl; delete anba;

    cout<<"============================================================"<<endl;

    Flavour x=Flavour(kf::b).Bar(); cout<<x<<endl;
    cout<<"More testing: "<<info.quark.d()<<" "<<info.antiquark.b()<<endl;
    cout<<"More testing: "<<info.gluon.g()<<" "<<info.antiquark.d()<<endl;
    cout<<"More testing: "<<info.antiquark.s()<<" "<<info.quark.t()<<endl;

    cout<<"============================================================"<<endl;

    Duo duo=Down; bool4::level lev=bool4::two;
    cout<<"duo="<<duo<<" lev="<<lev<<endl;
    double td=7.2431893289;
    cout<<td<<" ==> "<<sform(td)<<" !"<<endl;
    sform.status();
    cout<<td<<" ==> "<<sform75(td)<<endl;
    sform75.status();

  }

  cout<<"======================================"<<endl;
  cout<<" More testing: The dipole structures. "<<endl;
  cout<<"======================================"<<endl;

  {

    cout<<"Global number of dipoles: "<<Dipole::InStore<<endl;
    Dipole D;
    cout<<"Dipole "<<D.Name<<" has status     "<<D.Status()<<endl;
    cout<<"Dipole "<<D.Name<<" is copied from "<<D.CopyOf<<endl;
    cout<<"Dipole "<<D.Name<<" emerged from   "<<D.Source()<<endl;
    cout<<"Dipole "<<D.Name<<" is of type     "<<D.IsType()<<endl;
    cout<<"Pointer handling of Dipole "<<D.Name<<": "
	<<D.PointerHandling()<<endl;
    cout<<"TotNumber="<<D.InStore<<endl;
    cout<<"DipaNumber="<<Dipole_Particle::InStore<<endl;
    Dipole B(D);
    cout<<"Dipole "<<B.Name<<" has status     "<<B.Status()<<endl;
    cout<<"Dipole "<<B.Name<<" is copied from "<<B.CopyOf<<endl;
    cout<<"Dipole "<<B.Name<<" emerged from   "<<B.Source()<<endl;
    cout<<"Dipole "<<B.Name<<" is of type     "<<B.IsType()<<endl;
    cout<<"Pointer handling of Dipole "<<B.Name<<": "
	<<B.PointerHandling()<<endl;
    cout<<"TotNumber="<<Dipole::InStore<<endl;
    cout<<"DipaNumber="<<Dipole_Particle::InStore<<endl;
    Dipole C(B,true);
    cout<<"Dipole "<<C.Name<<" has status     "<<C.Status()<<endl;
    cout<<"Dipole "<<C.Name<<" is copied from "<<C.CopyOf<<endl;
    cout<<"Dipole "<<C.Name<<" emerged from   "<<C.Source()<<endl;
    cout<<"Dipole "<<C.Name<<" is of type     "<<C.IsType()<<endl;
    cout<<"Pointer handling of Dipole "<<C.Name<<": "
	<<C.PointerHandling()<<endl;
    cout<<"BranchNumber="<<Dipole_Particle::InStore<<endl;
    cout<<D<<endl<<B<<endl<<C<<endl;
    D.PrintTowers();
    D.GetTopBranchPointer()->ShowParticle();
    D.GetBotBranchPointer()->ShowParticle();
    B.PrintTowers();
    C.PrintTowers();
    cout<<"Global number of dipoles: "<<Dipole::InStore<<endl;

    //inline const Dipole::Branch*const& Dipole::GetBranch() const {
    //  return (p_top);
    //}
    //warning: returning reference to temporary

    cout<<"============================================================"<<endl;

    cout<<"Global number of dipoles: "<<Dipole::InStore<<endl;
    Dipole::Branch ban;
    Dipole::Antibranch ati;
    Dipole dipt(ban,ati,7,true);
    Dipole dipz(ban,ati,9,false);
    cout<<"Dipole "<<dipt.Name<<" has status     "<<dipt.Status()<<endl;
    cout<<"Dipole "<<dipt.Name<<" is copied from "<<dipt.CopyOf<<endl;
    cout<<"Dipole "<<dipt.Name<<" emerged from   "<<dipt.Source()<<endl;
    cout<<"Dipole "<<dipt.Name<<" is of type     "<<dipt.IsType()<<endl;
    cout<<"Dipole "<<dipt.Name<<" has vector     "<<dipt.TotP()<<endl;
    cout<<"Pointer handling of Dipole "<<dipt.Name<<": "
	<<dipt.PointerHandling()<<endl;
    cout<<"Dipole "<<dipz.Name<<" has status     "<<dipz.Status()<<endl;
    cout<<"Dipole "<<dipz.Name<<" is copied from "<<dipz.CopyOf<<endl;
    cout<<"Dipole "<<dipz.Name<<" emerged from   "<<dipz.Source()<<endl;
    cout<<"Dipole "<<dipz.Name<<" is of type     "<<dipz.IsType()<<endl;
    cout<<"Dipole "<<dipz.Name<<" has mass       "<<dipz.Mass()<<endl;
    cout<<"Dipole "<<dipz.Name<<" has sqrmass    "<<dipz.InvMass()<<endl;
    cout<<"Pointer handling of Dipole "<<dipz.Name<<": "
	<<dipz.PointerHandling()<<endl;
    dipz=dipt;
    cout<<"Dipole "<<dipz.Name<<" has status     "<<dipz.Status()<<endl;
    cout<<"Dipole "<<dipz.Name<<" is copied from "<<dipz.CopyOf<<endl;
    cout<<"Dipole "<<dipz.Name<<" emerged from   "<<dipz.Source()<<endl;
    cout<<"Dipole "<<dipz.Name<<" is of type     "<<dipz.IsType()<<endl;
    cout<<"Dipole "<<dipz.Name<<" has vector     "<<dipz.TotP()<<endl;
    cout<<"Dipole "<<dipz.Name<<" has sqrmass    "<<dipz.InvMass()<<endl;
    cout<<"Pointer handling of Dipole "<<dipz.Name<<": "
	<<dipz.PointerHandling()<<endl;
    cout<<"DipaNumber="<<Dipole_Particle::InStore<<endl;
    cout<<"Global number of dipoles: "<<Dipole::InStore<<endl;

  }

  cout<<"=============================================================="<<endl;

  cout<<endl; cin>>enter; cout<<endl;

  cout<<"======================================"<<endl;
  cout<<"  More: Testing the Dipole_Particle.  "<<endl;
  cout<<"======================================"<<endl;

  {

    Dipole_Particle utest;
    cout<<(&(*utest))<<endl;
    utest.Quarkize(info.quark.t);
    cout<<(&(*utest))<<endl;
    utest.Antiquarkize(info.antiquark.u);
    cout<<(&(*utest))<<endl;
    utest.Gluonize();
    cout<<(&(*utest))<<endl;
    utest.Antiquarkize(info.antiquark.b);
    cout<<(&(*utest))<<endl;

    cout<<"------------------------------------------------------------"<<endl;

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Particle no;
    //Dipole::Branch b0(no);    //as wished exits

    Dipole::Branch b1;    //d-quark
    b1.WhatIsIt();
    b1.SetMomentum(pl);
    cout<<"Branch "<<b1.Name<<" has type      "<<b1.OrgType()<<endl;
    cout<<"Branch "<<b1.Name<<" has tag       "<<b1.Tag()<<endl;
    cout<<"Branch "<<b1.Name<<" has flavour   "<<b1.Flav()<<endl;
    cout<<"Branch "<<b1.Name<<" has vector    "<<b1.Momentum()<<endl;
    cout<<"Branch "<<b1.Name<<" has particle:"<<endl;
    cout<<(&(*b1))<<endl;
    cout<<"Branch "<<b1.Name<<": Show the particle:";
    b1.ShowParticle();
    b1.ShowDipoles();
    //b1.Antiquarkize(info.antiquark.u);    //crashes due to d-quark
    //b1.Quarkize(info.quark.s);
    b1.Gluonize();
    cout<<"Branch "<<b1.Name<<" has type      "<<b1.OrgType()<<endl;
    cout<<"Branch "<<b1.Name<<" has tag       "<<b1.Tag()<<endl;
    cout<<"Branch "<<b1.Name<<" has flavour   "<<b1.Flav()<<endl;
    cout<<"Branch "<<b1.Name<<" has vector    "<<b1.Momentum()<<endl;
    cout<<"Branch "<<b1.Name<<" has particle:"<<endl;
    cout<<(&(*b1))<<endl;

    Dipole::Glubranch g1;
    cout<<"Branch "<<g1.Name<<" has type      "<<g1.OrgType()<<endl;
    cout<<"Branch "<<g1.Name<<" has tag       "<<g1.Tag()<<endl;
    cout<<"Branch "<<g1.Name<<" has flavour   "<<g1.Flav()<<endl<<endl;
    const Particle& newpa=g1.Quarkize(info.quark.u);
    cout<<"Branch "<<g1.Name<<" has type      "<<g1.OrgType()<<endl;
    cout<<"Branch "<<g1.Name<<" has tag       "<<g1.Tag()<<endl;
    cout<<"Branch "<<g1.Name<<" has flavour   "<<g1.Flav()<<endl;
    cout<<"Branch "<<g1.Name<<" has vector    "<<g1.Momentum()<<endl;
    cout<<"Branch "<<g1.Name<<" has particle:"<<endl;
    cout<<(&newpa)<<endl;
    g1.ShowDipoles();

    //Dipole::Branch b7(info.antiquark.b,pr);    //exits as wished

    Dipole::Branch b2(info.quark.d,pr);
    cout<<"Branch "<<b2.Name<<" has type      "<<b2.OrgType()<<endl;
    cout<<"Branch "<<b2.Name<<" has tag       "<<b2.Tag()<<endl;
    cout<<"Branch "<<b2.Name<<" has flavour   "<<b2.Flav()<<endl;
    cout<<"Branch "<<b2.Name<<" has vector    "<<b2.Momentum()<<endl;
    cout<<"Branch "<<b2.Name<<" has particle:"<<endl;
    cout<<(&(*b2))<<endl;
    cout<<"Branch "<<b2.Name<<": Show the particle:";
    b2.ShowParticle();
    b2.ShowDipoles();

    Dipole::Branch b3(b2);
    cout<<"Branch "<<b3.Name<<" has type      "<<b3.OrgType()<<endl;
    cout<<"Branch "<<b3.Name<<" has tag       "<<b3.Tag()<<endl;
    cout<<"Branch "<<b3.Name<<" has flavour   "<<b3.Flav()<<endl;
    cout<<"Branch "<<b3.Name<<" has vector    "<<b3.Momentum()<<endl;
    cout<<"Branch "<<b3.Name<<" has particle:"<<endl;
    cout<<(&(*b3))<<endl;
    b3.ShowDipoles();

    b3=b1;
    cout<<"Branch "<<b3.Name<<" has type      "<<b3.OrgType()<<endl;
    cout<<"Branch "<<b3.Name<<" has tag       "<<b3.Tag()<<endl;
    cout<<"Branch "<<b3.Name<<" has flavour   "<<b3.Flav()<<endl;
    cout<<"Branch "<<b3.Name<<" has vector    "<<b3.Momentum()<<endl;
    cout<<"Branch "<<b3.Name<<" has particle:"<<endl;
    cout<<(&(*b3))<<endl;

    //Dipole::Antibranch a1(b1);    //gives error as wished - no mixing allowed
    //b1==b2;
    //Dipole::Antibranch a0(info.quark.u,pr);

    Dipole::Antibranch a1(info.antiquark.u,pr);
    cout<<"Antibranch "<<a1.Name<<" has type      "<<a1.OrgType()<<endl;
    cout<<"Antibranch "<<a1.Name<<" has tag       "<<a1.Tag()<<endl;
    cout<<"Antibranch "<<a1.Name<<" has flavour   "<<a1.Flav()<<endl;
    cout<<"Antibranch "<<a1.Name<<" has vector    "<<a1.Momentum()<<endl;
    cout<<"Antibranch "<<a1.Name<<" has particle:"<<endl;
    cout<<(&(*a1))<<endl;
    cout<<"Antibranch "<<a1.Name<<": Show the particle:";
    a1.ShowParticle();
    a1.ShowDipoles();
    Dipole::Antibranch a2(a1);
    cout<<"Antibranch "<<a2.Name<<" has type      "<<a2.OrgType()<<endl;
    cout<<"Antibranch "<<a2.Name<<" has tag       "<<a2.Tag()<<endl;
    cout<<"Antibranch "<<a2.Name<<" has flavour   "<<a2.Flav()<<endl;
    cout<<"Antibranch "<<a2.Name<<" has vector    "<<a2.Momentum()<<endl;
    cout<<"Antibranch "<<a2.Name<<" has particle:"<<endl;
    cout<<(&(*a2))<<endl;
    //a2.Quarkize(info.quark.s);
    //a2.Antiquarkize(info.antiquark.s);
    //cout<<"Antibranch "<<a2.Name<<" has type      "<<a2.OrgType()<<endl;
    //cout<<"Antibranch "<<a2.Name<<" has tag       "<<a2.Tag()<<endl;
    //cout<<"Antibranch "<<a2.Name<<" has flavour   "<<a2.Flav()<<endl<<endl;
    a2.Gluonize();
    cout<<"Antibranch "<<a2.Name<<" has type      "<<a2.OrgType()<<endl;
    cout<<"Antibranch "<<a2.Name<<" has tag       "<<a2.Tag()<<endl;
    cout<<"Antibranch "<<a2.Name<<" has flavour   "<<a2.Flav()<<endl;
    cout<<"Antibranch "<<a2.Name<<" has vector    "<<a2.Momentum()<<endl;
    cout<<"Antibranch "<<a2.Name<<" has particle:"<<endl;
    cout<<(&(*a2))<<endl;

    Dipole Dip(b1,a2,99);
    cout<<Dip<<endl;
    Dip.PrintTowers();
    Dip.Print();
    Dipole Dup(b2,a2,77);
    cout<<Dup<<endl;
    Dup.PrintTowers();
    b1.ShowDipoles();
    b2.ShowDipoles();
    b3.ShowDipoles();
    a1.ShowDipoles();
    a2.ShowDipoles();
    b2.SetMomentum(pl); cout<<Dup<<endl;
    //a2.Antiquarkize(info.quark.d);    //no matching function
    Dip.SetStatus()=Blocked;
    //Dip.GetConstAntibranch().SetMomentum(pl);    //gives error as wished
    b1.SetMomentum(pr);
    cout<<Dip.TotP()<<endl;
    Dup.SetSource()=777;
    cout<<Dip<<endl<<Dup<<endl;
    //Dip==Dup;    //as wished exits
    cout<<"Dipole "<<Dip.Name<<"'s Branch momentum="
	<<Dip.GetTopBranchPointer()->Momentum()<<endl;
    cout<<"Dipole "<<Dup.Name<<"'s Antibranch flavour="
    	<<Dup.GetBotBranchPointer()->Flav()<<endl;
    Dip=Dup;
    cout<<Dip<<endl<<Dup<<endl;
    b1.ShowDipoles();
    b2.ShowDipoles();
    b3.ShowDipoles();
    a1.ShowDipoles();
    a2.ShowDipoles();
    cout<<"DipoleTotalNumber="<<Dip.InStore<<endl;
    cout<<"BranchTotalNumber="<<Dipole::Branch::InStore<<endl;

    cout<<"------------------------------------------------------------"<<endl;

    cout<<utest.OrgType()<<" : "
	<<g1.OrgType()<<" : "
	<<b1.OrgType()<<" : "
	<<b2.OrgType()<<" : "
	<<b3.OrgType()<<" : "
	<<a1.OrgType()<<" : "
	<<a2.OrgType()<<endl;

    utest.WhatIsIt();
    g1.WhatIsIt();
    b1.WhatIsIt();
    b2.WhatIsIt();
    b3.WhatIsIt();
    a1.WhatIsIt();
    a2.WhatIsIt();

    utest.ShowParticle();
    g1.ShowParticle();
    b1.ShowParticle();
    b2.ShowParticle();
    b3.ShowParticle();
    a1.ShowParticle();
    a2.ShowParticle();

    utest.ShowDipoles();
    g1.ShowDipoles();
    b1.ShowDipoles();
    b2.ShowDipoles();
    b3.ShowDipoles();
    a1.ShowDipoles();
    a2.ShowDipoles();

    //Fill test
    //=========

    //Dipole fa(utest,g1,41);    //no compile
    //Dipole fb(g1,g1,42);    //no run : single gluon
    //Dipole fc(a1,g1,43);    //no compile
    //Dipole fd(g1,b1,44);    //no compile
    //Dipole fe(b1,g1,45);    //no run : invalid dipole
    Dipole ff(b2,a2,46,true);
    Dipole fg(g1,a1,47);
    //Dipole fh(b3,g1,48);    //no run : invalid dipole
    //Dipole fi(a2,a1,49);    //no compile
    //Dipole fj(b3,b2,50);    //no compile

    ff.GetTopBranchPointer()->WhatIsIt();
    ff.GetBotBranchPointer()->WhatIsIt();
    fg.GetTopBranchPointer()->WhatIsIt();
    fg.GetBotBranchPointer()->WhatIsIt();

    cout<<ff<<endl<<fg<<endl;

    ff.RenewBranch(a1);
    ff.RenewBranch(b3);
    fg.RenewBranch(b3);
    ff.RenewBranch(true,g1);
    //fg.RenewBranch(false,g1);    //no run : invalid dipole

    cout<<ff<<endl<<fg<<endl;

    //Construction tests
    //==================

    //Dipole::Branch b4(utest);
    Dipole::Branch b4(b3);
    //Dipole::Branch b4(a2);
    //Dipole::Branch b4(g1);

    //Dipole::Antibranch a3(utest);
    //Dipole::Antibranch a3(b3);
    //Dipole::Antibranch a3(g1);
    Dipole::Antibranch a3(a2);

    //Dipole::Glubranch g2(utest);
    //Dipole::Glubranch g2(b3);
    //Dipole::Glubranch g2(a2);
    Dipole::Glubranch g2(g1);

    Dipole::Branch b5(info.quark.s,pr);
    //Dipole::Branch b6(info.antiquark.d,pl);
    //Dipole::Branch b7(info.gluon.g,pr);

    //Dipole::Antibranch a4(info.quark.s,pr);
    Dipole::Antibranch a4(info.antiquark.d,pl);
    //Dipole::Antibranch a6(info.gluon.g,pr);

    //Dipole::Glubranch g3(info.quark.s,pr);
    //Dipole::Glubranch g4(info.antiquark.d,pl);
    //Dipole::Glubranch g5(info.gluon.g,pr);    //not yet implemented
    Dipole::Glubranch g3(pr);

    //b4=utest;
    b4=b5;
    //b4=a4;
    //b4=g3;
    //a3=utest;
    //a3=b5;
    a3=a4;
    //a3=g3;
    //g2=utest;
    //g2=b5;
    //g2=a4;
    g2=g3;

    //Cast test
    //=========

    Dipole_Particle* P=&g3;
    Dipole_Particle* p=(Dipole::Antibranch*)P;    //has no effect
    cout<<typeid(P).name()<<endl;
    cout<<typeid(*P).name()<<endl;
    cout<<typeid(p).name()<<endl;
    cout<<typeid(*p).name()<<endl;

    //p=&(Dipole::Branch*)a3;    //no matching function
    p=(Dipole::Branch*)&a3;
    //Dipole D1(*p,a3);    //no matching function
    //Dipole D1(*(Dipole::Branch*)&a3,a3);    //this is bad, ext UpdTyp helps
    //Dipole D1(*static_cast<Dipole::Branch*>(&a3),a3);    //no compile
    //Dipole D1(*dynamic_cast<Dipole::Branch*>(&a3),a3);    //error new op
    p=(Dipole::Antibranch*)&a3;
    cout<<typeid(p).name()<<endl;
    cout<<typeid(*p).name()<<endl;
    //Dipole D2(g1,*p);    //no matching function
    Dipole D2(g1,*(Dipole::Antibranch*)p);

    a4.WhatIsIt();
    p=&a4;
    P->WhatIsIt();
    p->WhatIsIt();
    p=P;
    p->WhatIsIt();
    cout<<" ----"<<endl;
    b1.WhatIsIt();
    P=&b1;
    p=P;
    P=&a3;
    P->WhatIsIt();
    p->WhatIsIt();
    //Dipole D3(*p,a1);    //no matching function
    Dipole D3(*static_cast<Dipole::Branch*>(p),a1);
    //Dipole D4(*static_cast<Dipole::Branch*>(P),a1);

    Dipole::Branch* pb=static_cast<Dipole::Branch*>(P);
    cout<<pb->OrgType()<<endl;
    pb->WhatIsIt();
    //Dipole D4(*pb,a1);    //cast protector

    cout<<D3<<endl;//<<D4<<endl;
    D3.GetTopBranchPointer()->WhatIsIt();
    D3.GetBotBranchPointer()->WhatIsIt();

  }
 
  cout<<"=============================================================="<<endl;

  cout<<endl; cin>>enter; cout<<endl;

  cout<<"======================================"<<endl;
  cout<<"  More: Testing the PointerHandling.  "<<endl;
  cout<<"======================================"<<endl;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);

    Dipole D1;
    Dipole D2(D1);
    Dipole::Branch b1(info.quark.u,pl);
    Dipole::Antibranch a1(info.antiquark.u,pr);
    Dipole D3(b1,a1,77);
    Dipole D4(D3);
    cout<<D1<<endl<<D2<<endl<<D3<<endl<<D4<<endl;
    D1.PrintTowers(); D2.PrintTowers(); D3.PrintTowers(); D4.PrintTowers();
    cout<<endl;

    D3=D2;
    D1=D4;
    cout<<D1<<endl<<D2<<endl<<D3<<endl<<D4<<endl;
    D1.PrintTowers(); D2.PrintTowers(); D3.PrintTowers(); D4.PrintTowers();
    cout<<endl;

  }

  cout<<"=============================================================="<<endl;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole D1;
    Dipole::Branch b0(info.quark.s,pl);
    Dipole::Antibranch a0(info.antiquark.t,pl);
    Dipole D2(b0,a0,33);
    cout<<D1<<endl<<D2<<endl;
    D1.PrintTowers(); D2.PrintTowers();
    D1=D2;
    cout<<endl;

  }

  cout<<"=============================================================="<<endl;

  //goto n1;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole D1;
    Dipole D2(D1);
    Dipole::Branch b0(info.quark.u,pl);
    Dipole::Antibranch a0(info.antiquark.t,pl);
    Dipole D3(b0,a0,33);

    Dipole::Particle_Pointer b;
    Dipole::Particle_Pointer* bp;
    //bp=new Dipole::Particle_Pointer(b);
    //  error: `ADICIC::Dipole::Branch_Pointer::Branch_Pointer(const
    //  ADICIC::Dipole::Branch_Pointer&)' is private
    //bp=new Dipole::Particle_Pointer(b0);
    //  no matching function for call to `
    //  ADICIC::Dipole::Branch_Pointer::Branch_Pointer(
    //  ADICIC::Dipole_Particle<ADICIC::_Real>&)'
    //  ../Structure/Dipole.inl.hh:357: error: candidates are:
    //  ADICIC::Dipole::Branch_Pointer::Branch_Pointer()
    //  ../Structure/Dipole.inl.hh:382: error:
    //  ADICIC::Dipole::Branch_Pointer::Branch_Pointer(const
    //  ADICIC::Dipole::Branch_Pointer&)

    Dipole::Particle_Pointer test;
    //test->ShowParticle();    //gives defined abort
    cout<<"Branch_Pointer test 1 ? 0=="<<bool(test)<<endl;
    Dipole::Particle_Pointer test2;
    cout<<"Branch_Pointer test 2 ? 1=="<<(test==test2)<<endl;
    cout<<"Branch_Pointer test 3 ? 0=="
	<<(test==D1.GetTopBranchPointer())<<endl;
    cout<<"Branch_Pointer test 4 ? 1=="
	<<(D1.GetTopBranchPointer()==D2.GetTopBranchPointer())<<endl;
    cout<<"Branch_Pointer test 5 ? 0=="
	<<(D1.GetTopBranchPointer()==D3.GetTopBranchPointer())<<endl;

    cout<<(D1.GetTopBranchPointer()==D3.GetBotBranchPointer())<<endl;

    //Dipole::Particle_Pointer tt=D2.GetTopBranchPointer();
    //Dipole::Particle_Pointer tl=D2.GetBotBranchPointer();
    //  error: `ADICIC::Dipole::Branch_Pointer::Branch_Pointer(const
    //  ADICIC::Dipole::Branch_Pointer&)' is private

    cout<<"Branch_Pointer test 6"<<endl;
    //delete D2.GetTopBranchPointer();    //does not work :)
    Dipole::Particle_Pointer* ptest;
    //ptest=&(D2.GetTopBranchPointer);
    //  error: ISO C++ forbids taking the address of a bound member
    //  function to form a pointer to member function.  Say `&ADICIC::Dipole::
    //  GetBranchPointer'
    //  error: cannot convert `ADICIC::Dipole::Branch_Pointer
    //  (ADICIC::Dipole::*)()' to `ADICIC::Dipole::Branch_Pointer*' in assign.
    //(*ptest)->ShowParticle();    //gives segm. fault as expected
    //delete ptest;    //attention
    //delete (*ptest);    //does not work :)

    cout<<"Branch_Pointer test 7"<<endl;
    //delete D2.GetTopBranchPointer().operator->();
    //  Bug: TopBranch physically belongs to a Dipole!
    //Dipole::Branch* pB=D2.GetTopBranchPointer().operator->();
    //  error: invalid conversion from `ADICIC::Dipole_Particle*' to `
    //  ADICIC::Dipole_Branch*'
    Dipole_Particle* pB=D2.GetTopBranchPointer().operator->();
    Dipole_Particle* pA=D2.GetBotBranchPointer().operator->();
    //  who really wants to do that?

    D2.GetTopBranchPointer()->ShowParticle();
    D2.GetBotBranchPointer()->ShowParticle();
    D2.GetBotBranchPointer()->Antiquarkize(info.antiquark.s);
    D2.GetTopBranchPointer()->Quarkize(info.quark.s);
    D2.GetTopBranchPointer()->ShowParticle();
    D2.GetBotBranchPointer()->ShowParticle();

    cout<<"\n+++++++++++++++++"<<endl;

    D2=D3;
    D2.GetTopBranchPointer()->ShowParticle();
    pB->ShowParticle();
    pB->SetMomentum(pr);
    //delete pB;    //Bug: TopBranch physically belongs to a Dipole!
    //delete pA;    //Bug: Antibranch physically belongs to a Dipole!

    cout<<D1<<endl<<D2<<endl<<D3<<endl;
    D1.PrintTowers(); D2.PrintTowers(); D3.PrintTowers();

    //////////////////////////
    //  Pain in the ass(1)  //
    //////////////////////////
    cout<<"Branch_Pointer test 8"<<endl;
    double fake;
    double* pfake=&fake;
    //delete pfake;    //results in free(): invalid pointer ...! or segm. fault
    pB=D2.GetTopBranchPointer().operator->();
    //delete pB;    //results in a segmentation fault, since
    //  a Branch which hasn't been created with a new operator is destructed.
    cout<<"Attention."<<endl;

    cout<<endl;

    //non-smart pointer version:
    //gives dangling pointers since b0 and a0 are destroyed before D2
    //results in an segmentation fault while destroying D2 afterwards
    //more explicit: while the application of RemoveDipoleFromTowers
    //first two ShowParticle() give the same

    //new: smart pointer version:
    //pointer test is always up-to-date - Branch_Pointer
    //no segmentation fault
    //Dipoles seemingly never have dangling pointers

    //PLEASE AVOID IT:
    //================
    //the only crux is Dipole::Branch* pB=D2.GetBranchPointer().operator->();
    //or delete D2.GetBranchPointer().operator->();
    //but this is really malicious

  }
 
  cout<<"=============================================================="<<endl;

 n1:

  //goto n2;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole D1;
    Dipole D2(D1);
    D2.GetTopBranchPointer()->ShowParticle();
    Dipole_Particle* pB=D2.GetTopBranchPointer().operator->();
    Dipole_Particle* pC=NULL;
    pB->WhatIsIt();
    cout<<pB<<endl;    //(*)
    {
      Dipole::Branch b2(info.quark.u,pl);
      Dipole::Antibranch a2(info.antiquark.t,pl);
      Dipole D3(b2,a2,33);
      D2=D3;
      pC=D2.GetTopBranchPointer().operator->();
      pC->ShowParticle();
      pB->ShowParticle();
      //delete pC;
      //  proceeds destructor correctly
      //  but fails the destruction due to the fact that b2 was not created by
      //  the new operator: free(): invalid pointer 0x7fbfffdd10!
      //  pain in the ass(1)

      const int& name=D1.GetBotBranchPointer()->Name;
      const Flavour& fla=D1.GetBotBranchPointer()->Flav();
      cout<<"Antibranch before="<<name<<endl;
      cout<<"Flavour before="<<fla<<endl;
      D1.GetBotBranchPointer()->WhatIsIt();
      const Particle& pa1=
	D1.GetBotBranchPointer()->Antiquarkize(info.antiquark.u);
      D1.GetBotBranchPointer()->WhatIsIt();
      const Particle& pa2=
	D1.GetBotBranchPointer()->Quarkize(info.quark.t);
      D1.GetBotBranchPointer()->WhatIsIt();
      cout<<"Flavour before="<<fla<<endl;
      cout<<"Flavour before="<<pa1.Flav()<<endl;
      cout<<"Flavour before="<<pa2.Flav()<<endl;
      cout<<D1<<endl;

      D1=D3;

      //////////////////////////
      //  Pain in the ass(3)  //
      //////////////////////////
      cout<<"Antibranch after="<<name<<endl;    //here it shows strange behav.
      cout<<"Flavour after="<<fla<<endl;
      cout<<"Flavour after="<<pa1.Flav()<<endl;
      cout<<"Flavour after="<<pa2.Flav()<<endl;
      //Be careful with these things.

      cout<<"New Flavour="<<(D1.GetBotBranchPointer()->Flav())<<endl;

      D2.GetTopBranchPointer()->ShowParticle();
      pB->ShowParticle();    //Be careful! expect CRASH-pain in the ass(2) (**)
      cout<<(D2.GetTopBranchPointer().operator->())<<endl;
      cout<<(D1.GetTopBranchPointer().operator->())<<endl;
      cout<<pB<<endl;
      //  here luckily the last two addresses are the same, see also (*)

      cout<<D1<<endl<<D2<<endl<<D3<<endl;
      D1.PrintTowers(); D2.PrintTowers(); D3.PrintTowers();
    }
    D2.GetTopBranchPointer()->ShowParticle();    //no crash

    //////////////////////////
    //  Pain in the ass(2)  //
    //////////////////////////
    pB->ShowParticle();    //expected CRASH, compare (**) and (*)
    //pC->ShowParticle();    //expected CRASH - yes

    cout<<D1<<endl<<D2<<endl;    //no crash
    D1.PrintTowers(); D2.PrintTowers();    //no crash

    cout<<endl;

  }

  cout<<"=============================================================="<<endl;

 n2:

  cout<<"=============================================================="<<endl;
  cout<<endl;

}





//eof
