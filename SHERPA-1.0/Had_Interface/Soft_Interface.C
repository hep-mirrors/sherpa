#include "Soft_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"

using namespace MOCAIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


Soft_Interface::Soft_Interface() {
  mypartons = new Parton_List;
  myblobs   = new Blob_List;

  Data_Read dr(rpa.GetPath()+std::string("/Run.dat"));
  double lund_a     = dr.GetValue<double>("LUND_A",0.4);
  double lund_b     = dr.GetValue<double>("LUND_B",0.85);
  double lund_sigma = dr.GetValue<double>("LUND_SIGMA",0.36);
  msg.Out()<<" LUND_A = "<<lund_a<<endl;
  msg.Out()<<" LUND_B = "<<lund_b<<endl;
  msg.Out()<<" LUND_SIGMA = "<<lund_sigma<<endl;

  lund      = new Lund_Fortran_Interface(lund_a,lund_b,lund_sigma);
  if (rpa.gen.Beam1().ishadron() || rpa.gen.Beam2().ishadron()) {
    constituents = new Flavour *[2];
    for (int i=0;i<2;i++) constituents[i] = new Flavour[3];
    n_const      = new int[2];
    n_const[0]   = Constituents(rpa.gen.Beam1(),constituents[0]); 
    n_const[1]   = Constituents(rpa.gen.Beam2(),constituents[1]); 
  }
  else {
    constituents = 0;n_const = 0; 
  }
}
   
Soft_Interface::~Soft_Interface() {
  msg.Tracking()<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  msg.Tracking()<<"In Soft_Interface::~Soft_Interface :"<<std::endl;

  if (!(mypartons->empty())) {
    mypartons->erase(mypartons->begin(),mypartons->end());
  }
  if (mypartons)    delete mypartons;
  msg.Tracking()<<"   Deleted mypartons."<<std::endl;

  if (constituents) delete constituents;
  msg.Tracking()<<"   Deleted constituents."<<std::endl;

  if (n_const)      delete n_const;
  msg.Tracking()<<"   Deleted n_const."<<std::endl;

  msg.Tracking()<<"Out Soft_Interface::~Soft_Interface :"<<std::endl;
  msg.Tracking()<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
}; 




void Soft_Interface::EmptyMyLists() {
  if (!(mypartons->empty())) {
    mypartons->erase(mypartons->begin(),mypartons->end());
    msg.Debugging()<<"   Deleted partons of internal parton list."<<std::endl;
  }
  if (!(myblobs->empty())) {
    myblobs->erase(myblobs->begin(),myblobs->end());
    msg.Debugging()<<"   Deleted partons of internal parton list."<<std::endl;
  }
}

bool Soft_Interface::PerformFragmentation(APHYTOOLS::Blob_List * bl,
					  APHYTOOLS::Parton_List * pl) {

  EmptyMyLists();
  if (!ExtractSinglets(bl,pl)) return 0;
  /*
    myblob->BoostInCMS();
    lund->Hadronize(pl,myblob);
    myblob->BoostInLab();
  */
  for (Blob_Iterator biter = myblobs->begin();biter != myblobs->end();++biter) {
    if ( ((*biter)->Type() == std::string("Fragmentation")) && 
	 ((*biter)->Status() != 0) ) {
      (*biter)->BoostInCMS();
      lund->Hadronize(pl,(*biter));
      (*biter)->BoostInLab();
    }
  }
}

bool Soft_Interface::ExtractSinglets(Blob_List * bl,Parton_List * pl) {
  msg.Debugging()<<"Soft_Interface::ExtractSinglets :"<<std::endl<<(*pl)<<std::endl;
  bool found = 0;

  Blob   * newb=0;
  bool   use_one_blob=1;

  for (Parton_Iterator piter = pl->begin();piter != pl->end();++piter) {
    if ( ( ((*piter)->info() == 'F') || ((*piter)->info() == 'H') ) &&
	 ( (*piter)->status() == 1) ) {
      if ( ((((*piter)->flav().isquark()) && (!(*piter)->flav().isanti())) ||
	    (((*piter)->flav().isdiquark()) && (*piter)->flav().isanti())) &&
	   ((*piter)->flow(1) > 0)) {
	msg.Debugging()<<"   Start new singlet chain."<<std::endl
		       <<"   for parton : "<<(*piter)<<std::endl;

	if (use_one_blob==0 || newb==0) {
	  newb = new Blob();
	  newb->SetId(bl->size());
	  newb->SetStatus(1);
	  newb->SetType(std::string("Fragmentation"));
	  myblobs->push_back(newb);
	  bl->push_back(newb);
	}

	(*piter)->set_status(2);
	(*piter)->set_dec(newb);

	mypartons->push_back((*piter));
	newb->AddToInPartons(new Parton(*piter));

	found = 1;

	if (!(FindConnected(pl,(*piter),newb))) {
	  msg.Error()<<"Soft_Interface::ExtractSinglets :"
		     <<"Could not find connected parton for quark:"<<std::endl
		     <<(*piter)<<std::endl;
	  return 0;
	}
      }
    }
  }
  if (!(found)) {
    msg.Debugging()<<"Soft_Interface::ExtractSinglets :"
		   <<"Could not find any colour singlet in this process."<<std::endl;
    return 1;
  }
  return 1;
}; 


bool Soft_Interface::FindConnected(Parton_List * pl,Parton * compare,Blob * blob) {
  msg.Debugging()<<"Soft_Interface::FindConnected :"<<std::endl;

  for (Parton_Iterator piter = pl->begin();piter != pl->end();++piter) {
    if ((*piter) == compare) continue;
    if ( ((*piter)->info() != 'F') && ((*piter)->info() != 'H')) continue;
    if ((*piter)->status() != 1) continue; 
    if ((*piter)->flow(2) == compare->flow(1)) {
      msg.Debugging()<<"Colour match for "<<compare->flow(1)<<std::endl;
      Parton * newp = (*piter);
      newp->set_status(2);
      newp->set_dec(blob);

      mypartons->push_back(newp);
      blob->AddToInPartons(new Parton(newp));

      if ((*piter)->flow(1) ==0) {
	if (((*piter)->flav().isquark()) && ((*piter)->flav().isanti())) {
	  msg.Debugging()<<"Closed singlet list with an antiquark."<<std::endl;
	  return 1;
	}
	if (((*piter)->flav().isdiquark()) && (!((*piter)->flav().isanti()))) {
	  msg.Debugging()<<"Closed singlet list with a diquark."<<std::endl;
	  return 1;
	}
      }
      else { return FindConnected(pl,(*piter),blob); }
    }
  }
  msg.Error()<<"Soft_Interface::FindConnected : No closed singlet line !"<<std::endl;
  return 0;
}


/*
  bool Soft_Interface::ExtractSinglets(Blob_List * bl) {
  if (!(mypartons->empty())) {
  mypartons->erase(mypartons->begin(),mypartons->end());
  }
  }; 
*/




bool Soft_Interface::HadronsToPartons(Blob_List * bl,Parton_List * pl) {
  Flavour  flav;
  vec4d    vec;
  Parton * beam;
  Blob   * blob;

  Blob_Iterator endblob = bl->end();
  for (Blob_Iterator biter = bl->begin();biter != endblob;++biter) {
    if (((*biter)->Beam() == 0) || ((*biter)->Beam() == 1)) {
      Parton * inpart = (*biter)->InParton(0);

      blob = new APHYTOOLS::Blob();
      blob->SetId(bl->size());
      blob->SetType(std::string("Beam Remnant"));
      bl->push_back(blob);

      if ((*biter)->Beam()==0) {
	flav = rpa.gen.Beam1();
	vec  = vec4d(rpa.gen.Ecms()/2.,0.,0.,
		     sqrt(sqr(rpa.gen.Ecms()/2.)-sqr(flav.mass())));
	beam = new Parton(-1,flav,vec);
	blob->SetBeam(0);
      }
      if ((*biter)->Beam()==1) {
	flav = rpa.gen.Beam2();
	vec  = vec4d(rpa.gen.Ecms()/2.,0.,0.,
		     -sqrt(sqr(rpa.gen.Ecms()/2.)-sqr(flav.mass())));
	beam = new Parton(-1,flav,vec);
	blob->SetBeam(1);
      }      
      blob->AddToInPartons(beam);
      FillHadron(pl,blob,inpart);
    }
  }
}

int Soft_Interface::Constituents(Flavour had,Flavour * flavs) {
  int hadint = (had.kfcode() - had.kfcode()/10000)/10;

  if ((hadint > 100) && (hadint < 1000)) {
    flavs[0] = Flavour(kf::code(hadint)/100);
    flavs[1] = Flavour(kf::code((hadint-(hadint/100)*100)/10));
    flavs[2] = Flavour(kf::code(hadint-(hadint/10)*10));
    if (had.isanti()) {
      for(int i=0;i<3;i++) flavs[i]=flavs[i].bar();
    }
    msg.Events()<<"Soft_Interface::Constituents for "<<had<<" is baryon :"
		<<hadint<<std::endl<<"   "<<flavs[0]<<", "<<flavs[1]<<", "<<flavs[2]<<std::endl;

    return 3;
  }
  if ((hadint > 10) && (hadint < 100)) {
    flavs[0] = Flavour(kf::code(hadint)/10);
    flavs[1] = Flavour(kf::code(hadint-(hadint/10)*10));
    flavs[2] = Flavour(kf::none);
    if (had.isanti()) {
      for(int i=0;i<2;i++) flavs[i]=flavs[i].bar();
    }
    msg.Events()<<"Soft_Interface::Constituents for "<<had<<" is meson :"
		<<hadint<<std::endl<<"   "<<flavs[0]<<", "<<flavs[1]<<std::endl;

    return 2;
  }

  msg.Error()<<"Soft_Interface::HadronsToPartons :"
	     <<"No idea how to handle this case: "<<had<<std::endl;

  return 0;
}      

bool Soft_Interface::FillHadron(Parton_List * pl,Blob * blob,Parton * part) {
  msg.Tracking()<<"Soft_Interface::FillHadron : "
		<<part->flav()<<" in "<<blob->InParton(0)->flav()<<std::endl;


  Parton * newpart = part;
  blob->AddToOutPartons(newpart);

  vec4d vec1,vec2;
  Flavour fl,difl;
  int di[2];
  int pos;

  vec4d   vec = blob->InParton(0)->momentum() + (-1.)*(part->momentum());


  /*  
      Parton is one of the constituents.
      
      Colour flow:
      
      ------------- part (assumed quark, i.e hadron = baryon)   = (a,0)
      |
      |
      |
      ------------- newpart1 (assumed diquark)                  = (0,a)

  */
  for (int i=0;i<3;i++) {
    if (part->flav() == constituents[blob->Beam()][i]) {
      msg.Tracking()<<"   Parton is one of the constituents "
		    <<blob->InParton(0)->flav()<<std::endl;

      pos = 0;
      for (int j=0;j<3;j++) {
	if (i!=j) {
	  di[pos]  = constituents[blob->Beam()][j].kfcode();
	  pos++;
	}
      }
      if (di[0] != di[1]) {
	if (Ran.get()<0.25) difl = Flavour(kf::code(di[0]*1000+di[1]*100+1));
	else difl = Flavour(kf::code(di[0]*1000+di[1]*100+3));
      }
      else difl = Flavour(kf::code(di[0]*1100+3));
      if (constituents[blob->Beam()][0].isanti()) difl = difl.bar();

      Parton * newpart1 = new Parton(pl->size(),difl,vec); 
      newpart1->set_flow(1,part->flow(2));
      newpart1->set_flow(2,part->flow(1));
      newpart1->set_prod(blob);
      newpart1->set_info('F');
      blob->AddToOutPartons(newpart1);
      pl->push_back(newpart1);
      return 1;
    }
  }

  /*  
      Parton is a gluon
      
      Colour flow:
      
      ------------- newpart (assumed quark)                     = (a,0)
      |
      |
      ------------- 
                    part (gluon)                                = (a,b)
      ------------- 
      |
      ------------- newpart2 (assumed diquark)                  = (0,b)

  */
  if ((part->flav()).isgluon()) {
    int single = int(Ran.get()*3.); 
    fl  = constituents[blob->Beam()][single];
    pos = 0;
    for (int i=0;i<3;i++) {
      if (i!=single) {
	di[pos]  = constituents[blob->Beam()][i].kfcode();
	pos++;
      }
    }
    if (di[0] != di[1]) {
      if (Ran.get()<0.25) difl = Flavour(kf::code(di[0]*1000+di[1]*100+1));
      else difl = Flavour(kf::code(di[0]*1000+di[1]*100+3));
    }
    else difl = Flavour(kf::code(di[0]*1100+3));
    if (constituents[blob->Beam()][0].isanti()) difl = difl.bar();

    msg.Tracking()<<"   Parton is a gluon."<<std::endl<<"   Split "
		  <<blob->InParton(0)->flav()<<" into "<<fl<<" and "<<difl<<std::endl;
    
    vec1 = GetX(fl,difl,vec[0]) * vec;
    vec2 = vec + (-1.)*vec1;

    Parton * newpart1 = new Parton(pl->size(),fl,vec1); 
    newpart1->set_prod(blob);
    newpart1->set_info('F');
    blob->AddToOutPartons(newpart1);
    pl->push_back(newpart1);

    Parton * newpart2 = new Parton(pl->size(),difl,vec2); 
    newpart2->set_prod(blob);
    newpart2->set_info('F');
    blob->AddToOutPartons(newpart2);
    pl->push_back(newpart2);

    if (fl.isanti()) {
      newpart1->set_flow(1,0);
      newpart1->set_flow(2,part->flow(1));
      newpart2->set_flow(1,part->flow(2));
      newpart2->set_flow(2,0);
    }
    else {
      newpart1->set_flow(1,part->flow(2));
      newpart1->set_flow(2,0);
      newpart2->set_flow(1,0);
      newpart2->set_flow(2,part->flow(1));
    }
    return 1;
  }

  /*  
      Parton is a (sea-) quark
      
      Colour flow:
      
      ------------- newpart3 (assumed diquark)                   = (0,a)
      |
      |
      ------------- part     (quark Q)                          = (a,0)
                    
      ------------- newpart2 (anti quark Qbar)                  = (0,b)
      |
      ------------- newpart1 (assumed quark)                    = (b,0)

      Similarly for antiquarks Qbar as part. Just flip role of Q and Qbar, i.e.
      role of newpart1 and newpart3.
  */

  if ((part->flav()).isquark()) {
    int single = int(Ran.get()*3.); 
    fl  = constituents[blob->Beam()][single];
    pos = 0;
    for (int i=0;i<3;i++) {
      if (i!=single) {
	di[pos]  = constituents[blob->Beam()][i].kfcode();
	pos++;
      }
    }
    if (di[0] != di[1]) {
      if (Ran.get()<0.25) difl = Flavour(kf::code(di[0]*1000+di[1]*100+1));
      else difl = Flavour(kf::code(di[0]*1000+di[1]*100+3));
    }
    else difl = Flavour(kf::code(di[0]*1100+3));
    if (constituents[blob->Beam()][0].isanti()) difl = difl.bar();

    msg.Tracking()<<"   Parton is a quark."<<std::endl
		  <<"   Split "<<blob->InParton(0)->flav()
		  <<" into "<<fl<<" and "<<difl<<std::endl;
    
    vec1 = GetX(fl,difl,vec[0]) * vec;
    vec2 = GetX((part->flav()).bar(),difl,vec[0]) * vec;

    Parton * newpart1 = new Parton(pl->size(),fl,vec1); 
    newpart1->set_prod(blob);
    newpart1->set_info('F');
    blob->AddToOutPartons(newpart1);
    pl->push_back(newpart1);

    Parton * newpart2 = new Parton(pl->size(),(part->flav()).bar(),vec2); 
    newpart2->set_prod(blob);
    newpart2->set_info('F');
    blob->AddToOutPartons(newpart2);
    pl->push_back(newpart2);

    Parton * newpart3 = new Parton(pl->size(),difl,vec + (-1.)*(vec1+vec2)); 
    newpart3->set_prod(blob);
    newpart3->set_info('F');
    blob->AddToOutPartons(newpart3);
    pl->push_back(newpart3);

    if ( (fl.isanti() && !((part->flav()).isanti()) ) ||
	 (!(fl.isanti()) && (part->flav()).isanti() ) ) {
      /* 
	 newpart1, the quark from the hadron connects with the outgoing
	 hard parton part, since one is a quark and the other is an antiquark.
	 In this case, the antiflavour to part connects with the diquark
	 to form a baryon (baryonic cluster).
      */
      
      newpart1->set_flow(1,part->flow(2));
      newpart1->set_flow(2,part->flow(1));

      if (fl.isanti()) {
	newpart2->set_flow(1,0);
	newpart2->set_flow(2,-1);
	newpart3->set_flow(1,newpart2->flow(2));
	newpart3->set_flow(2,0);
      }
      else {
	newpart2->set_flow(1,-1);
	newpart2->set_flow(2,0);
	newpart3->set_flow(1,0);
	newpart3->set_flow(2,newpart2->flow(1));
      }
    }
    else {
      /* 
	 newpart3, the diquark from the hadron connects with the outgoing
	 hard parton part, since one is a (anti-) diquark and the other is 
	 a (anti-) quark. In this case, the antiflavour to part connects with the 
	 (anti-) quark from the hadron to form a meson (mesonic cluster).
      */
      newpart3->set_flow(1,part->flow(2));
      newpart3->set_flow(2,part->flow(1));

      if (fl.isanti()) {
	newpart2->set_flow(1,-1);
	newpart2->set_flow(2,0);
	newpart1->set_flow(1,0);
	newpart1->set_flow(2,newpart2->flow(1));
      }
      else {
	newpart2->set_flow(1,0);
	newpart2->set_flow(2,-1);
	newpart1->set_flow(1,newpart2->flow(2));
	newpart1->set_flow(2,0);
      }
    }
    return 1;
  }
};

double Soft_Interface::GetX(Flavour f1,Flavour f2,double E) {
  double ran, cut, mass; 
  if (f1.isquark() && f2.isdiquark()) {
    mass = 0.3;
    cut  = 2.*mass/E;
    for (;;) {
      ran = Ran.get();
      // Check the boundaries.
      if ((ran > cut) && (ran < (1.-cut))) {
	if (pow(1.-ran,3)/sqrt(ran*ran+cut*cut) > 1./cut*Ran.get()) return ran;
      }
    }
  }
  msg.Error()<<"Soft_Interface::GetX called with Flavours "
	     <<f1<<", "<<f2<<std::endl
	     <<"   Don't know, how to handle this. Return 0."<<std::endl;

  return 0.;
}






