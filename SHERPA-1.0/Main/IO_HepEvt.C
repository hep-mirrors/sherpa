#include "IO_HepEvt.H"
#include "Message.H"

using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

extern "C" {
  void f2parton_(int&, int*, double*);
  void f2hepevt_(int&, int*, int*, int*, int*, double*, double*);
  
  void pyhepc_(int&);
  void pylist_(int&);
};



// sort Colored Event according to "Color Flow" information
void IO_HepEvt::MakeColorChain() {
  if (mypl.size()==0) return;
  
  // check if colored at all / find a color seed
  int colored=0;
  Parton * seed =0;
  Parton_List::iterator seed_pit=mypl.begin();
  Parton_List::iterator start_pit=mypl.begin();


  for (;;) {
    // look for a seed
    colored=0;
    seed   =0;
    for (Parton_List::iterator pit=start_pit;pit!=mypl.end();++pit) {
      if ((*pit)->Flav().Strong()) {
	colored=1;
	if ((*pit)->GetFlow(1)!=0  && (*pit)->GetFlow(2)==0) {
	  // found a seed
	  seed     = (*pit);
	  seed_pit = pit;
	}
      }
    }

    if (!colored) break;

    if (!seed) {
      msg.Error()<<"WARNING in IO_HepEvt, no singlet found!"<<endl;
      break;
    }

    // start with seed
    if (seed_pit!=start_pit) {
      // swap if seed is not the begin of the list
      Parton * help  = (*(start_pit));
      (*(start_pit))= (*(seed_pit));
      (*(seed_pit))= help;
      seed_pit=start_pit;
      //      cout<<"+"<<flush<<endl;
    }

    // sort colors
    for (Parton_List::iterator qit=start_pit;qit!=mypl.end();) {
      //      cout<<(*qit)->Flav()<<" ("<<(*qit)->GetFlow(1)<<","<<(*qit)->GetFlow(2)<<")"<<endl;
      if ((*qit)->GetFlow(2)==(*seed_pit)->GetFlow(1)) {
	++start_pit;
	// found a color partner
	if (qit!=start_pit) {
	  // swap if not the begin of remaining list
	  Parton * help  = (*(start_pit));
	  (*(start_pit))=(*(qit));
	  (*(qit))= help;
	}
	seed_pit = start_pit;
	qit      = start_pit;
      }
      else {
	++qit;
      }
    }

    // check if finished; (perhaps a second color singlet)
    if (start_pit==--mypl.end()) break;
  }
}

void IO_HepEvt::Blob2HepEvt(Blob * blob, int mode) {
  // mode ==  0  add outpartons to  hepevt and transfer to FORTRAN
  // mode ==  1  add outpartons to internal list
  // mode == 10  add inpartons to  hepevt and transfer to FORTRAN
  // mode == 11  add inpartons to internal list
  // mode == 20  add in- and outpartons to  hepevt and transfer to FORTRAN
  // mode == 21  add in- and outpartons to internal list

  if (blob) {

    if (mode>=10) {
      for (int i=0;i<blob->NInP();++i) {
	Parton2HepEvt(blob->InParton(i),1);
      }
    }

    if (mode<=10 || mode>=20 ) {
      for (int i=0;i<blob->NOutP();++i) {
	Parton2HepEvt(blob->OutParton(i),1);
      }
    }

    for (int i=0;i<nhep;i++) {
      msg.Tracking()<<"  "<<i<<" : "<<std::endl
	       <<"  "<<idhep[i]<<" / ("
	       <<phep[5*i+3]<<","<<phep[5*i+0]<<","
	       <<phep[5*i+1]<<","<<phep[5*i+2]<<")"<<std::endl<<std::endl;
    }
  }

  if (mode%10==0) {
    MakeColorChain();
    Parton2HepEvt(&mypl);
    mypl.clear();

    f2hepevt_(nhep, isthep, idhep, jmohep, jdahep, phep, vhep);

    // test output
    /*
    int dummy=2;
    pyhepc_(dummy);
    dummy=1;
    pylist_(dummy);
    */


  }
}


void IO_HepEvt::Blob2HepEvt(int blob, int mode) {
  // mode ==  0/10/20  transfer internal list to  hepevt (FORTRAN)

  if (blob) {
    msg.Error()<<" UNEXPECTED void IO_HepEvt::Blob2HepEvt(int blob, int mode) call "<<endl;
  } 

  if (mode%10==0) {
    MakeColorChain();
    Parton2HepEvt(&mypl);
    mypl.clear();


    f2hepevt_(nhep, isthep, idhep, jmohep, jdahep, phep, vhep);

    // test output
    /*
    int dummy=2;
    pyhepc_(dummy);
    dummy=1;
    pylist_(dummy);
    */
  }
}

void IO_HepEvt::Blob2HepEvt(Blob_List * blobs, int mode) {
  for (Blob_List::const_iterator bit=blobs->begin(); bit!=blobs->end();++bit) {
    Blob2HepEvt((*bit),1);
  }

  if (mode==0) {
    MakeColorChain();
    Parton2HepEvt(&mypl);
    mypl.clear();
  }

}

void IO_HepEvt::Parton2HepEvt(Parton * parton, int mode) {
  // add parton to internal list

  mypl.push_back(parton);
}

void IO_HepEvt::Parton2HepEvt(Parton_List * pl, int mode) {
  if (mode==0) {
    nhep=0;

  }

  int ic=0;
   for (Parton_List::const_iterator pit=pl->begin(); pit!=pl->end();++pit) {
    Parton * parton=(*pit);

    if (nhep>=maxentries) {
      --nhep;
      msg.Error()<<" ERROR :in void IO_HepEvt::Parton2HepEvt()"<<endl
		 <<"   nehp>=maxentries "<<endl;
    }

    idhep[nhep] = parton->Flav().HepEvt();
        
    for (short int j=1; j<4; ++j) phep[(j-1)+nhep*5] = parton->Momentum()[j];
    phep[3+nhep*5] = parton->Momentum()[0];
    double pabs = (parton->Momentum()).Abs2();
    if (pabs<0) phep[4+nhep*5] = 0.;
    else phep[4+nhep*5] = sqrt(pabs);

    for (short int j=1; j<4; ++j) {
      vhep[(j-1)+nhep*4] = parton->XProd()[j];
      vhep[3+nhep*4] = parton->XProd()[j];
    }

    isthep[nhep] = 1;
    jmohep[nhep*2] = 0;
    jmohep[1+nhep*2] = 0;
  
    jdahep[nhep*2] = 0;
    jdahep[1+nhep*2] = 0;

    ++nhep;
  }

  if (mode==0) {

    // transfer to FORTRAN HEPEVT
    f2hepevt_(nhep, isthep, idhep, jmohep, jdahep, phep, vhep);

    // test output
    /*
    int dummy=2;
    pyhepc_(dummy);
    dummy=1;
    pylist_(dummy);
    */
  }
}



IO_HepEvt::~IO_HepEvt() {
  if (phep)   delete phep;
  if (vhep)   delete vhep;
  if (jmohep) delete jmohep;
  if (jdahep) delete jdahep;
  if (isthep) delete isthep;
  if (idhep)  delete idhep;
}

IO_HepEvt::IO_HepEvt()
{
  nhep = 0;
  phep   = new double[5*maxentries];
  vhep   = new double[4*maxentries];
  jmohep = new int[2*maxentries];
  jdahep = new int[2*maxentries];
  isthep = new int[maxentries];
  idhep  = new int[maxentries];
}
