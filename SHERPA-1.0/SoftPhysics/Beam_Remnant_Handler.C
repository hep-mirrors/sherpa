#include "Beam_Remnant_Handler.H"

#include "Message.H"
#include "Random.H"
#include "Vector.H"
#include "Data_Read.H"

using namespace SHERPA;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


Beam_Remnant_Handler::Beam_Remnant_Handler(std::string _dir,std::string _file,
					   PDF::ISR_Handler * _isr,
					   BEAM::Beam_Spectra_Handler * _beam) :
  m_dir(_dir), m_file(_file), p_isr(_isr), p_beam(_beam),
  p_constituents(NULL), p_numberofconstituents(NULL)
{
  if (p_isr->Flav(0).IsHadron() || p_isr->Flav(1).IsHadron()) {
    p_constituents         = new Flavour *[2];
    p_numberofconstituents = new int[2];
    for (int i=0;i<2;i++) {
      p_numberofconstituents[i] = 0;
      if (p_isr->Flav(i).IsHadron()) 
	p_numberofconstituents[i] = Constituents(p_isr->Flav(i),p_constituents[i]);
      else p_constituents[i] = NULL;
    }
    m_fill = 1;
  }
  else m_fill = 0;

  m_q2min = 1.;
}

int Beam_Remnant_Handler::Constituents(Flavour _had,Flavour *& _flavs) {
  int hadint = (_had.Kfcode() - _had.Kfcode()/10000)/10;
  
  if ((hadint > 100) && (hadint < 1000)) {
    _flavs    = new Flavour[3];
    _flavs[0] = Flavour(kf::code(hadint)/100);
    _flavs[1] = Flavour(kf::code((hadint-(hadint/100)*100)/10));
    _flavs[2] = Flavour(kf::code(hadint-(hadint/10)*10));
    if (_had.IsAnti()) {
      for(int i=0;i<3;i++) _flavs[i]=_flavs[i].Bar();
    }
    msg.Tracking()<<"Beam_Remnant_Handler::Constituents for "<<_had<<" is baryon :"
		  <<hadint<<std::endl<<"   "<<_flavs[0]<<", "<<_flavs[1]<<", "<<_flavs[2]<<std::endl;
    return 3;
  }
  if ((hadint > 10) && (hadint < 100)) {
    _flavs    = new Flavour[2];
    _flavs[0] = Flavour(kf::code(hadint)/10);
    _flavs[1] = Flavour(kf::code(hadint-(hadint/10)*10));
    if (_had.IsAnti()) {
      for(int i=0;i<2;i++) _flavs[i]=_flavs[i].Bar();
    }
    msg.Tracking()<<"Beam_Remnant_Handler::Constituents for "<<_had<<" is meson :"
		  <<hadint<<std::endl<<"   "<<_flavs[0]<<", "<<_flavs[1]<<std::endl;
    return 2;
  }
  
  msg.Error()<<"Beam_Remnant_Handler::Constituents :"
	     <<"No idea how to handle this case: "<<_had<<". Abort."<<std::endl;
  abort();
}      

Beam_Remnant_Handler::~Beam_Remnant_Handler() {  
  if (p_constituents) {
    for (int i=0;i<2;i++) {
      if (p_constituents[i]) delete [] p_constituents[i];
    }
    delete [] p_constituents; p_constituents = NULL;
  }
  
  if (p_numberofconstituents) delete [] p_numberofconstituents;
}


bool Beam_Remnant_Handler::FillBunchBlobs(Blob_List * _bloblist,Parton_List * _partonlist)
{
  Blob_Iterator endblob = _bloblist->end(); 
  Blob * blob;
  bool match;
  bool flag = 0;
  int pos;
  for (int i=0;i<2;i++) {
    for (Blob_Iterator biter = _bloblist->begin();biter != endblob;++biter) {
      match = 1;
      if (m_fill) match = ((*biter)->Type()==std::string("Beam Remnant"));
      if ((*biter)->Status()==1 && (*biter)->Beam()==i && match) {
	blob = new APHYTOOLS::Blob();
	blob->SetId(_bloblist->size());
	blob->SetType(std::string("Bunch"));
	blob->SetBeam(i);
	blob->AddToOutPartons((*biter)->InParton(0));
	(*biter)->InParton(0)->SetProductionBlob(blob);
	if ((*biter)->InParton(0)->Flav()==p_beam->GetBeam(i)->Flav() &&
	    (*biter)->InParton(0)->E()==p_beam->GetBeam(i)->Momentum()[0]) {
	  Parton * p = new Parton((*biter)->InParton(0));
	  p->SetDecayBlob(blob);
	  p->SetProductionBlob(NULL);
	  blob->AddToInPartons(p);
	}
	else {
	  blob->AddToInPartons(new Parton(-1,p_beam->GetBeam(i)->Flav(),p_beam->GetBeam(i)->Momentum()));
	  blob->InParton(0)->SetDecayBlob(blob);
	  blob->InParton(0)->SetStatus(2);
	  (*biter)->InParton(0)->SetProductionBlob(blob);
	}
	_bloblist->insert(_bloblist->begin(),blob);
	flag = 1;
      }
    }
  }
  return flag;
}

bool Beam_Remnant_Handler::FillBeamBlobs(Blob_List * _bloblist,Parton_List * _partonlist)
{
  if (!m_fill) return 0;
  Blob_Iterator endblob = _bloblist->end(); 
  Blob * blob;
  bool okay = 0;
  int pos;
  for (int i=0;i<2;i++) {
    bool flag=1;
    for (Blob_Iterator biter = _bloblist->begin();biter != endblob;++biter) {
      pos = (*biter)->Type().find(string("IS"));
      if ((*biter)->Beam()==i && (*biter)->Status()==1 && 
	  (p_numberofconstituents[i]>0) && pos>-1 ) {
	if (flag) {
	  blob = new APHYTOOLS::Blob();
	  blob->SetId(_bloblist->size());
	  blob->SetType(std::string("Beam Remnant"));
	  blob->SetBeam(i);
	  blob->SetStatus(1);
	  blob->AddToInPartons(new Parton(-1,p_isr->Flav(i),p_beam->GetBeam(i)->OutMomentum()));
	  blob->InParton(0)->SetStatus(2);
	  _bloblist->insert(_bloblist->begin(),blob);
	  flag = 0;
	  okay = 1;
	}
	blob->AddToOutPartons((*biter)->InParton(0));
	// Similar things are needed for photons, leptons, etc. !!!!
	if (p_isr->Flav(i).IsHadron()) okay == okay && FillHadron(blob,i,_partonlist);
      }
    }
  }
  return okay;
}
		   
		   
bool Beam_Remnant_Handler::FillHadron(Blob * blob,int _beam,Parton_List * pl)
{
  if (blob->NOutP()!=1) {
    msg.Error()<<"Error in Beam_Remnant_Handler::FillHadron("<<blob->NOutP()<<")"<<std::endl
	       <<"   This case is not implemented yet ! Abort."<<std::endl;
    abort();
  }

  Vec4D vec1,vec2;
  Flavour fl,difl;
  int di[2];
  int pos;

  Parton * part = blob->OutParton(0);
  Vec4D   vec   = blob->InParton(0)->Momentum() + (-1.)*(part->Momentum());


  /*  
      Parton is one of the constituents.
      
      Colour flow:
      
      ------------- part (assumed quark, i.e hadron = baryon)   = (a,0)
      |
      |
      |
      ------------- newpart1 (assumed diquark)                  = (0,a)

  */

  int number;
  for (int i=0;i<3;i++) {
    if (part->Flav() == p_constituents[blob->Beam()][i]) {
      pos = 0;
      for (int j=0;j<3;j++) {
	if (i!=j) {
	  di[pos]  = p_constituents[blob->Beam()][j].Kfcode();
	  pos++;
	}
      }
      if (di[0] != di[1]) {
	if (ran.Get()<0.25) difl = Flavour(kf::code(di[0]*1000+di[1]*100+1));
	else difl = Flavour(kf::code(di[0]*1000+di[1]*100+3));
      }
      else difl = Flavour(kf::code(di[0]*1100+3));
      if (p_constituents[blob->Beam()][0].IsAnti()) difl = difl.Bar();

      Parton * newpart1 = new Parton(-1,difl,vec); 
      if (pl) number = pl->size();
      else    number = int(newpart1);
      newpart1->SetNumber(number);
      newpart1->SetFlow(1,part->GetFlow(2));
      newpart1->SetFlow(2,part->GetFlow(1));
      newpart1->SetProductionBlob(blob);
      newpart1->SetInfo('F');
      blob->AddToOutPartons(newpart1);
      if (pl) pl->push_back(newpart1);
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
  if ((part->Flav()).IsGluon()) {
    int single = int(ran.Get()*3.); 
    fl  = p_constituents[blob->Beam()][single];
    pos = 0;
    for (int i=0;i<3;i++) {
      if (i!=single) {
	di[pos]  = p_constituents[blob->Beam()][i].Kfcode();
	pos++;
      }
    }
    if (di[0] != di[1]) {
      if (ran.Get()<0.25) difl = Flavour(kf::code(di[0]*1000+di[1]*100+1));
      else difl = Flavour(kf::code(di[0]*1000+di[1]*100+3));
    }
    else difl = Flavour(kf::code(di[0]*1100+3));
    if (p_constituents[blob->Beam()][0].IsAnti()) difl = difl.Bar();
    
    vec1 = GetX_Lund(fl,difl,vec[0]) * vec;
    vec2 = vec + (-1.)*vec1;

    Parton * newpart1 = new Parton(-1,fl,vec1); 
    if (pl) number = pl->size();
    else    number = int(newpart1);
    newpart1->SetNumber(number);
    newpart1->SetProductionBlob(blob);
    newpart1->SetInfo('F');
    blob->AddToOutPartons(newpart1);
    if (pl) pl->push_back(newpart1);

    Parton * newpart2 = new Parton(-1,difl,vec2); 
    if (pl) number = pl->size();
    else    number = int(newpart2);
    newpart2->SetNumber(number);
    newpart2->SetProductionBlob(blob);
    newpart2->SetInfo('F');
    blob->AddToOutPartons(newpart2);
    if (pl) pl->push_back(newpart2);

    if (fl.IsAnti()) {
      newpart1->SetFlow(1,0);
      newpart1->SetFlow(2,part->GetFlow(1));
      newpart2->SetFlow(1,part->GetFlow(2));
      newpart2->SetFlow(2,0);
    }
    else {
      newpart1->SetFlow(1,part->GetFlow(2));
      newpart1->SetFlow(2,0);
      newpart2->SetFlow(1,0);
      newpart2->SetFlow(2,part->GetFlow(1));
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

  if ((part->Flav()).IsQuark()) {
    int single = int(ran.Get()*3.); 
    fl  = p_constituents[blob->Beam()][single];
    pos = 0;
    for (int i=0;i<3;i++) {
      if (i!=single) {
	di[pos]  = p_constituents[blob->Beam()][i].Kfcode();
	pos++;
      }
    }
    if (di[0] != di[1]) {
      if (ran.Get()<0.25) difl = Flavour(kf::code(di[0]*1000+di[1]*100+1));
      else difl = Flavour(kf::code(di[0]*1000+di[1]*100+3));
    }
    else difl = Flavour(kf::code(di[0]*1100+3));
    if (p_constituents[blob->Beam()][0].IsAnti()) difl = difl.Bar();

    vec1 = GetX_Lund(fl,difl,vec[0]) * vec;
    vec2 = GetX_Lund((part->Flav()).Bar(),difl,vec[0]) * vec;

    Parton * newpart1 = new Parton(-1,fl,vec1); 
    if (pl) number = pl->size();
    else    number = int(newpart1);
    newpart1->SetNumber(number);
    newpart1->SetProductionBlob(blob);
    newpart1->SetInfo('F');
    blob->AddToOutPartons(newpart1);
    if (pl) pl->push_back(newpart1);

    Parton * newpart2 = new Parton(-1,(part->Flav()).Bar(),vec2); 
    if (pl) number = pl->size();
    else    number = int(newpart2);
    newpart2->SetNumber(number);
    newpart2->SetProductionBlob(blob);
    newpart2->SetInfo('F');
    blob->AddToOutPartons(newpart2);
    if (pl) pl->push_back(newpart2);

    Parton * newpart3 = new Parton(-1,difl,vec + (-1.)*(vec1+vec2)); 
    if (pl) number = pl->size();
    else    number = int(newpart3);
    newpart3->SetNumber(number);
    newpart3->SetProductionBlob(blob);
    newpart3->SetInfo('F');
    blob->AddToOutPartons(newpart3);
    if (pl) pl->push_back(newpart3);

    if ( (fl.IsAnti() && !((part->Flav()).IsAnti()) ) ||
	 (!(fl.IsAnti()) && (part->Flav()).IsAnti() ) ) {
      /* 
	 newpart1, the quark from the hadron connects with the outgoing
	 hard parton part, since one is a quark and the other is an antiquark.
	 In this case, the antiflavour to part connects with the diquark
	 to form a baryon (baryonic cluster).
      */
      
      newpart1->SetFlow(1,part->GetFlow(2));
      newpart1->SetFlow(2,part->GetFlow(1));

      if (fl.IsAnti()) {
	newpart2->SetFlow(1,0);
	newpart2->SetFlow(2,-1);
	newpart3->SetFlow(1,newpart2->GetFlow(2));
	newpart3->SetFlow(2,0);
      }
      else {
	newpart2->SetFlow(1,-1);
	newpart2->SetFlow(2,0);
	newpart3->SetFlow(1,0);
	newpart3->SetFlow(2,newpart2->GetFlow(1));
      }
    }
    else {
      /* 
	 newpart3, the diquark from the hadron connects with the outgoing
	 hard parton part, since one is a (anti-) diquark and the other is 
	 a (anti-) quark. In this case, the antiflavour to part connects with the 
	 (anti-) quark from the hadron to form a meson (mesonic cluster).
      */
      newpart3->SetFlow(1,part->GetFlow(2));
      newpart3->SetFlow(2,part->GetFlow(1));

      if (fl.IsAnti()) {
	newpart2->SetFlow(1,-1);
	newpart2->SetFlow(2,0);
	newpart1->SetFlow(1,0);
	newpart1->SetFlow(2,newpart2->GetFlow(1));
      }
      else {
	newpart2->SetFlow(1,0);
	newpart2->SetFlow(2,-1);
	newpart1->SetFlow(1,newpart2->GetFlow(2));
	newpart1->SetFlow(2,0);
      }
    }
    return 1;
  }
  return 0;
}

double Beam_Remnant_Handler::GetX_Lund(Flavour f1,Flavour f2,double E) {
  double rn, cut, mass; 
  if (f1.IsQuark() && f2.IsDiQuark()) {
    mass = 0.3;
    cut  = 2.*mass/E;
    for (;;) {
      rn = ran.Get();
      // Check the boundaries.
      if ((rn > cut) && (rn < (1.-cut))) {
	if (pow(1.-rn,3)/pow(rn*rn+cut*cut,0.25) > 1./cut*ran.Get()) return rn;
      }
    }
  }
  msg.Error()<<"Beam_Remnant_Handler::GetX_Lund called with Flavours "
	     <<f1<<", "<<f2<<std::endl
	     <<"   Don't know, how to handle this. Return 0."<<std::endl;

  return 0.;
}

double Beam_Remnant_Handler::GetX_PDF(Flavour f1,Flavour f2,double E,int beam) {
  double cut1,cut2,x;
  PDF::PDF_Base * pdf = p_isr->PDF(beam);
  if (f1.IsQuark() && f2.IsDiQuark()) {
    cut1 = f1.PSMass()/E;
    cut2 = 1.-f2.PSMass()/E;
    for (;;) {
      x = ran.Get();
      if (x>cut1 && x<cut2) {
	pdf->Calculate(x,m_q2min);
	if (pdf->GetXPDF(f1)/x>ran.Get()) return x;
      }
    }
  }
  msg.Error()<<"Beam_Remnant_Handler::GetX called with Flavours "
	     <<f1<<", "<<f2<<std::endl
	     <<"   Don't know, how to handle this. Return 0."<<std::endl;
  return 0.;
}






