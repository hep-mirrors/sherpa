/*  Flavour.C: Definitions for Flavour class
 *  A Class to define particle properties.
 */

#include <iostream>

#include "Flavour.H"
#include "MathTools.H"
#include "Message.H"
#include "MyStrStream.H"

namespace ATOOLS {
  Part_Info particles[MAX_PARTICLES];
  Kf_To_Int kf_table;
}

using namespace ATOOLS;
using namespace std;

int Kf_To_Int::is_initialised=0;

void Kf_To_Int::Init()
{
  for(anz=0;anz<MAX_PARTICLES;++anz) {
    kf_tab[anz] = particles[anz].kfc;
    if(kf_tab[anz]==kf::none) break;
  }

  if(anz==MAX_PARTICLES) {
    std::cerr<<"ERROR in Kf_To_Int::kftab(): Too many particle types !"<<std::endl;
    exit(1);
  }
  is_initialised = 1;
}

kf::code Kf_To_Int::FromInt(int code)
{
  return kf_tab[code];
}

int Kf_To_Int::ToInt(kf::code kfc)
{
  for(int i=0;i<anz+1;i++) { if (kf_tab[i]==kfc) return i; }
  std::cerr<<"ERROR in Kf_To_Int::to_int(): Particle type unknown ! "
	   <<kfc<<" in "<<anz<<std::endl;
  return -1;
}

kf::code Kf_To_Int::FromString(std::string st)
{
  for(int i=0;i<anz;i++) {
    if (std::string(particles[i].n)==st) return kf_tab[i];
  }

  return kf::none;
}




Part_Info::Part_Info(kf::code kfc_, 
		     double mass, double width, int icharge, int isoweak, 
		     bool strong, int spin, bool Majorana, 
		     bool Take, bool stable,bool massive,
		     char* name,int _masssign) : 
  kfc(kfc_), m(mass), w(width), yuk(mass), iq(icharge), isow(isoweak), sp(spin), masssign(_masssign),  
  str(strong), Maj(Majorana), on(Take), stbl(stable), msv(massive),hadron(0), n(name), group(0) { 
  n  = new char[strlen(name)+1];
  strcpy(n,name);
}

Part_Info::Part_Info(kf::code kfc_, 
		     double mass, double width, int icharge, int isoweak, 
		     int spin, bool Take, bool stable,
		     char* name) : 
  kfc(kfc_), m(mass), w(width), yuk(0.), iq(icharge), isow(isoweak), sp(spin), masssign(1), 
  str(0), Maj(0), on(Take), stbl(stable), msv(1), hadron(1), n(name), group(0) { 
  n  = new char[strlen(name)+1];
  strcpy(n,name);
}

void Part_Info::Add(const Flavour _fl) {
  if (_fl.Size()==1) {
    flavs.push_back(new Flavour(_fl));
    return;
  }
  for (int i=0;i<_fl.Size();i++) {
    flavs.push_back(new Flavour(_fl[i]));
  }
}

const int Part_Info::Size() {
  if (!group) return 1;
  return flavs.size();
}

const bool Part_Info::Includes(Flavour _fl) {
  bool isin;
  for (int j=0;j<_fl.Size();j++) {
    isin = 0;
    for (int i=0;i<flavs.size();i++) {
      if (flavs[i]->Kfcode() == _fl[j].Kfcode()) { isin = 1; break; }
    }
    if (isin==0) return 0;
  }
  return 1;
}

Flavour Part_Info::operator[](const int i) {
  if (!group) return Flavour(kfc);
  if ((i<0) || (i>flavs.size())) {
    std::cerr<<"Error in Part_Info::operator[]("<<i<<") . Delimiter out of bounds."<<std::endl;
    return Flavour(kfc);
  }
  return (*flavs[i]);
}

void Part_Info::SetGroup() { group = 1; }




Flavour Flavour::operator[](int idx) const  { 
  if (!(particles[ Index() ].group)) return (*this);
  if (anti) return particles[ Index() ][idx].Bar();
  return particles[ Index() ][idx];
} 



Flavour Flavour::Bar() { 
  Flavour flbar = *this;
  
  if (flbar==Flavour(kf::photon))   return flbar;
  if (flbar==Flavour(kf::Z))        return flbar;
  if (flbar==Flavour(kf::gluon))    return flbar;
  if (flbar==Flavour(kf::h))        return flbar;
  if (flbar==Flavour(kf::h0))       return flbar;
  if (flbar==Flavour(kf::H0))       return flbar;
  if (flbar==Flavour(kf::A0))       return flbar;
  if (flbar==Flavour(kf::graviton)) return flbar;
  if (flbar==Flavour(kf::gscalar))  return flbar;
  if (flbar==Flavour(kf::none))     return flbar;
  if (Majorana())                   return flbar;

  flbar.anti = (anti)?0:1;
  return flbar;
}

Flavour::Flavour(int kf) { 
  anti = (kf<0)? 1:0;
  kfc  = kf_table.FromInt(abs(kf));
}


int Flavour::Ctq() {
  int code;
  switch(kfc) {
  case kf::u    : code=1;break;
  case kf::d    : code=2;break;
  case kf::s    : code=3;break;
  case kf::c    : code=4;break;
  case kf::b    : code=5;break;
  case kf::t    : code=6;break;
  case kf::gluon: code=0;break;
  default       : code=7;
  }
  if(anti) code= -code;
  return code+6;
}

void Flavour::FromCtq(int code) {
  anti = 0;
  switch(abs(code)) {
  case 0 : kfc=kf::t;break;
  case 1 : kfc=kf::b;break;
  case 2 : kfc=kf::c;break;
  case 3 : kfc=kf::s;break;
  case 4 : kfc=kf::d;break;
  case 5 : kfc=kf::u;break;
  case 6 : kfc=kf::gluon;break;
  default: kfc=kf::none;
  }
}

int Flavour::HepEvt() {
  if (IsLepton() || IsQuark() || IsHadron()) return (anti)? -Kfcode():Kfcode();
  if (IsDiQuark())                           return (anti)? -Kfcode():Kfcode();
  if (IsGluon())                             return 21;
  if (IsPhoton())                            return 22;
  if (kfc==kf::Z)                            return 23;
  if (kfc==kf::W)                            return (anti)? 24:-24;
  if ((kfc==kf::h) || (kfc==kf::h0))         return 25;

  if (kfc==kf::H0)                           return 35;
  if (kfc==kf::A0)                           return 36;
  if (kfc==kf::Hmin)                         return (anti)? 37:-37;
  if (kfc==kf::graviton)                     return 39;
  if (kfc==kf::gscalar)                      return 89;

  if (IsSquark() || IsSlepton() || IsSneutrino() || IsIno()) {
    
    int pdgnum = -1;
    
    if (kfc==kf::sDownL)                     pdgnum = 1000001;
    if (kfc==kf::sUpL)                       pdgnum = 1000002;
    if (kfc==kf::sStrangeL)                  pdgnum = 1000003;
    if (kfc==kf::sCharmL)                    pdgnum = 1000004;
    if (kfc==kf::sBottom1)                   pdgnum = 1000005;
    if (kfc==kf::sTop1)                      pdgnum = 1000006;
    if (kfc==kf::sElectronL)                 pdgnum = 1000011;
    if (kfc==kf::sNu1)                       pdgnum = 1000012;
    if (kfc==kf::sMuL)                       pdgnum = 1000013;
    if (kfc==kf::sNu2)                       pdgnum = 1000014;
    if (kfc==kf::sTau1)                      pdgnum = 1000015;
    if (kfc==kf::sNu3)                       pdgnum = 1000016;
    
    if (kfc==kf::sDownR)                     pdgnum = 2000001;
    if (kfc==kf::sUpR)                       pdgnum = 2000002;
    if (kfc==kf::sStrangeR)                  pdgnum = 2000003;
    if (kfc==kf::sCharmR)                    pdgnum = 2000004;
    if (kfc==kf::sBottom2)                   pdgnum = 2000005;
    if (kfc==kf::sTop2)                      pdgnum = 2000006;
    if (kfc==kf::sElectronR)                 pdgnum = 2000011;
    if (kfc==kf::sMuR)                       pdgnum = 2000013;
    if (kfc==kf::sTau2)                      pdgnum = 2000015;
  
    if (kfc==kf::Chargino1)                  pdgnum = 1000024; 
    if (kfc==kf::Chargino2)                  pdgnum = 1000037; 
    
    if (pdgnum!=-1) return (anti) ? pdgnum: -pdgnum;
  }

  if (IsGluino())                            return 1000021;
  if (IsIno()) {
    if (kfc==kf::Neutralino1)                return 1000022; 
    if (kfc==kf::Neutralino2)                return 1000023; 
    if (kfc==kf::Neutralino3)                return 1000025; 
    if (kfc==kf::Neutralino4)                return 1000035; 
  }
  cerr<<"Error in Flavour::HepEvt() : No HepEvt number for "<<Flavour(kfc)<<endl;
  return 0;
}

void Flavour::FromHepEvt(int code) {
  anti = (code<0);
  code = abs(code);
  
  if ((code<23) || (code>100 && code<1000000) || (code>9000000)) {
    kfc = kf::code(code);
    return;
  }
  switch (code) {
  case 23:      kfc = kf::Z; return; 
  case 24:      kfc = kf::W; anti = 1-anti; return; 
  case 25: 
    if (Flavour(kf::h0).IsOn()) kfc = kf::h0; 
    else kfc = kf::h; 
    return;
  case 35:      kfc = kf::H0; return;
  case 36:      kfc = kf::A0; return;
  case 37:      kfc = kf::Hmin; anti = 1-anti; return; 
  case 39:      kfc = kf::graviton; return;
  case 89:      kfc = kf::gscalar; return;
  case 91:      kfc = kf::cluster; return;   // pythia cluster ....
  case 92:      kfc = kf::string; return;    // pythia string ....
  case 1000001: kfc = kf::sDownL; return;
  case 1000002: kfc = kf::sUpL; return;
  case 1000003: kfc = kf::sStrangeL; return;
  case 1000004: kfc = kf::sCharmL; return;
  case 1000005: kfc = kf::sBottom1; return;
  case 1000006: kfc = kf::sTop1; return;
  case 1000011: kfc = kf::sElectronL; return;
  case 1000012: kfc = kf::sNu1; return;
  case 1000013: kfc = kf::sMuL; return;
  case 1000014: kfc = kf::sNu2; return;
  case 1000015: kfc = kf::sTau1; return;
  case 1000016: kfc = kf::sNu3; return;
  case 2000001: kfc = kf::sDownR; return;
  case 2000002: kfc = kf::sUpR; return;
  case 2000003: kfc = kf::sStrangeR; return;
  case 2000004: kfc = kf::sCharmR; return;
  case 2000005: kfc = kf::sBottom2; return;
  case 2000006: kfc = kf::sTop2; return;
  case 2000011: kfc = kf::sElectronR; return;
  case 2000013: kfc = kf::sMuR; return;
  case 2000015: kfc = kf::sTau2; return;
  case 1000021: kfc = kf::Gluino; return;
  case 1000022: kfc = kf::Neutralino1; return;
  case 1000023: kfc = kf::Neutralino2; return;
  case 1000024: kfc = kf::Chargino1; return;
  case 1000025: kfc = kf::Neutralino3; return;
  case 1000035: kfc = kf::Neutralino4; return;
  case 1000037: kfc = kf::Chargino2; return;
  default: cerr<<"Error in Flavour::FromHepEvt() : No flavour for "<<code<<endl;
  }
  return;
}

std::string Flavour::TexName() 
{
  std::string name;

  if (anti) name = std::string("\\bar ");

  switch(kfc) {
  case kf::d: {name+=std::string("d");break;}
  case kf::u: {name+=std::string("u");break;}
  case kf::s: {name+=std::string("s");break;}
  case kf::c: {name+=std::string("c");break;}
  case kf::b: {name+=std::string("b");break;}
  case kf::t: {name+=std::string("t");break;}
  case kf::e: {name=std::string("e^\\m");break;}
  case kf::nue: {name += std::string("\\nu_e");break;}
  case kf::mu: {name=std::string("\\mu^\\m");break;}
  case kf::numu: {name+=std::string("\\nu_\\mu");break;}
  case kf::tau: {name= std::string("\\tau^\\m");break;}
  case kf::nutau: {name+= std::string("\\nu_\\tau");break;}
  case kf::gluon:  {name=std::string("g");break;}
  case kf::photon: {name= std::string("\\gamma");break;}
  case kf::W: {name=std::string("W^\\m");break;}
  case kf::Z: {name=std::string("Z^0");break;}
  case kf::h: {name=std::string("h");break;}
  case kf::h0: {name=std::string("h^0");break;}
  case kf::H0: {name=std::string("H^0");break;}
  case kf::A0: {name=std::string("A^0");break;}
  case kf::Hmin: {name=std::string("H^\\m");break;}
  case kf::Chargino1 :{name=std::string("\\chi^\\p_1");break;}
  case kf::Chargino2 :{name=std::string("\\chi^\\p_2");break;}
  case kf::Neutralino1 :{name=std::string("\\chi^0_1");break;}
  case kf::Neutralino2 :{name=std::string("\\chi^0_2");break;}
  case kf::Neutralino3 :{name=std::string("\\chi^0_3");break;}
  case kf::Neutralino4 :{name=std::string("\\chi^0_4");break;}
  case kf::Gluino :{name=std::string("\\tilde g");break;}
  case kf::sUpL :{name=std::string("\\tilde u_L");break;}
  case kf::sUpR :{name=std::string("\\tilde u_R");break;}
  case kf::sCharmL :{name=std::string("\\tilde c_L");break;}
  case kf::sCharmR :{name=std::string("\\tilde c_R");break;}
  case kf::sTop2 :{name=std::string("\\tilde t_2");break;}
  case kf::sTop1 :{name=std::string("\\tilde t_1");break;}
  case kf::sDownL :{name=std::string("\\tilde d_L");break;}
  case kf::sDownR :{name=std::string("\\tilde d_R");break;}
  case kf::sStrangeL :{name=std::string("\\tilde s_L");break;}
  case kf::sStrangeR :{name=std::string("\\tilde s_R");break;}
  case kf::sBottom2 :{name=std::string("\\tilde b_2");break;}
  case kf::sBottom1 :{name=std::string("\\tilde b_1");break;}
  case kf::sElectronL :{name=std::string("\\tilde e_L");break;}
  case kf::sElectronR :{name=std::string("\\tilde e_R");break;}
  case kf::sMuL :{name=std::string("\\tilde\\mu_L");break;}
  case kf::sMuR :{name=std::string("\\tilde\\mu_R");break;}
  case kf::sTau2 :{name=std::string("\\tilde\\tau_2");break;}
  case kf::sTau1 :{name=std::string("\\tilde\\tau_1");break;}
  case kf::sNu1 :{name=std::string("\\tilde\\nu_1");break;}
  case kf::sNu2 :{name=std::string("\\tilde\\nu_2");break;}
  case kf::sNu3 :{name=std::string("\\tilde\\nu_3");break;}
 
  default : break;
  }

  //Nomenklatur fuer Anti-Tilde-Teilchen nochmal ueberpruefen

  if (anti) {
  switch(kfc) {
  case kf::e: {name=std::string("e^\\p");break;}
  case kf::mu: {name=std::string("\\mu^\\p");break;}
  case kf::tau: {name= std::string("\\tau^\\p");break;}
  case kf::W: {name=std::string("W^\\p");break;}
  case kf::Hmin: {name=std::string("H^\\p");break;}
  case kf::Chargino1 :{name=std::string("\\chi^\\m_1");break;}
  case kf::Chargino2 :{name=std::string("\\chi^\\m_2");break;}
  case kf::Neutralino1 :{name=std::string("\\chi^0_1");break;}
  case kf::Neutralino2 :{name=std::string("\\chi^0_2");break;}
  case kf::Neutralino3 :{name=std::string("\\chi^0_3");break;}
  case kf::Neutralino4 :{name=std::string("\\chi^0_4");break;}
  case kf::Gluino :{name=std::string("\\tilde g^\\ti");break;}
  case kf::sUpL :{name=std::string("\\tilde u_L^\\ti");break;}
  case kf::sUpR :{name=std::string("\\tilde u_R^\\ti");break;}
  case kf::sCharmL :{name=std::string("\\tilde c_L^\\ti");break;}
  case kf::sCharmR :{name=std::string("\\tilde c_R^\\ti");break;}
  case kf::sTop2 :{name=std::string("\\tilde t_2^\\ti");break;}
  case kf::sTop1 :{name=std::string("\\tilde t_1^\\ti");break;}
  case kf::sDownL :{name=std::string("\\tilde d_L^\\ti");break;}
  case kf::sDownR :{name=std::string("\\tilde d_R^\\ti");break;}
  case kf::sStrangeL :{name=std::string("\\tilde s_L^\\ti");break;}
  case kf::sStrangeR :{name=std::string("\\tilde s_R^\\ti");break;}
  case kf::sBottom2 :{name=std::string("\\tilde b_2^\\ti");break;}
  case kf::sBottom1 :{name=std::string("\\tilde b_1^\\ti");break;}
  case kf::sElectronL :{name=std::string("\\tilde e_L^\\ti");break;}
  case kf::sElectronR :{name=std::string("\\tilde e_R^\\ti");break;}
  case kf::sMuL :{name=std::string("\\tilde\\mu_L^\\ti");break;}
  case kf::sMuR :{name=std::string("\\tilde\\mu_R^\\ti");break;}
  case kf::sTau2 :{name=std::string("\\tilde\\tau_2^\\ti");break;}
  case kf::sTau1 :{name=std::string("\\tilde\\tau_1^\\ti");break;}
  case kf::sNu1 :{name=std::string("\\tilde\\nu_1^\\ti");break;}
  case kf::sNu2 :{name=std::string("\\tilde\\nu_2^\\ti");break;}
  case kf::sNu3 :{name=std::string("\\tilde\\nu_^\\ti3");break;}  
  default : break;}
  }
  
  return name;
}

std::ostream& ATOOLS::operator<<(std::ostream& os, const Flavour& f)
{
  if(f.anti) {
    if (f==Flavour(kf::e).Bar())      return os<<"e+";
    if (f==Flavour(kf::mu).Bar())     return os<<"mu+";
    if (f==Flavour(kf::tau).Bar())    return os<<"tau+";
    if (f==Flavour(kf::W).Bar())      return os<<"W+";

    if (f==Flavour(kf::p_plus).Bar())  return os<<"P-";
    if (f==Flavour(kf::pi_plus).Bar()) return os<<"pi-";
    if (f==Flavour(kf::K_plus).Bar())  return os<<"K-";
    
    bool found = 1;
    string tmp = string(f.Name());
    int pos = tmp.find("+");
    if (pos>-1) {
      while (found) {
	pos = tmp.find("+");
	if (pos>-1) tmp.replace(pos,pos,string("-"));
	       else found = 0;
      }
      return os<<tmp.c_str();
    }

    pos = tmp.find("-");
    if (pos>-1) {
      while (found) {
	pos = tmp.find("-");
	if (pos>-1) tmp.replace(pos,pos,string("+"));
	       else found = 0;
      }
      return os<<tmp.c_str();
    }

    return os<<(string("anti-")+string(f.Name())).c_str();
  }
  return os<<f.Name();
}

// Definitions of all particles

void ATOOLS::ParticleInit(std::string path)
{
  Part_Info * pi = particles;

  int    kfc,charge,icharge,spin;
  int    kfcold=-1;
  bool   strong,Majorana,Take,stable,massive;
  double mass,width;
  char   name[20];
  char   buffer[150];

  int pc=0;

  *(pi++) = Part_Info( kf::start, -1,0,100000,0,0,0,0,0,1,0,
		       "Fatal error: particle type does not exist!",1 );
  ++pc;


  std::string filename = (path+std::string("/Particle.dat"));
  std::ifstream part(filename.c_str());
  if (!part) {
    std::cerr<<"Error in Particle_Init : File "<<filename<<" not found !"<<std::endl;
    return;
  }
  
  part.getline(buffer,150);


  for(;part;) {
    part>>kfc>>mass>>width>>charge>>icharge>>strong>>spin
	>>Majorana>>Take>>stable>>massive>>name;
    if (kfc!=kfcold) { // read last line only once!
      *(pi++)=Part_Info( kf::code(kfc), mass, width, charge, icharge, 
			 strong, spin, Majorana, Take, stable, massive, name, 1);
      //isrphoton
      if (kfc==22)
	*(pi++)=Part_Info( kf::code(26), mass, width, charge, icharge, 
			   strong, spin, Majorana, Take, stable, massive, "isrphoton",1);
      kfcold = kfc;
    }
  }
  part.close();

  filename=(path+std::string("/Hadron.dat"));
  std::ifstream part2(filename.c_str());
  if (!part2) {
    std::cerr<<"Error in Particle_Init : File "<<filename<<" not found !"<<std::endl;
  }
  else {
    part2.getline(buffer,150);
  
    for(;part2;) {
      part2>>kfc>>mass>>width>>charge>>icharge>>spin>>Take>>stable>>name;
      if (kfc!=kfcold) {  // read last line only once!	
	//	if (Take) {
	  ++pc;
	  *(pi++)=Part_Info( kf::code(kfc), mass, width, charge, icharge, 
			     spin, Take, stable, name);
	  kfcold=kfc;
	  //	}
      }
    }
    part2.close();
  }
  // kfcode,mass,width,charge,icharge,strong,spin,majorana,take,stable,massive,name,1
  *(pi++) = Part_Info( kf::pol,      0.,0., 0, 0,0, 0,0,0,1,0,"polarisation",1);
  *(pi++) = Part_Info( kf::lepton,   0.,0.,-3,-1,0, 1,0,1,1,0,"lepton",1);
  *(pi++) = Part_Info( kf::neutrino, 0.,0., 0, 1,0, 1,0,1,1,0,"neutrino",1);
  *(pi++) = Part_Info( kf::fermion,  0.,0., 0, 0,0, 1,0,1,1,0,"fermion",1);
  *(pi++) = Part_Info( kf::jet,      0.,0., 0, 0,1, 2,0,1,1,0,"jet",1);
  *(pi++) = Part_Info( kf::quark,    0.,0., 0, 0,1, 1,0,1,1,0,"Quark",1);

  *(pi++) = Part_Info( kf::cluster,    0.,0., 0, 0,0, 0,0,1,1,0,"cluster",1);
  *(pi++) = Part_Info( kf::string,    0.,0., 0, 0,0, 0,0,1,1,0,"string",1);

  *(pi++)=Part_Info( kf::none, -1,0,0,0,0,0,0,0,1,0, "no particle",1);

  kf_table.Init(); 


  Flavour addit;
  Flavour test = Flavour(kf::jet);
  int idx      = kf_table.ToInt(test.Kfcode());
  particles[idx].SetGroup();
  for (int i=1;i<7;i++) {
    addit = Flavour(kf::code(i));
    if (addit.Mass()==0.) {
      particles[idx].Add(addit);
      particles[idx].Add(addit.Bar());
    }
  }
  particles[idx].Add(Flavour(kf::gluon));

  test = Flavour(kf::quark);
  idx      = kf_table.ToInt(test.Kfcode());
  particles[idx].SetGroup();
  for (int i=1;i<7;i++) {
    addit = Flavour(kf::code(i));
    if (addit.Mass()==0.) {
      particles[idx].Add(addit);
      particles[idx].Add(addit.Bar());
    }
  }

  test = Flavour(kf::lepton);
  idx      = kf_table.ToInt(test.Kfcode());
  particles[idx].SetGroup();
  for (int i=11;i<17;i+=2) {
    addit = Flavour(kf::code(i));
    if (addit.Mass()==0.) {
      particles[idx].Add(addit);
      particles[idx].Add(addit.Bar());
    }
  }
  test = Flavour(kf::neutrino);
  idx      = kf_table.ToInt(test.Kfcode());
  particles[idx].SetGroup();
  for (int i=12;i<18;i+=2) {
    addit = Flavour(kf::code(i));
    if (addit.Mass()==0.) {
      particles[idx].Add(addit);
      particles[idx].Add(addit.Bar());
    }
  }
}
  
// Unique Identifier of Particle.dat
// aaaaaa-bbbbbb-ccccccc 
// a   SM   - particles
// b    MSSM - particles
// c    everything else
// each
// (n)nn  number of particles in this section 
//      xxxx 
// double 0.mmmmm00000 * 10^ee => int ((ee+mmmmm)&0xFFFF)rot(kf::code * number of flags + flag)

// Switch off masses for SM fermions

void SetMassless()
{
  particles[ kf_table.ToInt(kf::d) ].msv=0;
  particles[ kf_table.ToInt(kf::u) ].msv=0;
  particles[ kf_table.ToInt(kf::s) ].msv=0;
  particles[ kf_table.ToInt(kf::c) ].msv=0;
  particles[ kf_table.ToInt(kf::b) ].msv=0;
  //particles[ kf_table.ToInt(kf::t) ].msv=0;

  particles[ kf_table.ToInt(kf::e) ].msv=0;
  particles[ kf_table.ToInt(kf::mu) ].msv=0;
  particles[ kf_table.ToInt(kf::tau) ].msv=0;
  particles[ kf_table.ToInt(kf::nue) ].msv=0;
  particles[ kf_table.ToInt(kf::numu) ].msv=0;
  particles[ kf_table.ToInt(kf::nutau) ].msv=0;
}

int Flavour::GetIntID(double value) {
  // double 0.mmmmm00000 * 10^ee
  // negative exponent possible
  // negative mantisse is not included!!!!
  if (value==0) return 0;

  int length = 4;
  int exp_10=int(log(value)/log(10.));
  // quick
  double nref=pow(10.,exp_10)*pow(2.,-length*4+4);
  if (exp_10>=0) exp_10=exp_10&0x7F;
  else exp_10=0x80+((-exp_10)&0x7F);

  int result = int (value/nref+0.5)*0x100+exp_10;
  return result;
}


int Flavour::PropertiesID() {
// IntID(mass) + Rot(IntID(width),8,16)  xor flags
  int pid=GetIntID(PSMass());
  int wid=GetIntID(Width());
  pid^=Rot(wid,12,24);

  pid=pid^FlagID();
  return pid;
}

int Flavour::FlagID() {
  //ch+4 -3, -1, 0, 2       
  //    1   3   4  6
  //
  //y+1 0 1 2
  //s   0 1 2
  //
  //mj   0 1
  //on  0 1
  //t   0 1
  //M   0 1 
  //
  // su3 0 1
  int ix=Index();

  int flag=0;
  flag+= IntCharge() +4;
  flag*=16;
  flag+= (particles[ix].isow +1)*4 + particles[ix].sp;
  flag*=16;
  flag+=  particles[ix].Maj*8+ particles[ix].on*4+ particles[ix].stbl*2 + particles[ix].msv;
  flag*=16;
  flag+=  particles[ix].str;
  return flag;
}


int Flavour::Rot( int value , int rotby, int length) {
  int cut=1;
  rotby=rotby%length;
  cut=cut<<length; 
  int mask=cut-1;
  value=value&mask;
  for (int i=0; i<rotby; ++i) {
    value=value<<1;
    if (value&cut) 
      value=(value|1)&mask;
  }
  return value;
}

int Flavour::ID_SM() {
  Fl_Iter fli;
  int count=0;
  int sm_id=0;

  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.IsOn() && !flav.IsHadron() && !flav.IsSusy()) {
      count++;
      int local_id=flav.PropertiesID();
      local_id=Rot(local_id,flav.Kfcode(),24);
      sm_id=sm_id^local_id;
    } 
  }
  return sm_id+(count*0x1000000);
}

int Flavour::ID_MSSM() {
  Fl_Iter fli;
  int count=0;
  int mssm_id=0;

  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.IsOn() && !flav.IsHadron() && flav.IsSusy()) {
      count++;
      int local_id=flav.PropertiesID();
      local_id=Rot(local_id,flav.Kfcode(),24);
      mssm_id=mssm_id^local_id;
    } 
  }
  return mssm_id+(count*0x1000000);

}

int Flavour::ID_Had() {
  Fl_Iter fli;
  int count=0;
  int had_id=0;

  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.IsOn() && flav.IsHadron() && !flav.IsSusy()) {
      count++;
      int local_id=flav.PropertiesID();
      local_id=Rot(local_id,flav.Kfcode(),24);
      had_id=had_id^local_id;
    } 
  }
  return had_id+(count*0x1000000);
}

int Flavour::WriteOut() {
  std::string sm_name, mssm_name, had_name;
  MyStrStream str;

  str<<ID_SM();
  str>>sm_name;
  str<<ID_MSSM();
  str>>mssm_name;
  str<<ID_Had();
  str>>had_name;

  std::ofstream ofile;

  ofile.open((std::string("save/sm/")+sm_name+std::string(".dat")).c_str());

  Fl_Iter fli;
  ofile<<"kf      Mass   Width    3*e     Y     SU(3)  2*Spin majorana   ON    stabil   massive   Name"<<std::endl;
  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.IsOn() && !flav.IsHadron() && !flav.IsSusy() && !flav.Kfcode()!=26) {
      ofile<<flav.Kfcode()<<"\t"
	  <<flav.PSMass()<<"\t"
	  <<flav.Width()<<"\t"
	  <<flav.IntCharge()<<"\t"
	  <<flav.IsoWeak()*2.<<"\t"
	  <<flav.Strong()<<"\t"
	  <<flav.IntSpin()<<"\t"
	  <<flav.Majorana()<<"\t"
	  <<flav.IsOn()<<"\t"
	  <<flav.IsStable()<<"\t"
	  <<flav.IsMassive()<<"\t"
	  <<flav.Name()<<std::endl;
    } 
  }
  ofile.close();

  return 1;
  


}

