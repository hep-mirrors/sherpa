/*  Flavour.C: Definitions for Flavour class
 *  A Class to define particle properties.
 */

#include <iostream>

#include "Flavour.H"
//#include "Run_Parameter.H"
#include "MathTools.H"
#include "Message.H"
#include "MyStrStream.H"

namespace APHYTOOLS {
  part_info Particles[MAX_PARTICLES];
  kf_to_int KF_table;
}

using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;

int APHYTOOLS::kf_to_int::is_initialised=0;

void kf_to_int::init()
{
  for(anz=0;anz<MAX_PARTICLES;++anz) {
    kf_tab[anz] = Particles[anz].kfc;
    if(kf_tab[anz]==kf::none) break;
  }

  if(anz==MAX_PARTICLES) {
    std::cerr<<"ERROR in kf_to_int::kftab(): Too many particle types !"<<std::endl;
    exit(1);
  }
  is_initialised = 1;
}

kf::code kf_to_int::from_int(int code)
{
  return kf_tab[code];
}

int kf_to_int::to_int(kf::code kfc)
{
  for(int i=0;i<anz+1;i++)
    if(kf_tab[i]==kfc) return i;

  std::cerr<<"ERROR in kf_to_int::to_int(): Particle type unknown ! "
	   <<kfc<<" in "<<anz<<std::endl;
  return -1;
}

kf::code kf_to_int::from_string(std::string st)
{
  for(int i=0;i<anz;i++) {
    if(std::string(Particles[i].n)==st) return kf_tab[i];
  }

  std::cerr<<"ERROR in kf_to_int::from_string("<<st<<"): Particle type unknown !"<<std::endl;
  return kf::none;
}




part_info::part_info(kf::code kfc_, 
		     double mass, double width, int icharge, int isoweak, 
		     bool strong, int spin, bool Majorana, 
		     bool Take, bool stable,bool massive,
		     char* name,int _masssign) : 
  kfc(kfc_), m(mass), yuk(mass), w(width), iq(icharge), isow(isoweak), 
  str(strong), sp(spin), Maj(Majorana), on(Take), stbl(stable), msv(massive),hadron(0),
  n(name), masssign(_masssign), group(0) { 
  n  = new char[strlen(name)+1];
  strcpy(n,name);
}

part_info::part_info(kf::code kfc_, 
		     double mass, double width, int icharge, int isoweak, 
		     int spin, bool Take, bool stable,
		     char* name) : 
  kfc(kfc_), m(mass), yuk(0.),  w(width), iq(icharge), isow(isoweak), 
  str(0), sp(spin), Maj(0), on(Take), stbl(stable), msv(1), hadron(1),
  n(name), masssign(1), group(0) { 
  n  = new char[strlen(name)+1];
  strcpy(n,name);
}

void part_info::Add(const Flavour _fl) {
  if (_fl.Size()==1) {
    flavs.push_back(new Flavour(_fl));
    return;
  }
  for (int i=0;i<_fl.Size();i++) {
    flavs.push_back(new Flavour(_fl[i]));
  }
}

const int part_info::Size() {
  if (!group) return 1;
  return flavs.size();
}

const bool part_info::Includes(Flavour _fl) {
  bool isin;
  for (int j=0;j<_fl.Size();j++) {
    isin = 0;
    for (int i=0;i<flavs.size();i++) {
      if (flavs[i]->kfcode() == _fl[j].kfcode()) { isin = 1; break; }
    }
    if (isin==0) return 0;
  }
  return 1;
}

Flavour part_info::operator[](const int i) {
  if (!group) return Flavour(kfc);
  if ((i<0) || (i>flavs.size())) {
    std::cerr<<"Error in part_info::operator[]("<<i<<") . Delimiter out of bounds."<<std::endl;
    return Flavour(kfc);
  }
  return (*flavs[i]);
}

void part_info::SetGroup() { group = 1; }




Flavour Flavour::operator[](int idx) const  { 
  if (!(Particles[ index() ].group)) return (*this);
  if (anti) return Particles[ index() ][idx].bar();
  return Particles[ index() ][idx];
} 



Flavour Flavour::bar() { 
  Flavour flbar = *this;
  
  if (flbar==Flavour(kf::photon)) return flbar;
  if (flbar==Flavour(kf::Z))      return flbar;
  if (flbar==Flavour(kf::gluon))  return flbar;
  if (flbar==Flavour(kf::h))      return flbar;
  if (flbar==Flavour(kf::h0))     return flbar;
  if (flbar==Flavour(kf::H0))     return flbar;
  if (flbar==Flavour(kf::A0))     return flbar;
  if (flbar==Flavour(kf::none))   return flbar;
  if (Majorana())                 return flbar;

  flbar.anti = (anti)?0:1;
  return flbar;
}

Flavour::Flavour(int kf) { 
  anti = (kf<0)? 1:0;
  kfc  = KF_table.from_int(abs(kf));
}


int Flavour::ctq() {
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

void Flavour::from_ctq(int code) {
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

int Flavour::hepevt() {
  int code;
  switch(kfc) {
  case kf::d      : code=1;break;
  case kf::u      : code=2;break;
  case kf::s      : code=3;break;
  case kf::c      : code=4;break;
  case kf::b      : code=5;break;
  case kf::t      : code=6;break;
  case kf::gluon  : code=21;break;
  case kf::photon : code=22;break;
  default         : code=0;
  }
  if(anti) code= -code;
  return code;
}

void Flavour::from_hepevt(int code) {
  anti=0;
  if (code<0) {code=-code;anti=1;}
  switch(abs(code)) {
  case 1   : kfc=kf::d;break;
  case 2   : kfc=kf::u;break;
  case 3   : kfc=kf::s;break;
  case 4   : kfc=kf::c;break;
  case 5   : kfc=kf::b;break;
  case 7   : kfc=kf::photon;break;
  case 21  : kfc=kf::gluon;break;
  default  : kfc=kf::none;
  }
}

std::string Flavour::texname() 
{
  std::string name;

  if (anti) name+=std::string("\\bar ");

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
  case kf::Chargino1 :{name=std::string("\\chi^\\m_1");break;}
  case kf::Chargino2 :{name=std::string("\\chi^\\m_2");break;}
  case kf::Neutralino1 :{name=std::string("\\chi^0_1");break;}
  case kf::Neutralino2 :{name=std::string("\\chi^0_2");break;}
  case kf::Neutralino3 :{name=std::string("\\chi^0_3");break;}
  case kf::Neutralino4 :{name=std::string("\\chi^0_4");break;}
  case kf::Gluino :{name=std::string("\\tilde g");break;}
  case kf::sUpL :{name=std::string("\\tilde u_L");break;}
  case kf::sUpR :{name=std::string("\\tilde u_R");break;}
  case kf::sCharmL :{name=std::string("\\tilde c_L");break;}
  case kf::sCharmR :{name=std::string("\\tilde c_R");break;}
  case kf::sTopL :{name=std::string("\\tilde t_L");break;}
  case kf::sTopR :{name=std::string("\\tilde t_R");break;}
  case kf::sDownL :{name=std::string("\\tilde d_L");break;}
  case kf::sDownR :{name=std::string("\\tilde d_R");break;}
  case kf::sStrangeL :{name=std::string("\\tilde s_L");break;}
  case kf::sStrangeR :{name=std::string("\\tilde s_R");break;}
  case kf::sBottomL :{name=std::string("\\tilde b_L");break;}
  case kf::sBottomR :{name=std::string("\\tilde b_R");break;}
  case kf::sElectronL :{name=std::string("\\tilde e_L");break;}
  case kf::sElectronR :{name=std::string("\\tilde e_R");break;}
  case kf::sMuL :{name=std::string("\\tilde\\mu_L");break;}
  case kf::sMuR :{name=std::string("\\tilde\\mu_R");break;}
  case kf::sTauL :{name=std::string("\\tilde\\tau_L");break;}
  case kf::sTauR :{name=std::string("\\tilde\\tau_R");break;}
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
  case kf::Chargino1 :{name=std::string("\\chi^\\p_1");break;}
  case kf::Chargino2 :{name=std::string("\\chi^\\p_2");break;}
  case kf::Neutralino1 :{name=std::string("\\chi^0_1");break;}
  case kf::Neutralino2 :{name=std::string("\\chi^0_2");break;}
  case kf::Neutralino3 :{name=std::string("\\chi^0_3");break;}
  case kf::Neutralino4 :{name=std::string("\\chi^0_4");break;}
  case kf::Gluino :{name=std::string("\\tilde g^\\ti");break;}
  case kf::sUpL :{name=std::string("\\tilde u_L^\\ti");break;}
  case kf::sUpR :{name=std::string("\\tilde u_R^\\ti");break;}
  case kf::sCharmL :{name=std::string("\\tilde c_L^\\ti");break;}
  case kf::sCharmR :{name=std::string("\\tilde c_R^\\ti");break;}
  case kf::sTopL :{name=std::string("\\tilde t_L^\\ti");break;}
  case kf::sTopR :{name=std::string("\\tilde t_R^\\ti");break;}
  case kf::sDownL :{name=std::string("\\tilde d_L^\\ti");break;}
  case kf::sDownR :{name=std::string("\\tilde d_R^\\ti");break;}
  case kf::sStrangeL :{name=std::string("\\tilde s_L^\\ti");break;}
  case kf::sStrangeR :{name=std::string("\\tilde s_R^\\ti");break;}
  case kf::sBottomL :{name=std::string("\\tilde b_L^\\ti");break;}
  case kf::sBottomR :{name=std::string("\\tilde b_R^\\ti");break;}
  case kf::sElectronL :{name=std::string("\\tilde e_L^\\ti");break;}
  case kf::sElectronR :{name=std::string("\\tilde e_R^\\ti");break;}
  case kf::sMuL :{name=std::string("\\tilde\\mu_L^\\ti");break;}
  case kf::sMuR :{name=std::string("\\tilde\\mu_R^\\ti");break;}
  case kf::sTauL :{name=std::string("\\tilde\\tau_L^\\ti");break;}
  case kf::sTauR :{name=std::string("\\tilde\\tau_R^\\ti");break;}
  case kf::sNu1 :{name=std::string("\\tilde\\nu_1^\\ti");break;}
  case kf::sNu2 :{name=std::string("\\tilde\\nu_2^\\ti");break;}
  case kf::sNu3 :{name=std::string("\\tilde\\nu_^\\ti3");break;}  
  default : break;}
  }
  
  return name;
}

std::ostream& APHYTOOLS::operator<<(std::ostream& os, const Flavour& f)
{
  if(f.anti) {
    if (f==Flavour(kf::e).bar())      return os<<"e+";
    if (f==Flavour(kf::mu).bar())     return os<<"mu+";
    if (f==Flavour(kf::tau).bar())    return os<<"tau+";
    if (f==Flavour(kf::W).bar())      return os<<"W+";

    if (f==Flavour(kf::p_plus).bar()) return os<<"P-";
    if (f==Flavour(kf::pi_plus).bar()) return os<<"pi-";
    if (f==Flavour(kf::K_plus).bar()) return os<<"K-";
    os<<"anti-";
  }
  return os<<f.name();
}

// Definitions of all particles

void APHYTOOLS::particle_init(std::string path)
{
  part_info * pi = Particles;

  int    kfc,charge,icharge,spin;
  int    kfcold=-1;
  bool   strong,Majorana,Take,stable,massive;
  double mass,width;
  char   name[20];
  char   buffer[150];

  int pc=0;

  *(pi++) = part_info( kf::start, -1,0,100000,0,0,0,0,0,1,0,
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
      *(pi++)=part_info( kf::code(kfc), mass, width, charge, icharge, 
			 strong, spin, Majorana, Take, stable, massive, name, 1);
      //isrphoton
      if (kfc==22)
	*(pi++)=part_info( kf::code(26), mass, width, charge, icharge, 
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
 
	  *(pi++)=part_info( kf::code(kfc), mass, width, charge, icharge, 
			   spin, Take, stable, name);
	  kfcold=kfc;
	  //	}
      }
    }
    part2.close();
  }
  // kfcode,mass,width,charge,icharge,strong,spin,majorana,take,stable,massive,name,1
  *(pi++) = part_info( kf::pol,      0.,0., 0, 0,0, 0,0,0,1,0,"polarisation",1);
  *(pi++) = part_info( kf::lepton,   0.,0.,-3,-1,0, 2,0,1,1,0,"lepton",1);
  *(pi++) = part_info( kf::neutrino, 0.,0., 0, 1,0, 2,0,1,1,0,"neutrino",1);
  *(pi++) = part_info( kf::fermion,  0.,0., 0, 0,0, 2,0,1,1,0,"fermion",1);
  *(pi++) = part_info( kf::jet,      0.,0., 0, 0,1, 0,0,1,1,0,"jet",1);
  *(pi++) = part_info( kf::quark,    0.,0., 0, 0,1, 2,0,1,1,0,"quark",1);

  *(pi++)=part_info( kf::none, -1,0,0,0,0,0,0,0,1,0, "no particle",1);

  KF_table.init(); 


  Flavour test = Flavour(kf::jet);
  int idx      = KF_table.to_int(test.kfcode());
  Particles[idx].SetGroup();
  Particles[idx].Add(Flavour(kf::d));
  Particles[idx].Add(Flavour(kf::u));
  Particles[idx].Add(Flavour(kf::s));
  Particles[idx].Add(Flavour(kf::c));
  Particles[idx].Add(Flavour(kf::b));
  Particles[idx].Add(Flavour(kf::d).bar());
  Particles[idx].Add(Flavour(kf::u).bar());
  Particles[idx].Add(Flavour(kf::s).bar());
  Particles[idx].Add(Flavour(kf::c).bar());
  Particles[idx].Add(Flavour(kf::b).bar());
  Particles[idx].Add(Flavour(kf::gluon));

  test = Flavour(kf::quark);
  idx      = KF_table.to_int(test.kfcode());
  Particles[idx].SetGroup();
  Particles[idx].Add(Flavour(kf::d));
  Particles[idx].Add(Flavour(kf::u));
  Particles[idx].Add(Flavour(kf::s));
  Particles[idx].Add(Flavour(kf::c));
  Particles[idx].Add(Flavour(kf::b));
  Particles[idx].Add(Flavour(kf::d).bar());
  Particles[idx].Add(Flavour(kf::u).bar());
  Particles[idx].Add(Flavour(kf::s).bar());
  Particles[idx].Add(Flavour(kf::c).bar());
  Particles[idx].Add(Flavour(kf::b).bar());

  test = Flavour(kf::lepton);
  idx      = KF_table.to_int(test.kfcode());
  Particles[idx].SetGroup();
  Particles[idx].Add(Flavour(kf::e));
  Particles[idx].Add(Flavour(kf::mu));
  Particles[idx].Add(Flavour(kf::tau));
  Particles[idx].Add(Flavour(kf::e).bar());
  Particles[idx].Add(Flavour(kf::mu).bar());
  Particles[idx].Add(Flavour(kf::tau).bar());

  std::cout<<"End of Particle_Init. Initialised the KF_table. "<<std::endl;
}

// Unique Identifier of Particle.dat
// aaaaaa-bbbbbb-ccccccc 
// a   SM   - partons
// b    MSSM - partons
// c    everyting els
// each
// (n)nn  numer of particles in this section 
//      xxxx 
// double 0.mmmmm00000 * 10^ee => int ((ee+mmmmm)&0xFFFF)rot(kf::code * number of flags + flag)

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
  //  cout<<" doubleID="<<result<<std::endl;
  return result;
}


int Flavour::PropertiesID() {
// IntID(mass) + Rot(IntID(width),8,16)  xor flags
  int pid=GetIntID(PSmass());
  int wid=GetIntID(width());
  pid^=Rot(wid,12,24);

  pid=pid^FlagID();
  //  cout<<" pid="<<pid;
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
  int ix=index();

  int flag=0;
  flag+= icharge() +4;
  flag*=16;

  flag+= (Particles[ix].isow +1)*4 + Particles[ix].sp;
  flag*=16;

  flag+=  Particles[ix].Maj*8+ Particles[ix].on*4+ Particles[ix].stbl*2 + Particles[ix].msv;
  flag*=16;

  flag+=  Particles[ix].str;

  //  cout<<" flags="<<flag<<std::endl;
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
  fl_iter fli;
  int count=0;
  int sm_id=0;

  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.ison() && !flav.ishadron() && !flav.issusy()) {
      count++;
      //      cout<<flav<<":"<<std::endl;
      int local_id=flav.PropertiesID();
      // Rot ( all single IDs rot by kfcode )
      local_id=Rot(local_id,flav.kfcode(),24);
      // combine all results by XOR
      sm_id=sm_id^local_id;
      //      cout<<" sm_id="<<sm_id<<std::endl;
    } 
  }
  return sm_id+(count*0x1000000);
}

int Flavour::ID_MSSM() {
  fl_iter fli;
  int count=0;
  int mssm_id=0;

  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.ison() && !flav.ishadron() && flav.issusy()) {
      count++;
      //      cout<<flav<<":"<<std::endl;
      int local_id=flav.PropertiesID();
      local_id=Rot(local_id,flav.kfcode(),24);
      mssm_id=mssm_id^local_id;
      //      cout<<" mssm_id="<<mssm_id<<std::endl;
    } 
  }
  return mssm_id+(count*0x1000000);

}

int Flavour::ID_Had() {
  fl_iter fli;
  int count=0;
  int had_id=0;

  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.ison() && flav.ishadron() && !flav.issusy()) {
      count++;
      //      cout<<flav<<":"<<std::endl;
      int local_id=flav.PropertiesID();
      local_id=Rot(local_id,flav.kfcode(),24);
      had_id=had_id^local_id;
      //      cout<<" had_id="<<had_id<<std::endl;
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

  fl_iter fli;
  ofile<<"kf      Mass   Width    3*e     Y     SU(3)  2*Spin majorana   ON    stabil   massive   Name"<<std::endl;
  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.ison() && !flav.ishadron() && !flav.issusy() && !flav.kfcode()!=26) {
      ofile<<flav.kfcode()<<"\t"
	  <<flav.PSmass()<<"\t"
	  <<flav.width()<<"\t"
	  <<flav.icharge()<<"\t"
	  <<flav.isoweak()*2.<<"\t"
	  <<flav.strong()<<"\t"
	  <<flav.ispin()<<"\t"
	  <<flav.Majorana()<<"\t"
	  <<flav.ison()<<"\t"
	  <<flav.isstable()<<"\t"
	  <<flav.ismassive()<<"\t"
	  <<flav.name()<<std::endl;
    } 
  }
  ofile.close();

  return 1;
  


}

