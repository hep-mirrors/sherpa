#ifndef AHADIC_Tools_Proto_Particle_H
#define AHADIC_Tools_Proto_Particle_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include <list>
#include <iostream>

namespace AHADIC {
  class Proto_Particle;
  class Proto_Particle {
    static std::list<Proto_Particle * > s_actives;
  public:
    ATOOLS::Flavour  m_flav;
    ATOOLS::Vec4D    m_mom;
    char             m_info;
    double           m_mass, m_kt2max;
    Proto_Particle * p_partner;
  public:
    Proto_Particle(const Proto_Particle &);  
    Proto_Particle(ATOOLS::Flavour flav,ATOOLS::Vec4D mom,char info);  
    ~Proto_Particle();

    bool CheckConsistency(std::ostream & s,std::string method);

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Proto_Particles -----"<<std::endl;
      for (std::list<Proto_Particle * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->m_flav<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };
  std::ostream & operator<<(std::ostream & s, const Proto_Particle &);



  class Proto_Particle_List;
  class Proto_Particle_List : public std::list<Proto_Particle * > {
  private:
    static std::list<Proto_Particle_List * > s_actives;
  public:
    Proto_Particle_List() {
      s_actives.push_back(this);
    }

    Proto_Particle_List(const Proto_Particle_List & pl) {
      s_actives.push_back(this);
      for (std::list<Proto_Particle * >::const_iterator piter=pl.begin();
	   piter!=pl.end();piter++)
	push_back((*piter));
    }

    ~Proto_Particle_List() {
      s_actives.remove(this);
    }

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active Proto_Particle_Lists -----"<<std::endl;
      for (std::list<Proto_Particle_List * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->size()<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };

  typedef Proto_Particle_List::iterator       PPL_Iterator; 
  typedef Proto_Particle_List::const_iterator PPL_Const_Iterator; 
  std::ostream & operator<<(std::ostream & s, const Proto_Particle_List & pl);


  class ListOfPPLs;
  class ListOfPPLs : public std::list<Proto_Particle_List * > {
  private:
    static std::list<ListOfPPLs * > s_actives;
  public:
    ListOfPPLs() {
      s_actives.push_back(this);
    }

    ~ListOfPPLs() {
      s_actives.remove(this);
    }

    static void PrintActives(std::ostream & s) {
      s<<"   ----- Active ListOfPPLss -----"<<std::endl;
      for (std::list<ListOfPPLs * >::iterator sit=s_actives.begin();
	   sit!=s_actives.end(); sit++) {
	s<<"   ("<<(*sit)->size()<<"/"<<(*sit)<<")"<<std::endl;
      }
      s<<"   ----- Total: "<<s_actives.size()<<" -----"<<std::endl;
    };
    
    static long int RemainingActives() {
      return s_actives.size();
    };
    
    static void DeleteActives() {
      if (s_actives.empty()) return;
      while (!s_actives.empty()) delete s_actives.front();
      s_actives.clear();
    };
  };

  typedef ListOfPPLs::iterator                LPPL_Iterator;

  struct control {
  public:
    static long int s_AHAparticles;
    static long int s_AHAprotoparticles;
    static long int s_AHAblobs;
  };

  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;
}

#endif
