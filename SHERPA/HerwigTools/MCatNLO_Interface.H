#ifndef SHERPA_HerwigTools_MCatNLO_Interface_H
#define SHERPA_HerwigTools_MCatNLO_Interface_H

#include "SHERPA/HerwigTools/Herwig_Interface.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "PDF/Main/ISR_Handler.H"
#include <map>

namespace SHERPA {
  struct ptp {
    enum code {
      Unspecified  =  0,
      Single_Boson =  1,
      Lepton_Pair  =  2,
      Higgs        = 11,
      Heavy_Quark  = 21,
      Boson_Pair   = 31
    };
  };
  std::ostream& operator<<(std::ostream&,const ptp::code);


  class MCatNLO_Interface {
  private:
    std::string        m_path, m_file;
    ptp::code          m_proc;
    Herwig_Interface * p_herwig;
    bool ReadInTheParameters();
    void Shorten(std::string&);
  public:
    MCatNLO_Interface(std::string _m_path,std::string _m_file,bool sherpa);
    ~MCatNLO_Interface();
    bool Initialize();

    bool OneEvent(ATOOLS::Blob_List *const bloblist,double &weight);
    void Terminate();
  };

  //////////////////
  // Doxygen part //
  //////////////////

  /*!
    \file MCatNLO_Interface.H
    \brief Declares the class WRAPPER::MCatNLO_Interface 
  */

  /*!
    \class MCatNLO_Interface 
    \brief Provides steering of \b MC@NLO
  */

  /*!
    \fn MCatNLO_Interface::Init(std::string,std::string)
    \brief Initializes \b MC@NLO
  */

  /*!
    \fn MCatNLO_Interface::GenerateEvent(bool)
    \brief Generates one event with \b MC@NLO
  */

  /*!
    \fn MCatNLO_Interface::Terminate()
    \brief Terminates \b MC@NLO
  */

}

#endif
