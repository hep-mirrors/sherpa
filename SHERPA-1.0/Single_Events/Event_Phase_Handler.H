#ifndef Event_Phase_Handler_h
#define Event_Phase_Handler_h

#include <set> 
#include <iterator> 
#include <iostream> 

#include "Parton_List.H"
#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Event_Phase_Handler {
  protected :
    std::string m_type,m_name;
  public :
    Event_Phase_Handler();
    virtual ~Event_Phase_Handler();
    virtual bool Treat(APHYTOOLS::Blob_List *) = 0;
    virtual void CleanUp()                     = 0;
    virtual std::string              Type();
    virtual void                     SetType(std::string);
    virtual std::string              Name();
    virtual void                     SetName(std::string);
  };
}

#endif

