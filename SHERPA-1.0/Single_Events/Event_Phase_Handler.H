#ifndef Event_Phase_Handler_h
#define Event_Phase_Handler_h

#include <set> 
#include <iterator> 
#include <iostream.h> 

#include "Parton_List.H"
#include "Parton.H"
#include "Blob_List.H"
#include "Blob.H"

namespace SHERPA {
  class Event_Phase_Handler {
  protected :
    std::string m_type,m_name;

    APHYTOOLS::Blob_List    m_myblobs;
    APHYTOOLS::Parton_List  m_mypartons;
    virtual void EmptyMyLists();
  public :
    Event_Phase_Handler();
    virtual ~Event_Phase_Handler();
    virtual bool Treat(APHYTOOLS::Blob_List *) = 0;
    virtual void CleanUp()                     = 0;

    virtual int                      NumberOfBlobs();
    virtual APHYTOOLS::Blob        * GetBlob(int);         
    virtual APHYTOOLS::Blob_List   * GetBlobs();
    virtual void                     AddBlob(APHYTOOLS::Blob *);
    virtual int                      NumberOfPartons();
    virtual APHYTOOLS::Parton      * GetParton(int);         
    virtual APHYTOOLS::Parton_List * GetPartons();
    virtual void                     AddParton(APHYTOOLS::Parton *);
    virtual std::string              Type();
    virtual void                     SetType(std::string);
    virtual std::string              Name();
    virtual void                     SetName(std::string);
  };
}

#endif

