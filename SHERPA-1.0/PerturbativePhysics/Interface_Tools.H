#ifndef InterfaceTools_h
#define InterfaceTools_h

#include <set> 
#include <iterator> 
#include <iostream> 


#include "XS_Base.H"

#include "Tree.H"
#include "Knot.H"
#include "Parton_List.H"
#include "Parton.H"
#include "Flow.H"
#include "Blob.H"


namespace SHERPA {
  class Interface_Tools {
    APACIC::Tree ** p_initrees; 
    APACIC::Tree  * p_fintree; 

    bool IsColourConnected(APHYTOOLS::Parton *,APHYTOOLS::Parton *) ;
  public :
    Interface_Tools(APACIC::Tree **,APACIC::Tree *); 
    ~Interface_Tools(); 

    void InitializeIncoming(APHYTOOLS::Blob *,double,double,double);
    void InitializeOutGoing(APHYTOOLS::Blob *,double,double,double);

    double ColourAngle(APHYTOOLS::Parton *,APHYTOOLS::Blob *);
  };
}
#endif




