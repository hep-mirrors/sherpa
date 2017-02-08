#include "ATOOLS/Org/Smart_Pointer.H"

#include "ATOOLS/Org/Shell_Tools.H"

template <class Class_Type> const ATOOLS::Smart_Pointer<Class_Type> *
ATOOLS::Smart_Pointer<Class_Type>::FindOwner() const
{
  if (p_owner) return p_owner->FindOwner();
  return this;
}

template <class Class_Type>
void ATOOLS::Smart_Pointer<Class_Type>::PrintForward
(std::ostream &str,const bool all) const
{
  if (all) {
    str<<"("<<this<<")["<<Demangle(typeid(p_this).name())
       <<"]: p_this = "<<p_this<<" {\n";
    FindOwner()->PrintForward(str,false); str<<"}";
    return;
  }
  str<<"   ("<<this<<"): { p_owner = "<<p_owner
     <<", p_copy = "<<p_copy<<" }\n";
  if (p_copy!=NULL) p_copy->PrintForward(str);
}
