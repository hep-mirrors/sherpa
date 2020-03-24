#include "ATOOLS/Phys/Blob_Data.H"

using namespace ATOOLS;

template <class Type>
Blob_Data<Type>::~Blob_Data() 
{
}

template <class Type>
Blob_Data_Base* Blob_Data<Type>::ClonePtr()
{
  return new Blob_Data(m_data);
}
