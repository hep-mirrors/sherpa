#include "ATOOLS/Phys/Spin_Density_Matrix.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace std;

Spin_Density_Matrix::Spin_Density_Matrix( size_t size ) :
    m_size (size)
{
  m_entries.clear();
  size_t nentr = size * size;
  for( size_t i=0; i<nentr; i++ ) m_entries.push_back( 0. );
}
 
Spin_Density_Matrix::Spin_Density_Matrix( std::vector<Complex> entr )
{
  m_entries.clear();
  for( size_t i=0; i<entr.size(); i++ ) m_entries.push_back( entr[i] );
  switch( entr.size() ) {
    case 1 : m_size = 1; break;
    case 4 : m_size = 2; break;
    case 9 : m_size = 3; break;
    case 16: m_size = 4; break;
    case 25: m_size = 5; break;
    default: m_size = 0;
  }
}
 

void Spin_Density_Matrix::SetNoCorrelation()
{
  for( size_t i=0; i<m_entries.size(); i++ ) 
    m_entries[i] = Complex(0.,0.);
  for( size_t i=0; i<m_size; i++ ) 
    m_entries[ (m_size+1) * i ] = Complex( 1./m_size, 0. );
  return;
}

void Spin_Density_Matrix::SetUnitMatrix()
{
  for( size_t i=0; i<m_entries.size(); i++ ) 
    m_entries[i] = Complex(0.,0.);
  for( size_t i=0; i<m_size; i++ ) 
    m_entries[ (m_size+1) * i ] = Complex( 1., 0. );
  return;
}

Complex Spin_Density_Matrix::Trace()
{
  Complex sum(0.,0.);
  for( size_t i=0; i<m_size; i++ ) 
    sum += m_entries[ (m_size+1) * i ];
  return sum;
}

void Spin_Density_Matrix::Normalise()
{
  Complex tr = Trace();
  for( size_t i=0; i<m_entries.size(); i++ ) 
    m_entries[i] /= tr;
}


Spin_Density_Matrix& Spin_Density_Matrix::operator+=(Spin_Density_Matrix SDM)
{
  if (m_entries.size() != SDM.m_entries.size()) {
    msg_Error()<<"ERROR in Spin_Density_Matrix::operator+=:"<<endl
	       <<"SDMs with unequal numbers ("<<m_entries.size()<<","<<SDM.m_entries.size()<<") of entries passed over."<<endl
	       <<"Operation will be ignored."<<endl;
    return *this;
  }
    
  for (size_t i=0; i<m_entries.size(); ++i) m_entries[i]+=SDM.m_entries[i];
  return *this;
}
 
Spin_Density_Matrix& Spin_Density_Matrix::operator*=(Complex factor)
{
  for (size_t i=0; i<m_entries.size(); ++i) m_entries[i]*=factor;
  return *this;
}

Complex& Spin_Density_Matrix::operator[](size_t entry)
{
  if (entry<m_entries.size()) {
    return m_entries[entry];
  } else {
    msg_Error()<<"ERROR in Spin_Density_Matrix::operator[]:"<<endl
	       <<"Rquest for entry "<<entry<<" exceeds length of data."<<endl;
    abort();
  }
}

Complex& Spin_Density_Matrix::operator()(size_t a, size_t b)
{
  return (*this)[m_size*b+a];
}
 
void Spin_Density_Matrix::Print()
{
  for( unsigned int i=0; i<m_size; i++ ) {
    for( unsigned int j=0; j<m_size; j++ ) {
      cout<<m_entries[ i*m_size+j ]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;
}
 
bool Spin_Density_Matrix::operator== (const Spin_Density_Matrix & sdm )
{
  if( m_entries.size() != sdm.m_entries.size() ) return false;
  for( size_t i=0; i<m_entries.size(); ++i ) {
    if( !ATOOLS::IsEqual(m_entries[i],sdm.m_entries[i]) ) return false;
  }
  return true;
}
