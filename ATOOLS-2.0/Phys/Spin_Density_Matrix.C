#include "Spin_Density_Matrix.H"
#include "Blob.H"
#include "Message.H"

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

double Spin_Density_Matrix::Contract( vector<Complex> * ampls, vector<int> * ind )
{
  // find the position of the mother
  int pos_mother (0);           // position of the decayer
  for( unsigned int i=0; i<ind->size(); ++i ) {
    if( (*ind)[i]==0 ) { pos_mother = i; break; }
  }
  // contract with SDM
  int n = (1<<ind->size());        // = 2^{number of indices}
  Complex ret (0., 0.);
  Complex contr_da (0.,0.);         // value over contracted daughters
  int l0  (0);                      // lambda_0
  int l0p (0);                      // lambda_0'
  int h1_a (0);                     // helicity combination for M
  int h1_b (0);                     // helicity combination for M*
  int ref (0);                      // a reference for calculating h1_a,b
  for( int h0=0; h0<4; ++h0 ) {     // for all combinations (lambda_0,lambda_0')
    l0  = h0>>1;
    l0p = h0 & 1;
    contr_da = Complex(0.,0.);        
    for( int h1=0; h1<(n>>1); ++h1 ) { // for all daughters helicity comb.
      // "insert" mother helicity
      ref = n - (1<<pos_mother);    
      h1_a = ((h1 & ref)<<1) + (l0<<pos_mother)  + (h1 & (~ref));
      h1_b = ((h1 & ref)<<1) + (l0p<<pos_mother) + (h1 & (~ref));
      // sum up
      contr_da += (*ampls)[h1_a] * conj( (*ampls)[h1_b] );
    }
    ret += m_entries[h0] * contr_da;
  }
  return ret.real();
} // the result is exactly the same as using SCT Contract Methods

Spin_Density_Matrix& Spin_Density_Matrix::operator+=(Spin_Density_Matrix SDM)
{
  if (m_entries.size() != SDM.m_entries.size()) {
    msg.Error()<<"ERROR in Spin_Density_Matrix::operator+=:"<<endl
	       <<"SDMs with unequal number of entries passed over."<<endl
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
    msg.Error()<<"Error in Spin_Density_Matrix::operator[]:"<<endl
	       <<"Rquest for entry exceeds length of data."<<endl;
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
 
std::ostream& operator<<(std::ostream &ostr, Spin_Density_Matrix &sdm)
{
    ostr<<"{";
    ostr<<sdm[0]<<",";   
    ostr<<sdm[1]<<",";      
    ostr<<sdm[2]<<",";   
    ostr<<sdm[3]<<"}";
  return ostr;
}

bool Spin_Density_Matrix::operator== (const Spin_Density_Matrix & sdm )
{
  if( m_entries.size() != sdm.m_entries.size() ) return false;
  for( size_t i=0; i<m_entries.size(); ++i ) {
    if( !ATOOLS::IsEqual(m_entries[i],sdm.m_entries[i]) ) return false;
  }
  return true;
}
