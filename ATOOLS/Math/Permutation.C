#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

Permutation::Permutation(int n) : m_n(n)     
{
//   std::cout<<"Init Permutation: "<<m_n<<std::endl;
  p_per = new int[m_n];
  p_snum = new int[m_n];
  p_an  = new int*[m_n];
  for (int i=0;i<m_n;i++) p_an[i]= new int[m_n];
  for (int i=0;i<m_n;i++) p_an[0][i]=i;
}

Permutation::~Permutation()
{
  for (int i=0;i<m_n;i++) delete[] p_an[i];
  delete[] p_an;
  delete[] p_snum;
  delete[] p_per;
}

int Permutation::MaxNumber()
{
  int pn=1;
  for(int i=2;i<=m_n;i++) pn*=i;
//   std::cout<<"call MaxNumber: "<<pn<<std::endl;
  return pn;
}

int* Permutation::Get(int n) 
{
  int x = MaxNumber();
  for(int i=m_n;i>0;i--) {
    n = n%x;
    x/=i;
    p_snum[m_n-i]=n/x;
//       std::cout<<p_snum[m_n-i]<<" "<<m_n-i<<"! ";
  }
//     std::cout<<std::endl;
  p_per[0]=p_snum[0];
  for(int i=1;i<m_n;i++) {
    int j=0; int k=0;
    while (j<m_n-i) {
      if (p_an[i-1][k]==p_per[i-1]) k++;
      p_an[i][j]=p_an[i-1][k];
      j++;k++;
    }
//     std::cout<<i<<": ";
//     for(int l=0;l<m_n-i;l++)std::cout<<p_an[i][l]<<" ";
//     std::cout<<std::endl;
    p_per[i]=p_an[i][p_snum[i]];
  }
  return p_per;
}
