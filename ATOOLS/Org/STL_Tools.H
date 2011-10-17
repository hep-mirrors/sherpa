#ifndef ATOOLS_Org_STL_Tools_H
#define ATOOLS_Org_STL_Tools_H

#include "ATOOLS/Org/CXXFLAGS.H"

#include <iostream>
#include <vector>
#ifdef USING__Threading
#include <pthread.h>
#endif

namespace std {

  template <typename __Tp>
  std::ostream &operator<<(std::ostream &str,const std::vector<__Tp> &v);

}

namespace ATOOLS {

  class String_Sort {
  private:

    std::less<std::string> m_less;

  public:

    bool operator()(const std::string &a,const std::string &b) const;
    
  };// end of class String_Sort

  std::vector<int> ID(size_t id);

  size_t IdCount(size_t id);

  template <class __Tp>
  class AutoDelete_Vector: public std::vector<__Tp*> {
  private:

#ifdef USING__Threading
    pthread_mutex_t m_mtx;
#endif

  public:

    AutoDelete_Vector()
    {
#ifdef USING__Threading
      pthread_mutex_init(&m_mtx,NULL);
#endif
    }

    virtual ~AutoDelete_Vector() 
    {
      while (!this->empty()) {
	delete this->back();
	this->pop_back();
      }
#ifdef USING__Threading
      pthread_mutex_destroy(&m_mtx);
#endif
    }

    inline void MtxLock()
    {
#ifdef USING__Threading
      pthread_mutex_lock(&m_mtx);
#endif
    }

    inline void MtxUnLock()
    {
#ifdef USING__Threading
      pthread_mutex_unlock(&m_mtx);
#endif
    }

  };// end of class AutoDelete_Vector

}

#endif

