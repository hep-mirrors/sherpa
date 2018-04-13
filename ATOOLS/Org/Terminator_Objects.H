#ifndef ATOOLS_Org_Terminator_Objects_H
#define ATOOLS_Org_Terminator_Objects_H

#include <vector>
#include <string>

namespace ATOOLS {

  class Exception;

  class Terminator_Object {
    
  protected:
    
    virtual bool ReadInStatus(const std::string &path)
    {return true;}
    
    virtual void PrepareTerminate() {}
    
    friend class Terminator_Object_Handler;
    
  public:

    virtual ~Terminator_Object() {}

  };

  class Terminator_Object_Handler {

  private:

    bool m_noremove;
    
    unsigned int m_nbus, m_nsegv;

    typedef void (*Terminator_Function)(void);

    std::vector<Terminator_Function> m_terminatorfunctions;
    std::vector<Terminator_Object*>  m_terminatorobjects;

    friend class Exception;
    
  public:

    Terminator_Object_Handler();

    ~Terminator_Object_Handler();

    bool ReadInStatus(const std::string &path);
    void PrepareTerminate();
    
    void HandleSignal(int signal);
    void Terminate(unsigned int excode);

    void AddTerminatorFunction(void (*function)(void));
    void AddTerminatorObject(Terminator_Object *const object);

    void RemoveTerminatorFunction(void (*function)(void));
    void RemoveTerminatorObject(Terminator_Object *const object);
    
  };

  extern Terminator_Object_Handler *exh;

  void HandleSignal(int signal); 
  void Terminate();
  
}

#endif
