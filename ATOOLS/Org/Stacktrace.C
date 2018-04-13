#include "ATOOLS/Org/Stacktrace.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"

#include <iomanip>

#define USING_Stack_Trace
#ifndef __USE_GNU
#ifdef __GNUC__
#define __USE_GNU
#ifdef ARCH_DARWIN
#undef USING_Stack_Trace
#endif
#else 
#undef USING_Stack_Trace
#endif
#endif

#ifdef USING_Stack_Trace
#include <execinfo.h>
#include <dlfcn.h>
#define MAX_BACKTRACE_DEPTH 128
#endif

void ATOOLS::GenerateStackTrace(std::ostream &ostr,
				const bool endline,
				const std::string &comment)
{
#ifdef USING_Stack_Trace
    ostr<<comment<<om::bold<<"GenerateStackTrace(..): "
	<<om::reset<<om::blue<<"Generating stack trace "
	<<om::reset<<om::bold<<"\n{"<<om::reset<<std::endl;
    // adapted from root version 3.10 TUnixSystem.cxx
    void *trace[MAX_BACKTRACE_DEPTH];
    int depth=backtrace(trace,MAX_BACKTRACE_DEPTH);
    for (int n=0; n<depth;++n) {
      unsigned long addr=(unsigned long)trace[n];
      Dl_info info;
      if (dladdr(trace[n],&info) && info.dli_fname && info.dli_fname[0]) {
	unsigned long symaddr=(unsigned long)info.dli_saddr;
	if (symaddr==(unsigned long)NULL) continue;
	const char *symname=info.dli_sname;
	if (!info.dli_sname || !info.dli_sname[0]) symname="<unknown function>";
	// if (!msg->LevelIsDebugging()) {
	//   if (std::string(symname).find
	//       ("Exception_Handler")!=std::string::npos ||
	//       std::string(symname).find
	//       ("HandleSignal")!=std::string::npos) continue;
	// }
	std::string linfo;
	unsigned long libaddr=(unsigned long)info.dli_fbase;
	unsigned long offset=(addr>=libaddr)?addr-libaddr:libaddr-libaddr;
	char cmd[4096];
	sprintf(cmd,"addr2line -se %s 0x%016lx 2>/dev/null",
		info.dli_fname,offset);
	if (FILE *pf=popen(cmd,"r")) {
	  char buf[2048];
	  if (fgets(buf,2048,pf)) {
	    linfo=buf;
	    linfo=linfo.substr(0,linfo.length()-1);
	  }
	  if (linfo=="??:0") {
	    pclose(pf);
	    sprintf(cmd,"addr2line -se %s 0x%016lx 2>/dev/null",
		    info.dli_fname,addr);
	    pf=popen(cmd,"r");
	    if (fgets(buf,2048,pf)) {
	      linfo=buf;
	      linfo=linfo.substr(0,linfo.length()-1);
	    }
	    if (linfo=="??:0") linfo="";
	  }
	  pclose(pf);
	}
	ostr<<comment<<"  "<<std::setiosflags(std::ios::left)
	    <<std::setw(15)<<trace[n]<<std::dec
	    <<" in '"<<om::red<<Demangle(symname)<<om::reset<<"' ";
	if (linfo!="") ostr<<"("<<om::lblue<<linfo<<om::reset<<")";
	ostr<<"\n";
	if (msg->LevelIsDebugging()) ostr<<"                from '"<<
				       om::brown<<info.dli_fname<<om::reset<<"'\n";
	ostr<<std::flush;
	if (std::string(info.dli_sname)=="main") break;
      } 
      else {
	ostr<<comment<<"   "<<addr<<" in   <unknown function>"<<std::endl;
      }
    }
    ostr<<comment<<om::bold<<"}"<<om::reset;
    if (endline) ostr<<std::endl;
#endif
}
