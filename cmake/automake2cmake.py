import sys
import os
import re
########################################################################
def get_suffixes(text):
    if (text=="xbae"):
      return ["_static"]
    return ["_static", ""]
def if_install_library(text):
    if (text=="xbae"):
      return False
    return True
def link_static(text):
    if (text=="pythia"):
        return "pythia"
    if (text=="geant321"):
        return "geant"
    if (text=="herwig"):
        return "herwig"
    if (text=="kernlib"):
        return "kernlib-${SHIFTSUFFIX}"
    if (text=="packlib"):
        return "packlib-${SHIFTSUFFIX}"
    return "no"
def output_name_static(text):
    if (text=="pythia"):
        return "pythia6205"
    if (text=="higz"):
        return "grafX11"
    if (text=="code_motif"):
        return "packlib-lesstif"
    if (text=="paw_motif"):
        return "pawlib-lesstif"
    if (text=="herwig"):
        return "herwig59"
    if (text=="photos"):
        return "photos202"
    if (text=="pdf"):
        return "pdflib804"
    if (text=="isajet"):
        return "isajet758"
    if (text=="jetset"):
        return "jetset74"
    if (text=="lepto63"):
        return "lepto651"
    return text
########################################################################    
def output_name(text):
    if (text=="pythia"):
        return "pythia6205"
    if (text=="higz"):
        return "grafX11"
    if (text=="code_motif"):
        return "packlib-lesstif"
    if (text=="paw_motif"):
        return "pawlib-lesstif"
    if (text=="herwig"):
        return "herwig59"
    if (text=="photos"):
        return "photos202"
    if (text=="pdf"):
        return "pdflib804"
    if (text=="isajet"):
        return "isajet758"
    if (text=="jetset"):
        return "jetset74" 
    if (text=="lepto63"):
        return "lepto651"               
    return text
########################################################################    
def get_simple_so_version(text):
    if (text=="jetset"):
        return "3_${COMPSUFFIX}"
    if (text=="fritiof"):
        return "1_${COMPSUFFIX}"
    if (text=="ariadne"):
        return "1_${COMPSUFFIX}"
    if (text=="higz"):
        return "1_${COMPSUFFIX}"
    if (text=="code_motif"):
        return "1_${COMPSUFFIX}"
    if (text=="paw_motif"):
        return "3_${COMPSUFFIX}"
    if (text=="eurodec"):
        return "1_${COMPSUFFIX}"
    if (text=="graflib"):
        return "1_${COMPSUFFIX}"
    if (text=="kernlib"):
        return "1_${COMPSUFFIX}"
    if (text=="packlib"):
        return "1_${COMPSUFFIX}"
    if (text=="photos"):
        return "1_${COMPSUFFIX}"
    if (text=="isajet"):
        return "3_${COMPSUFFIX}"        
    if (text=="pythia"):
        return "3_${COMPSUFFIX}"  
    if (text=="lepto63"):
        return "3_${COMPSUFFIX}"
    return "2_${COMPSUFFIX}"
########################################################################
def get_full_so_version(text):
     return get_simple_so_version(text)+".${SHERPA_VERSION_MAJOR}"
########################################################################
def comment_remover(text):
    def replacer(match):
        s = match.group(0)
        if s.startswith('/'):
            return " " # note: a space and not an empty string
        else:
            return s
    pattern = re.compile(r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',re.DOTALL | re.MULTILINE)
    return re.sub(pattern, replacer, text)
########################################################################
def write_header(f,ldirs):
   f.write("""
########################################################################
#
#  Automatically or semiautomaticall generated, do not edit.
#
########################################################################
# The following input was used
""")
   for a in ldirs:
     f.write("# "+a +"\n")
   f.write("""
########################################################################
""")
########################################################################
def get_if_condition(a):
  if len(a)>0:
    if re.match(r'.*packlib_cspack_sysreq_.*',a):
      return [
"""
#Inserted by get_if_condition()->
if ( (SHERPA_VAXVMS OR SHERPA_UNIX)  AND  ( (NOT SHERPA_SHIFT) AND ( NOT SHERPA_WINNT) ) )
""",
"""
else()
  set(packlib_cspack_sysreq_CSRC )
endif()
#Inserted by get_if_condition()<-
"""]
  return ["",""]
########################################################################
def transform_pilots(argv):
    path_to_file=argv
    file_contents=''
    with open(path_to_file) as f:
      file_contents = f.readlines()
    newlist =file_contents
    newlist = [x.replace("\t"," ") for x in newlist]
    newlist = [x.replace("    "," ") for x in newlist]
    newlist = [x.replace("   "," ") for x in newlist]
    newlist = [x.replace("  "," ") for x in newlist]
    newlist = [x.replace("# ","#") for x in newlist]
    newlist = [x.replace("# ","#") for x in newlist]
    newlist = [x.strip() for x in newlist]
    newlist = filter(lambda st: st != '' , newlist)
    #newlist = [  "\n"+x.replace("#define ","set(")+"1)\n  add_definitions(-D"+x.replace("#define ","")+")\n"  if re.match(r'^#define .*',x)  else x for x in newlist]
    newlist = [  "\n"+x.replace("#define ","set(")+" 1)\n"  if re.match(r'^#define .*',x)  else x for x in newlist]
    newlist = [  x.replace("#ifndef ","if ( NOT (")+") )\n" if re.match(r'^#ifndef .*',x)  else x for x in newlist]
    newlist = [  x.replace("#ifdef ","if ( (")+") )\n" if re.match(r'^#ifdef .*',x)  else x for x in newlist]
    newlist = [  x.replace("#endif","endif()") if re.match(r'^#endif.*',x)  else x for x in newlist]
    newlist = [  x.replace("&&"," AND ")  for x in newlist]
    newlist = [  x.replace("||"," OR ")  for x in newlist]

    newlist = [x.strip() for x in newlist]
    #newlist = filter(lambda st: st != '' , newlist)
    newlist = [x.replace("set(","  set(") if not re.match(r'^  .*',x)  else x for x in newlist]
    newlist = [  x.replace("#if","if (")+" )" if re.match(r'^#if.*',x)  else x for x in newlist]
    newlist = [  x.replace("!defined","NOT defined") for x in newlist]
    newlist = [  x.replace("*-","#") for x in newlist]
    newlist = [  x.replace("NOT defined(","(NOT") for x in newlist]
    newlist = [  x.replace("defined(","(") for x in newlist]
    
    
    #for ai in file_contents:
    #   print(ai)
    #for ai in newlist:
    #   print(ai)
    return newlist


def transform_imake_source(argv, dbg):
    path_to_file=argv
    with open(path_to_file) as f:
      file_contents = f.readlines()
    my_lst_str = ''.join(file_contents)
    x=comment_remover(my_lst_str)

    filtered = x.split('\n')
    filtered= filter(lambda st: st != '\n' , filtered)
    filtered= filter(lambda st: st != '' , filtered)
    filtered= filter(lambda st: st != ' ' , filtered)
#    filtered= filter(lambda st: not re.match(".*MotifDependantMakeVar.*",st), filtered)
#    filtered= filter(lambda st: not re.match(".*SubdirLibraryTarget.*",st), filtered)
#    filtered= filter(lambda st: not re.match(".*SHERPACcProgramTarget.*",st), filtered)
#    filtered= filter(lambda st: not re.match(".*TopOfPackage.*",st), filtered)
#    filtered= filter(lambda st: not re.match(".*NormalFortranProgramTarget.*",st), filtered)
#    filtered= filter(lambda st: not re.match("^NormalProgramTarget.*",st), filtered)
#    filtered= filter(lambda st: not re.match("^InstallNonExecFileTarget.*",st), filtered)    
#    filtered= filter(lambda st: not re.match("^InstallSharedLibrary.*",st), filtered)
#    filtered= filter(lambda st: not re.match("^InstallLibrary.*",st), filtered)
#    filtered= filter(lambda st: not re.match("^AllTarget.*",st), filtered)
#    filtered= filter(lambda st: not re.match("^TestTarget.*",st), filtered)
#    filtered= filter(lambda st: not re.match("^DoIncludePackage.*",st), filtered)
#    filtered= filter(lambda st: not re.match("^SubdirLibraryTarget.*",st), filtered)    
#    filtered= filter(lambda st: not re.match(".*LinkFileFromDir.*",st), filtered)
#    filtered= filter(lambda st: not re.match(".*FortranCmd.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*VMS_OPT_FILES.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*IMAKE_INCLUDES.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*CERNDEFINES.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*CLIBS.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*nstallProgram.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*eedTcpipLib.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*NeedSysexe.*",st), filtered)
    #filtered= filter(lambda st: not re.match(".*CCOPTIONS.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*\+Z \+DA1.*",st), filtered)
    filtered= filter(lambda st: not re.match("^gxint321.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*cd \$\(.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*gxint\.f.*",st), filtered)
    filtered= filter(lambda st: not re.match("^FC=mpxlf.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*install\.lib.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*clean::.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*Makefile::.*",st), filtered)
    ###filtered= filter(lambda st: not re.match(".*LibraryTargetName.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*RemoveFile.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*CopyFile.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*geant321_parallel.*",st), filtered)
    ###filtered= filter(lambda st: not re.match(".*CppFileTarget.*",st), filtered)
    #filtered= filter(lambda st: not re.match("^Install.*",st), filtered)
    #filtered= filter(lambda st: not re.match("^EXTRA_.*",st), filtered)
    #filtered= filter(lambda st: not re.match("^#include.*",st), filtered)
    filtered= filter(lambda st: not re.match("^IMAKE_DEFINES.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*FORTRANSAVEOPTION.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*FORTRAN.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]FORTRAN.*",st), filtered)
    filtered= filter(lambda st: not re.match("^FORTRANOPTIONS.*",st), filtered)
    filtered= filter(lambda st: not re.match("^SpecialFortranLibObjectRule.*",st), filtered)
    filtered= filter(lambda st: not re.match("^SpecialFortranObjectRule.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*SpecialFortranObjectRule.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*DefinePackage.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*INCLUDES.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*SHERPAFortran.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*PACKAGE.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*FDEBUGFLAGS.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*EXTRA_DEFINES.*",st), filtered)
    filtered= filter(lambda st: not re.match("^[:blank:]*CDEBUGFLAGS.*",st), filtered)
    filtered= filter(lambda st: not re.match("^#define.*",st), filtered)

    filtered = [x.replace(" : ",": ") for x in filtered]
    filtered = [x.replace("\t"," ") for x in filtered]
    filtered = [x.replace("="," = ") for x in filtered]
    filtered = [x.replace("\\ @@\\","") for x in filtered]
    filtered = [x.replace("\\@@\\","") for x in filtered]
    filtered = [x.replace("        "," ") for x in filtered]
    filtered = [x.replace("      "," ") for x in filtered]
    filtered = [x.replace("    "," ") for x in filtered]
    filtered = [x.replace("  "," ") for x in filtered]
    filtered = [x.replace("  \n","\n") for x in filtered]
    filtered = [x.replace("  \n","\n") for x in filtered]
    filtered = [x.replace("  \n","\n") for x in filtered]
    filtered = [x.replace("  \n","\n") for x in filtered]
    filtered = [x.replace("  \n","\n") for x in filtered]
    filtered = [x.replace("  "," ") for x in filtered]
    filtered = [x.strip(' ') for x in filtered]
    
    newlist=[]
    temp=''
    for a in filtered:

     if re.match("^#.*",a):
        newlist.append(temp)
        newlist.append(a)
        temp=''
     else:
      if re.match("^SOURCES.*",a):
        newlist.append(temp)
        temp=a
      else:
        temp+=' '
        temp+=a
     
     LISTNAME=path_to_file.replace("//","/").replace("/","_")
     LISTNAME=LISTNAME.replace("_Makefile.am","")
     LISTNAME=LISTNAME.replace(".","0")#OH!
     
    incfiles = [ x  for x in newlist if re.match(r'#include.*',x)]
    incfiles =  [ x.replace("#include","") for x in incfiles]
    incfiles =  [ x.replace("\"","") for x in incfiles]
    incfiles =  [ x.replace(" ","") for x in incfiles]
    newlist.append(temp)
    newlist = filter(lambda st: st != '' , newlist)
    newlist = [x.replace(" \\ "," ") for x in newlist]
    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*pkglib_LTLIBRARIES.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*CppTarget.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*FDEBUGFLAGS.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*CDEBUGFLAGS.*',x) else x for x in newlist]
    newlist = [ ("\n#"+x) if re.match(r'^[:blank:]*[^#]*dist_bin_SCRIPTS.*',x) else x for x in newlist]
    newlist = [ ("\n#"+x) if re.match(r'^[:blank:]*[^#]*SUBDIRS.*',x) else x for x in newlist]
#    newlist = [ ("\n#"+x) if re.match(r'^[:blank:]*[^#]*localinc.*',x) else x for x in newlist]

#    newlist = [x.replace("(PACKAGE_LIB)","{PACKAGE_LIB}") for x in newlist]
#    newlist = [x.replace("MOTIF_","") for x in newlist]
    newlist = [x.replace("$(GZIPSTREAMSOURCES)","${GZIPSTREAMSOURCES}") for x in newlist]
    newlist = [x.replace("GITTAG","\n#GITTAG") for x in newlist]
    newlist = [x.replace("bin_PROGRAMS","\n#bin_PROGRAMS") for x in newlist]
    newlist = [x.replace("localinc","\n#localinc") for x in newlist]
    newlist = [x.replace("Git_Info.C","") for x in newlist]
#    newlist = [x.replace(" : =",": =") for x in newlist]
#    newlist = [x.replace("_F+","_F +") for x in newlist]
#    newlist = [x.replace("_C+","_C +") for x in newlist]

    #newlist = [ ("#ORIGINAL "+x+"\n"+"\n"+x) if re.match(r'.*defined.*',x) else x for x in newlist]
    l=[]
    for x in newlist:
      if re.match(r'.*defined.*',x):
        l.append("#ORIGINAL "+x)
      l.append(x)

    newlist=l    
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("! defined(","!defined(") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#if (defined(SHERPA_LINUX) && (!defined(SHERPA_GFORTRAN)))","if (SHERPA_LINUX AND NOT SHERPA_GFORTRAN)") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#if defined(SHERPA_DECS) || (defined(SHERPA_LINUX) && !defined(SHERPA_PPC)) || defined(SHERPA_WINNT)","if (SHERPA_DECS OR (SHERPA_LINUX AND NOT SHERPA_PPC) OR SHERPA_WINNT)") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#if (!defined(SHERPA_ASSEMB) && defined(SHERPA_OLD))","if ( NOT SHERPA_ASSEMB AND SHERPA_OLD)") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#if (!defined(SHERPA_NTC)) && (!defined(SHERPA_X11))","if (NOT SHERPA_NTC AND NOT SHERPA_X11)") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#if !defined(","if (NOT ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("F_ARCHITECTURE = ","set(F_ARCHITECTURE ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("(F_ARCHITECTURE)","{F_ARCHITECTURE}  ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#else","else()") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#if defined","if ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#ifdef ","if (") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#ifndef ","if (NOT ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#endif","endif()") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("#endif","endif()") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace(") || defined("," OR ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace(") && !defined("," AND NOT ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace(") && defined("," AND ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace(") && defined("," AND ") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace(") && (!defined(SHERPA_PHIGS)) && (!defined(SHERPA_MSDOS)"," AND NOT SHERPA_PHIGS AND NOT SHERPA_MSDOS") for x in newlist] 
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("(SHERPA_UNIX) && (!defined(SHERPA_WINNT))","(SHERPA_UNIX AND NOT SHERPA_WINNT)") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("if (SHERPA_SHL) && ( defined(SHERPA_SUN OR SHERPA_SGI OR SHERPA_IBMRT OR SHERPA_QMVAOS OR SHERPA_LINUX) )","if (SHERPA_SHL AND ( SHERPA_SUN OR SHERPA_SGI OR SHERPA_IBMRT OR SHERPA_QMVAOS OR SHERPA_LINUX) )") for x in newlist]
#    newlist = [x if re.match(r'^#ORIGINAL .*$',x) else x.replace("SHERPA_LNX AND NOT SHERPA_QMLXIA64) && (!defined(SHERPA_GFORTRAN)"," SHERPA_LNX AND NOT SHERPA_QMLXIA64 AND NOT SHERPA_GFORTRAN") for x in newlist]
    newlist = [x.replace("_SOURCES =","\n  set("+LISTNAME+"_SOURCES") for x in newlist]
#    newlist = [x.replace("SRCS_CDF: = $(SRCS_CDF)","  list(APPEND "+LISTNAME+"_CDFSRC ") for x in newlist]
#    newlist = [x.replace("SRCS_C =","  set("+LISTNAME+"_CSRC") for x in newlist]
#    newlist = [x.replace("SRCS_C: = $(SRCS_C)","  list(APPEND "+LISTNAME+"_CSRC ") for x in newlist]
#    newlist = [x.replace("SRCS_C + = ","  list(APPEND "+LISTNAME+"_CSRC ") for x in newlist]
#    newlist = [x.replace("SRCS_F =","  set("+LISTNAME+"_FSRC") for x in newlist]
#    newlist = [x.replace("SRCS_F: = $(SRCS_F)","  list(APPEND "+LISTNAME+"_FSRC ") for x in newlist]
#    newlist = [x.replace("SRCS_F + = ","  list(APPEND "+LISTNAME+"_FSRC ") for x in newlist]
#    newlist = [x.replace("SRCS_S =","  set("+LISTNAME+"_SSRC") for x in newlist]
#    newlist = [x.replace("SRCS_S: = $(SRCS_S)","  list(APPEND "+LISTNAME+"_SSRC ") for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*EXTRA_DEFINES.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*EXTRA_INCLUDES.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*EXTRA_LDOPTIONS.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*test:.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*CCOPTIONS.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*CERNDEFINES: =.*',x) else x for x in newlist]    
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*DEFINES: =.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*SpecialCObjectRule.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*CCOPTIONS + =.*',x) else x for x in newlist]    
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*FCLDOPTIONS.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*InstallProgram.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*InstallScr.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*MotifD.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*NeedTcpipLib.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*SpecialO.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*SpecialFortran.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*SQUEEZE.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*FORTRANSAVEOPTION.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*FDEBUGFLAGS.*',x) else x for x in newlist]
#    newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*FDEBUGFLAGS.*',x) else x for x in newlist]
#    newlist = [x.replace(") || (defined(SHERPA_LINUX) && (!defined(SHERPA_PPC))"," OR (SHERPA_LINUX AND (NOT SHERPA_PPC))") for x in newlist] 
#    newlist = [x.replace("#if ( defined(SHERPA_UNIX OR SHERPA_VAXVMS) ) && (!defined(SHERPA_NOCIO))","if ( SHERPA_UNIX OR SHERPA_VAXVMS  AND (NOT SHERPA_NOCIO))") for x in newlist] 
#    newlist = [x.replace("(SRCS_F)","{SRCS_F}") for x in newlist]
    newlist=[ x.replace(" \n","\n") for x in newlist]
    gnlist=[]
    for a in newlist:
      gnlist+=a.split('\n')
    newlist=gnlist
    fin=["#"+argv, " "] 
    for a in newlist:
      if (len(a) > 0):
       if a[-1] != ')' and a[0] !='#':
        fin.append( (a+")"))
       else:
        fin.append(a)
    dirf=path_to_file.replace("Makefile.am","")
    sr=[]
    FSRC=0

    for a in fin:
     if re.match(".*_SOURCES.*",a):
       FSRC=1
    print(fin,FSRC)
    
    diractual="../"+dirf
    diractual=diractual.replace("//","/")
    if FSRC==1:
      fin.append("  list(TRANSFORM "+LISTNAME+"_SOURCES PREPEND \"${CMAKE_CURRENT_SOURCE_DIR}/"+diractual+"\")")
      sr.append("${"+LISTNAME+"_SOURCES}")
    lev=0
    fin.append(" ")
    return [[LISTNAME],fin,sr,incfiles]
########################################################################
def write_to_file_with_breaks(f, lin, n):
   p=0
   while p<len(lin) and p>=0:
      np=lin.find(" ",p+n)
      f.write(lin[p:np])
      f.write("\n")
      p=np
########################################################################
def create_library(ldirsI,lname,includes,installincludes,linklibs=[], cdff=[],pat=[ "\"*makefile*\"", "\"*\.c\"" ]):
   ldirs=ldirsI
   ldirs.sort()
   includes.sort()
   installincludes.sort()
   cdff.sort()
   t=[]
   f = open(lname+"/CMakeLists.txt", "w")
   write_header(f,ldirs)
#   f.write("set_package_flags("+lname+")\n")
   for a in ldirs:
     x=transform_imake_source(lname+"/"+a,0)
     t+=x[2]
     if len(x[3])>0:
        plt=x[3][0]
        f.write("#The original Imake file below included files:"+plt+"\n#Those were NOT processed.\n")
        U=[]
        #print(lname+"/"+lname+"/"+plt)
        #if os.path.exists(lname+"/"+lname+"/"+plt):
        #  U=transform_pilots(lname+"/"+lname+"/"+plt)
        #print(U)
        #for  tf in U:
        #  f.write(tf+"\n")
     
     for a in x[1]: 
          if a[0] =='#':
            f.write(a+"\n")
          else:
           if re.match(r'.*set.*',a):
              ifcondition = get_if_condition(a)
              if len(ifcondition[0])>0:
                f.write(ifcondition[0]+"\n")
              write_to_file_with_breaks(f, a+"\n", 200)
              if len(ifcondition[1])>0:
                f.write(ifcondition[1]+"\n")
           else:
             f.write(a+"\n")
   f.write("set("+lname+"_esources )\n")
   for cdf in cdff:
      f.write("cdf_compile(${CMAKE_CURRENT_SOURCE_DIR}/"+cdf+" ${CMAKE_CURRENT_BINARY_DIR}/"+cdf.split("/")[-1]+".c)\n")
      f.write("list(APPEND "+lname+"_esources ${CMAKE_CURRENT_BINARY_DIR}/"+cdf.split("/")[-1]+".c)\n")
   suffixes=get_suffixes(lname)
   for suff in suffixes:
     if suff=="": 
       f.write("if (SHERPA_BUILD_SHARED)\n")
       f.write("add_library("+lname+" SHARED ${"+lname+"_esources}\n")
       for z in t: f.write("                             "+z + " \n")
       f.write(")\n")
     if suff=="_static": 
       f.write("if (SHERPA_BUILD_STATIC)\n")
       f.write("add_library("+lname+"_static STATIC ${"+lname+"_esources}\n")
       for z in t: f.write("                             "+z + " \n")
       f.write(")\n")
     f.write("target_include_directories("+lname+suff+" PRIVATE ${PROJECT_SOURCE_DIR}/include)\n")
     f.write("target_include_directories("+lname+suff+" PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})\n")
     if os.path.exists(lname+"/"+lname+"/gen"):
       f.write("target_include_directories("+lname+suff+" PRIVATE \"${CMAKE_CURRENT_SOURCE_DIR}/"+lname+"/gen\"  )\n")
     if os.path.exists(lname+"/"+"gen"):
       f.write("target_include_directories("+lname+suff+" PRIVATE \"${CMAKE_CURRENT_SOURCE_DIR}/gen/\")\n")
     for inc in includes:
       if os.path.exists(lname+"/"+inc):
         u="/"+inc
         u=u.replace("//","/")
         f.write("target_include_directories("+lname+suff+" PRIVATE \"${CMAKE_CURRENT_SOURCE_DIR}"+u+"\")\n")
#     f.write("target_include_directories("+lname+suff+" PRIVATE ${FREETYPE_INCLUDE_DIRS})\n")
     for ll in linklibs:
       if ll=="packlib":
         ll=ll+suff
       f.write("target_link_libraries("+lname+suff+" PRIVATE "+ll+")\n")
     if if_install_library(lname):
       f.write("install(TARGETS "+lname+suff+" DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libs)\n")   
     if suff=="_static": 
       f.write("set_target_properties("+lname+"_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME "+output_name_static(lname)+")\n")
       if (link_static(lname)!="no"): 
         if if_install_library(lname):
           f.write("install_symlink(lib"+output_name_static(lname)+".a "+" ${CMAKE_INSTALL_LIBDIR}/lib"+link_static(lname)+".a)\n")
     if suff=="":
       f.write("set_target_properties("+lname+"        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME "+output_name(lname)+" SOVERSION "+get_full_so_version(lname)+")\n")
#       if if_install_library(lname):
#         f.write("install_symlink(lib"+output_name(lname)+".so."+get_full_so_version(lname)+" "+"${CMAKE_INSTALL_LIBDIR}/lib"+output_name(lname)+".so."+get_simple_so_version(lname)+")\n")
     f.write("endif()\n")
   for inc in installincludes:
     f.write("install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/"+inc+" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}  COMPONENT devel ")
     for pt in pat:
      f.write(" PATTERN "+ pt +" EXCLUDE ")
     f.write(")\n")
   f.close()
########################################################################
if __name__ == '__main__':

########################################################################
   ldirs =  "Channels/Makefile.am Decays/Makefile.am Enhance/Makefile.am Main/Makefile.am Makefile.am Process/Makefile.am Scales/Makefile.am Selectors/Makefile.am".split(" ")
   lname="PHASIC++"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)
########################################################################
   ldirs =  "Makefile.am Main/Makefile.am".split(" ")
   lname="BEAM"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)
########################################################################
   ldirs =  "Current_Library/Makefile.am ME_Library/Makefile.am Main/Makefile.am Makefile.am PS_Library/Makefile.am Run/Makefile.am".split(" ")
   ldirs =  "Current_Library/Makefile.am ME_Library/Makefile.am Main/Makefile.am Makefile.am PS_Library/Makefile.am".split(" ")
   lname="HADRONS++"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)
   
########################################################################
#   ldirs =  "Current_Library/Makefile.am ME_Library/Makefile.am Main/Makefile.am Makefile.am PS_Library/Makefile.am Run/Makefile.am".split(" ")
   ldirs =  "Decays/Makefile.am Formation/Makefile.am Main/Makefile.am Makefile.am Tools/Makefile.am".split(" ")
   lname="AHADIC++"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)


   ldirs =  "Amplitude/Makefile.am Amplitude/Zfunctions/Makefile.am DipoleSubtraction/Makefile.am Main/Makefile.am Makefile.am Phasespace/Makefile.am String/Makefile.am".split(" ")
   lname="AMEGIC++"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)


   ldirs =  "Main/Makefile.am Makefile.am Perturbative/Makefile.am Tools/Makefile.am".split(" ")
   lname="AMISIC++"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)

#  ./AHADIC++/Decays/Makefile.am
#./AHADIC++/Formation/Makefile.am
#./AHADIC++/Main/Makefile.am
#./AHADIC++/Makefile.am
#./AHADIC++/Tools/Makefile.am
 
   
########################################################################
   ldirs =  "Amplitude/Makefile.am Main/Makefile.am Makefile.am Phasespace/Makefile.am".split(" ")
   lname="COMIX"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)
########################################################################

   ldirs =  "Makefile.am Math/Makefile.am Org/Makefile.am Phys/Makefile.am YAML/Makefile.am YAML/yaml-cpp/Makefile.am YAML/yaml-cpp/node/Makefile.am YAML/yaml-cpp/node/detail/Makefile.am".split()
   lname="ATOOLS"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)

