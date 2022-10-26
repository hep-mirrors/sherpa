import sys
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
    return output_name(text)
########################################################################    
#"Electron", "GRV", "GRS", "SASG", "SAL", "Main", "NNPDF", "LHAPDF"
def output_name(text):
    if (text=="Electron"):
        return "PDFESherpa"
    if (text=="GRV"):
        return "GRVSherpa"
    if (text=="GRS"):
        return "GRSSherpa"
    if (text=="SAL"):
        return "SALSherpa"
    if (text=="SASG"):
        return "SASGSherpa"
    if (text=="NNPDF"):
        return "NNPDFSherpa"
    if (text=="Main"):
        return "PDF"
    if (text=="AmisicMain"):
        return "Amisic"
    if (text=="ComixMain"):
        return "Comix"

    if (text=="SherpaMain"):
        return "ModelMain"

    if (text=="SherpaUFO"):
        return "ModelUFO"


    if (text=="HadronsCurrent_Library"):
        return "HadronsCurrents"
    if (text=="HadronsME_Library"):
        return "HadronsMEs"
    if (text=="HadronsPS_Library"):
        return "HadronsPSs"

    if (text=="ShrimpsBeam_Remnants"):
        return "ShrimpsBeamRemnants"


    if (text=="ShrimpsCross_Sections"):
        return "ShrimpsXsecs"

    if (text=="ShrimpsEvent_Generation"):
        return "ShrimpsEvents"

    if (text=="EXTAMP"):
        return "ExtAmp"
#"Ladders",   "Beam_Remnants",  "Cross_Sections",  "Eikonals",  "Event_Generation",  "Main",  "Tools"
    if (text=="RECONNECTIONS"):
        return "Reconnections"
    if (text=="RemnantsMain"):
        return "Remnants"
        
    if (text=="CT14"):
        return "CT14Sherpa"
    if (text=="LHAPDF"):
        return "LHAPDFSherpa"
#Current_Library", "ME_Library", "Main", "PS_Library             
    return text
########################################################################    
def get_simple_so_version(text):
    return "0.0.0"
########################################################################
def get_full_so_version(text):
     return get_simple_so_version(text)
     #+".${SHERPA_VERSION_MAJOR}"
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
    newlist = [x.replace("if ENABLE_UFO","if (ENABLE_UFO)") for x in newlist]
    newlist = [x.replace("if PYTHIA_SUPPORT","if (PYTHIA_SUPPORT)") for x in newlist]
    newlist = [x.replace("if GZIP_SUPPORT","if (GZIP_SUPPORT)") for x in newlist]
    newlist = [x.replace("PYTHIAHEADERS =","set(PYTHIAHEADERS ") for x in newlist]
    newlist = [x.replace("GZIPSTREAMHEADERS =","set(GZIPSTREAMHEADERS ") for x in newlist]
    newlist = [x.replace("GZIPSTREAMSOURCES =","set(GZIPSTREAMSOURCES ") for x in newlist]
    newlist = [x.replace("PYTHIASOURCES =","set(PYTHIASOURCES ") for x in newlist]
    newlist = [x.replace("GZIPEXTRADIST =","set(GZIPEXTRADIST ") for x in newlist]
    newlist = [x.replace("$(GZIPSTREAMSOURCES)","${GZIPSTREAMSOURCES}") for x in newlist]
    newlist = [x.replace("$(PYTHIASOURCES)","${PYTHIASOURCES}") for x in newlist]
    newlist = [x.replace("else","else()") if re.match(r'^[:blank:]*else.*',x)  else x for x in newlist]
    newlist = [x.replace("endif","endif()") if re.match(r'^[:blank:]*endif.*',x)  else x for x in newlist]





    newlist = [x.replace("\\n"," ") for x in newlist]
    newlist = [x.replace("    "," ") for x in newlist]
    newlist = [x.replace("   "," ") for x in newlist]
    newlist = [x.replace("  "," ") for x in newlist]
    newlist = [x.replace("# ","#") for x in newlist]
    newlist = [x.replace("# ","#") for x in newlist]
    newlist = [x.strip() for x in newlist]
    newlist = filter(lambda st: st != '' , newlist)
    
    
    
    
    newlist = [  x.replace("&&"," AND ")  for x in newlist]
    newlist = [  x.replace("||"," OR ")  for x in newlist]

    newlist = [x.strip() for x in newlist]
    #newlist = filter(lambda st: st != '' , newlist)
    newlist = [x.replace("set(","  set(") if not re.match(r'^  .*',x)  else x for x in newlist]
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


    #filtered = ''.join(file_contents).split('\n')
    filtered = [""]
    for a in file_contents:
      if len(a)>2:
       print(a,"->",a[-2],"<-")
      if  len(a)>2 and a[-2]=='\\':
       if  len(a)>2 :
         filtered[-1]+=a[0:-2]
      else:
        filtered+=[a]  
    
    filtered= filter(lambda st: st != '\n' , filtered)
    filtered= filter(lambda st: st != '' , filtered)
    filtered= filter(lambda st: st != ' ' , filtered)
    filtered= filter(lambda st: not re.match(".*_la_LIBADD.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*_la_LIBADD.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*_CPPFLAGS.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*_CXXFLAGS.*",st), filtered)
    filtered= filter(lambda st: not re.match(".*_LIBADD.*",st), filtered)
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
    
    print("===========", LISTNAME)
    print(newlist) 
    incfiles = [ x  for x in newlist if re.match(r'#include.*',x)]
    incfiles =  [ x.replace("#include","") for x in incfiles]

    incfiles =  [ x.replace("\"","") for x in incfiles]
    incfiles =  [ x.replace(" ","") for x in incfiles]
    newlist.append(temp)
   # newlist = [ ("#"+x) if re.match(r'^.*_LIBADD.*',x) else x for x in newlist]
   # newlist = [ ("#"+x) if re.match(r'^.*_CPPFLAGS.*',x) else x for x in newlist]
    print("===========")
    print(newlist)
    newlist = filter(lambda st: st != '' , newlist)
  #  newlist = [x.replace(" \\ "," ") for x in newlist]

    newlist = [x.replace("if ENABLE_UFO","if (ENABLE_UFO)") for x in newlist]	
    newlist = [x.replace("if PYTHIA_SUPPORT","if (PYTHIA_SUPPORT)") for x in newlist]
    newlist = [x.replace("if GZIP_SUPPORT","if (GZIP_SUPPORT)") for x in newlist]
    newlist = [x.replace("NNPDF3archive =","set(NNPDF3archive  ") for x in newlist]
    newlist = [x.replace("SOBOL_EXT =","set(SOBOL_EXT  ") for x in newlist]
    newlist = [x.replace("PYTHIAHEADERS =","set(PYTHIAHEADERS ") for x in newlist]
    newlist = [x.replace("GZIPSTREAMHEADERS =","set(GZIPSTREAMHEADERS ") for x in newlist]
    newlist = [x.replace("GZIPSTREAMSOURCES =","set(GZIPSTREAMSOURCES ") for x in newlist]
    newlist = [x.replace("PYTHIASOURCES =","set(PYTHIASOURCES ") for x in newlist]
    newlist = [x.replace("GZIPEXTRADIST =","set(GZIPEXTRADIST ") for x in newlist]
    newlist = [x.replace("$(GZIPSTREAMSOURCES)","${GZIPSTREAMSOURCES}") for x in newlist]
    newlist = [x.replace("$(PYTHIASOURCES)","${PYTHIASOURCES}") for x in newlist]
    newlist = [x.replace("else","ELSE()") if re.match(r'^[:blank:]*else.*',x)  else x for x in newlist]
    newlist = [x.replace("endif","ENDIF()") if re.match(r'^[:blank:]*endif.*',x)  else x for x in newlist]
    newlist = [x.replace("include","#include") if re.match(r'^[:blank:]*include',x)  else x for x in newlist]
    newlist = [x.replace("SYSLIBS","#SYSLIBS") if re.match(r'^[:blank:]*SYSLIBS',x)  else x for x in newlist]
    newlist = [x.replace("dist-hook","#dist-hook") for x in newlist]
    
    
    
    newlist = [x.replace("AM_FFLAGS","#AM_FFLAGS")  for x in newlist]
    newlist = [x.replace("dist_ufo_PYTHON","#dist_ufo_PYTHON")  for x in newlist]
    newlist = [x.replace("ufodir =","#ufodir =")  for x in newlist]
    newlist = [x.replace("localinc","#localinc")  for x in newlist]
    newlist = [x.replace("pkglib_LTLIBRARIES","#pkglib_LTLIBRARIES")  for x in newlist]
    newlist = [x.replace("GITTAG","#GITTAG")  for x in newlist]
    newlist = [x.replace("EXTRA_DIST","#EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("rm -f","#rm -f")  for x in newlist]
    newlist = [x.replace("MAKE =","#MAKE =")  for x in newlist]
    newlist = [x.replace("nobase_","#nobase_")  for x in newlist]
    newlist = [x.replace("dist_bin_SCRIPTS","#dist_bin_SCRIPTS")  for x in newlist]
    newlist = [x.replace("dist_pkgdata","#dist_pkgdata")  for x in newlist]
    newlist = [x.replace("DIST_SUBDIRS","#DIST_SUBDIRS")  for x in newlist]
    newlist = [x.replace("SUBDIRS","#SUBDIRS")  for x in newlist]
    newlist = [x.replace("SYSLIBS","#SYSLIBS")  for x in newlist]
    newlist = [x.replace("include","#include")  for x in newlist]
    newlist = [x.replace("MD5_EXCLUDE","#MD5_EXCLUDE")  for x in newlist]
    newlist = [x.replace("MD5_EXCLUDE","#MD5_EXCLUDE")  for x in newlist]
    newlist = [x.replace("uninstall-hook","#uninstall-hook")  for x in newlist]
    newlist = [x.replace("install-data-hook","#install-data-hook")  for x in newlist]
    newlist = [x.replace("tar xjf","#tar xjf")  for x in newlist]
    newlist = [x.replace("-rm -rf","#-rm -rf")  for x in newlist]
    newlist = [x.replace("Sobol/%.gz","#Sobol/%.gz")  for x in newlist]
    newlist = [x.replace("mkdir -p Sobol","#mkdir -p Sobol")  for x in newlist]
    newlist = [x.replace("gzip < $< > $@","#gzip < $< > $@")  for x in newlist]


    newlist=[ x.replace(" \n","\n") for x in newlist]
    newlist=[ x.replace(" \n","\n") for x in newlist]


    print("===========2")
    print(newlist)
    newlistx = ''.join(newlist)
    newlistx=[ x.replace("\\\\\n","\n") for x in newlistx]
    newlistx=[ x.replace("\\\n","\n") for x in newlistx]
    print(newlistx) 
    newlist =''.join(newlistx).split('\n')

    print("===========3")
    print(newlist)    
    
    #newlist = [ ("#"+x) if re.match(r'^[:blank:]*[^#]*pkglib_LTLIBRARIES.*',x) else x for x in newlist]
    #newlist = [ ("\n#"+x) if re.match(r'^[:blank:]*[^#]*dist_bin_SCRIPTS.*',x) else x for x in newlist]
    #newlist = [ ("\n#"+x) if re.match(r'^[:blank:]*[^#]*SUBDIRS.*',x) else x for x in newlist]
    #newlist = [x.replace("GITTAG","\n#GITTAG") for x in newlist]
    #newlist = [x.replace("bin_PROGRAMS","\n#bin_PROGRAMS") for x in newlist]
    #newlist = [x.replace("nobase_dist_pkgdata","\n#nobase_dist_pkgdata") for x in newlist]
    #newlist = [x.replace("localinc","\n#localinc") for x in newlist]
    newlist = [x.replace("Git_Info.C","") for x in newlist]
    l=[]
    for x in newlist:
      if re.match(r'.*defined.*',x):
        l.append("#ORIGINAL "+x)
      l.append(x)

    newlist=l    
  #  newlist = [x.replace("_SOURCES =","set("+LISTNAME+"_SOURCES") for x in newlist]
    newlist = ["set("+LISTNAME+"_SOURCES" + x.split("=")[1] if re.match(r'.*_SOURCES =.*',x) else x  for x in newlist]
    newlist=[ x.replace(" \n","\n") for x in newlist]
    
    newlist = [x.replace("endif)","endif()")  for x in newlist]
    newlist = [x.replace("else)","else()")  for x in newlist]

    
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
    fin = [ "#"+x  if re.match(r'.*la_LIBADD.*',x)  else x for x in fin]
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
def create_library(ldirsI,lname,includes,installincludes,linklibs=[], cdff=[],pat=[ "\"*makefile*\"", "\"*\.c\"" ],subdir=""):
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
          
          a=a.replace("else)","else()")
          a=a.replace("endif)","endif()")
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
   suffixes=get_suffixes(lname)
   lname=subdir+lname
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
       f.write("target_link_libraries("+lname+suff+" PRIVATE "+ll+")\n")
     if if_install_library(lname):
       f.write("install(TARGETS "+lname+suff+" DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)\n")   
     if suff=="_static": 
       f.write("set_target_properties("+lname+"_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME "+output_name_static(lname)+")\n")
#       if (link_static(lname)!="no"): 
#         if if_install_library(lname):
#           f.write("install_symlink(lib"+output_name_static(lname)+".a "+" ${CMAKE_INSTALL_LIBDIR}/lib"+link_static(lname)+".a)\n")
     if suff=="":
       f.write("set_target_properties("+lname+"        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME "+output_name(lname)+" SOVERSION "+get_full_so_version(lname)+")\n")
#       if if_install_library(lname):
#         f.write("install_symlink(lib"+output_name(lname)+".so."+get_full_so_version(lname)+" "+"${CMAKE_INSTALL_LIBDIR}/lib"+output_name(lname)+".so."+get_simple_so_version(lname)+")\n")
     f.write("endif()\n")
   for inc in installincludes:
     if len(inc)!=0:
       f.write("install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/"+inc+" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SHERPA-MC/"+lname+"  COMPONENT devel  FILES_MATCHING PATTERN \"*.h\" PATTERN \"*.H\"  PATTERN \"*.deps*\" EXCLUDE   PATTERN \"*.lib*\" EXCLUDE")
       for pt in pat:
        f.write(" PATTERN "+ pt +" EXCLUDE ")
       f.write(")\n")
   f.close()
########################################################################
if __name__ == '__main__':
########################################################################
#   ldirs =  "Decays/Makefile.am Formation/Makefile.am Main/Makefile.am Makefile.am Tools/Makefile.am".split(" ")
#   lname="AHADIC++"
#   includes=" ".split(" ") 
#   installincludes=" ".split(" ")
#   create_library(ldirs,lname,includes,installincludes)

   os.chdir("AHADIC++")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "), 
             "Makefile.am".split(" ")
            ] 
   MClname= ["Decays", "Formation", "Main" , "Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Ahadic")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")

#   ldirs =  "Main/Makefile.am Makefile.am Perturbative/Makefile.am Tools/Makefile.am".split(" ")
#   lname="AMISIC++"
#   includes=" ".split(" ") 
#   installincludes=" ".split(" ") 
#   create_library(ldirs,lname,includes,installincludes)


   os.chdir("AMISIC++")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Main", "Perturbative",  "Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Amisic")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")
########################################################################
   ldirs =  "Amplitude/Makefile.am Amplitude/Zfunctions/Makefile.am DipoleSubtraction/Makefile.am Main/Makefile.am Makefile.am Phasespace/Makefile.am String/Makefile.am".split(" ")
   lname="AMEGIC++"
   includes=" ".split(" ") 
   installincludes=" ".split(" ")
   create_library(ldirs,lname,includes,installincludes)




########################################################################
   ldirs =  "Channels/Makefile.am Decays/Makefile.am Enhance/Makefile.am Main/Makefile.am Makefile.am Process/Makefile.am Scales/Makefile.am Selectors/Makefile.am".split(" ")
   lname="PHASIC++"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
  # create_library(ldirs,lname,includes,installincludes)

   os.chdir("PHASIC++")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Channels", "Decays", "Main", "Enhance", "Process", "Scales", "Selectors"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Phasic")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")


########################################################################
   #ldirs =  "Makefile.am Main/Makefile.am Spectra/Makefile.am".split(" ")
   #lname="BEAM"
   #includes=" ".split(" ") 
   #installincludes=" ".split(" ") 
   #create_library(ldirs,lname,includes,installincludes)
   
   os.chdir("BEAM")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Main", "Spectra"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Beam")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")
   
   
########################################################################
#   ldirs =  "Current_Library/Makefile.am ME_Library/Makefile.am Main/Makefile.am Makefile.am PS_Library/Makefile.am".split(" ")
#   lname="HADRONS++"
#   includes=" ".split(" ") 
#   installincludes=" ".split(" ") 
#   create_library(ldirs,lname,includes,installincludes)

   os.chdir("HADRONS++")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Current_Library", "ME_Library", "Main", "PS_Library"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Hadrons")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")


########################################################################

   ldirs =  "Calculators/Makefile.am Main/Makefile.am Makefile.am Showers/Makefile.am Tools/Makefile.am".split(" ")
   lname="CSSHOWER++"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)

   os.chdir("CSSHOWER++")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Main", "Calculators", "Showers", "Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"CS")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")   

   #ldirs =  "Gauge/Makefile.am Lorentz/Makefile.am Main/Makefile.am Makefile.am Shower/Makefile.am Tools/Makefile.am".split(" ")
   #lname="DIM"
   #includes=" ".split(" ") 
   #installincludes=" ".split(" ") 
   #create_library(ldirs,lname,includes,installincludes)
 
   os.chdir("DIM")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Gauge", "Lorentz", "Main", "Shower", "Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"DIM")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")

 
 
 
#   ldirs =  "Gauge/Makefile.am Lorentz/Makefile.am Main/Makefile.am Makefile.am Shower/Makefile.am Tools/Makefile.am".split(" ") 
   lname="DIRE"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
#   create_library(ldirs,lname,includes,installincludes)

   os.chdir("DIRE")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Gauge", "Lorentz", "Main", "Shower", "Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Dire")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")


   ldirs =  "Makefile.am".split(" ") 
   lname="EXTAMP"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)

   ldirs =  "Main/Makefile.am Makefile.am".split(" ")
   lname="RECONNECTIONS"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)

   ldirs =  "Main/Makefile.am Tools/Makefile.am".split(" ")
   lname="REMNANTS"
   includes=" ".split(" ") 
   installincludes="Main Tools".split(" ") 
   installincludes=" ".split(" ")
   create_library(ldirs,lname,includes,installincludes)


   os.chdir("REMNANTS")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Main","Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Remnants")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")  


   ldirs =  "Main/Makefile.am Makefile.am NLO/Makefile.am One2Three/Makefile.am One2Two/Makefile.am Two2Two/Makefile.am".split(" ")
   lname="EXTRA_XS"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)


#   ldirs =  "Calculators/Makefile.am Main/Makefile.am Makefile.am Showers/Makefile.am Tools/Makefile.am".split(" ") 
#   lname="MCATNLO"
#   includes=" ".split(" ") 
#   installincludes=" ".split(" ") 
#   create_library(ldirs,lname,includes,installincludes)

   os.chdir("MCATNLO")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Calculators", "Showers", "Main",  "Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"MCatNLO")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")




   #ldirs =  "Colors/Makefile.am Currents/Makefile.am Explicit/Makefile.am Loops/Makefile.am Main/Makefile.am Makefile.am SpinCorrelations/Makefile.am Vertices/Makefile.am".split(" ") 
  # lname="METOOLS"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
 #  create_library(ldirs,lname,includes,installincludes)
 
   os.chdir("METOOLS")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Colors",  "Currents",  "Explicit",  "Loops",  "Main",  "SpinCorrelations",  "Vertices"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"METools")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")  
 
 
      
  # ldirs =  "HEFT/Makefile.am Main/Makefile.am Makefile.am SM/Makefile.am SMEHC/Makefile.am TauPi/Makefile.am UFO/Makefile.am".split(" ") 
  # ldirs =  "Main/Makefile.am Makefile.am SM/Makefile.am TauPi/Makefile.am UFO/Makefile.am".split(" ") 
  # lname="MODEL"
  # includes=" ".split(" ") 
  # installincludes=" ".split(" ") 
  # create_library(ldirs,lname,includes,installincludes)


   os.chdir("MODEL")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["SM", "Main", "TauPi", "UFO", "HEFT", "SMGold", "SMEHC", "SMDM"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Sherpa")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")  

         ########################################################################
 #  ldirs =  "Amplitude/Makefile.am Main/Makefile.am Makefile.am Phasespace/Makefile.am".split(" ")
 #  lname="COMIX"
 #  includes=" ".split(" ") 
 #  installincludes=" ".split(" ")  
 #  create_library(ldirs,lname,includes,installincludes)
   
   os.chdir("COMIX")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Main", "Amplitude", "Phasespace"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Comix")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")   
   
########################################################################
   ldirs =  "Makefile.am Math/Makefile.am Org/Makefile.am Phys/Makefile.am YAML/Makefile.am YAML/yaml-cpp/Makefile.am YAML/yaml-cpp/node/Makefile.am YAML/yaml-cpp/node/detail/Makefile.am".split()
   lname="ATOOLS"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)
########################################################################
#   ldirs =  "MEs/Makefile.am Main/Makefile.am Makefile.am PhaseSpace/Makefile.am Tools/Makefile.am".split(" ") 
#   lname="PHOTONS++"
#   includes=" ".split(" ") 
#   installincludes=" ".split(" ")
#   create_library(ldirs,lname,includes,installincludes)


   os.chdir("PHOTONS++")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["MEs", "Main", "PhaseSpace", "Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Photons")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")  


   os.chdir("PDF")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "), 
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["CT14", "Electron", "GRV", "GRS", "SASG", "SAL", "Main", "NNPDF", "LHAPDF"]
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[])
   os.chdir("../")

   ldirs =  "Beam_Remnants.old/Makefile.am Beam_Remnants/Makefile.am Cross_Sections/Makefile.am Eikonals/Makefile.am Event_Generation/Makefile.am Main/Makefile.am Makefile.am Tools/Makefile.am".split()
   ldirs =  "Ladders/Makefile.am  Beam_Remnants/Makefile.am Cross_Sections/Makefile.am Eikonals/Makefile.am Event_Generation/Makefile.am Main/Makefile.am Makefile.am Tools/Makefile.am".split()
   lname="SHRiMPS"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
#   create_library(ldirs,lname,includes,installincludes)   


   os.chdir("SHRiMPS")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Ladders",   "Beam_Remnants",  "Cross_Sections",  "Eikonals",  "Event_Generation",  "Main",  "Tools"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Shrimps")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../") 

   ldirs =  "Initialization/Makefile.am LundTools/Makefile.am Main/Makefile.am Makefile.am PerturbativePhysics/Makefile.am Single_Events/Makefile.am SoftPhysics/Makefile.am Tools/Makefile.am".split()
   lname="SHERPA"
   includes=" ".split(" ") 
   installincludes=" ".split(" ") 
   create_library(ldirs,lname,includes,installincludes)   
