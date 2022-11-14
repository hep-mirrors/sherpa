import sys
import sys
import os
import re
########################################################################
def get_suffixes(text):
#    return ["_static", ""]
    return [""]
def if_install_library(text):
    return True
def output_name_static(text):
    return output_name(text)
########################################################################    
#"Electron", "GRV", "GRS", "SASG", "SAL", "Main", "NNPDF", "LHAPDF"
def output_name(text):
    if (text=="ToolsYAML"):
        return "ToolsYaml"
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

    if (text=="AmegicString"):
        return "String"

    if (text=="AmegicPhasespace"):
        return "AmegicPSGen"

    if (text=="AmegicMain"):
        return "Amegic"

    if (text=="AmegicDipoleSubtraction"):
        return "DipoleSubtraction"
        
    if (text=="AmegicAmplitude"):
        return "Amplitude"

    if (text=="AmegicPhaseSpace"):
        return "AmegicPSGen"

    if (text=="ComixMain"):
        return "Comix"

#    if (text=="SherpaMain"):
#        return "ModelMain"

    if (text=="ExtraXSMain"):
       return "ExtraXS"

    if (text=="ExtraXSOne2Three"):
       return "ExtraXS1_3"

    if (text=="ExtraXSOne2Two"):
       return "ExtraXS1_2"

    if (text=="ExtraXSTwo2Two"):
       return "ExtraXS2_2"
#"NLO", "One2Three", "One2Two", "Two2Two"
    if (text=="SherpaSingle_Events"):
      return "SherpaSingleEvents"
    if (text=="SherpaLundTools"):
      return "LundTools"
#SM", "Main", "TauPi", "UFO", "HEFT", "SMGold", "SMEHC", "SMDM"
    if (text=="ModelTauPi"):
        return "SherpaTauPi"
    if (text=="ModelHEFT"):
        return "SherpaHEFT"
    if (text=="ModelSM"):
        return "SherpaSM"
    if (text=="ModelSMGold"):
        return "SherpaSMGold"
    if (text=="ModelSMEHC"):
        return "SherpaSMEHC"
    if (text=="ModelSMDM"):
        return "SherpaSMDM"
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
    if (text=="SherpaLH_OLE"):
        return "SherpaLHOLE"
    if (text=="SherpaHZTool"):
        return "SherpaHZToolAnalysis"

    if (text=="SherpaNNLO"):
        return "NNLOqT"

    if (text=="SherpaAnalysisAnalyses"):
        return "SherpaAnalyses"

    if (text=="SherpaAnalysisMain"):
        return "SherpaAnalysis"

    if (text=="SherpaAnalysisTriggers"):
        return "SherpaAnalysisTrigger"

    if (text=="SherpaAnalysisObservables"):
        return "SherpaObservables"

    if (text=="SherpaRivet"):
        return "SherpaRivetAnalysis"

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
#  Automatically or semiautomaticaly generated, do not edit.
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
  return ["",""]
########################################################################
def transform_pilots(argv):
    path_to_file=argv
    file_contents=''
    with open(path_to_file) as f:
      file_contents = f.readlines()
    newlist =file_contents
    newlist = [x.replace("if ENABLE_DIHIGGS","if (SHERPA_ENABLE_DIHIGGS)") for x in newlist]
    newlist = [x.replace("if ENABLE_UFO","if (SHERPA_ENABLE_UFO)") for x in newlist]
    newlist = [x.replace("if PYTHIA_SUPPORT","if (SHERPA_ENABLE_PYTHIA)") for x in newlist]
    newlist = [x.replace("if GZIP_SUPPORT","if (SHERPA_ENABLE_GZIP)") for x in newlist]
    newlist = [x.replace("PYTHIAHEADERS =","set(PYTHIAHEADERS ") for x in newlist]
    newlist = [x.replace("GZIPSTREAMHEADERS =","set(GZIPSTREAMHEADERS ") for x in newlist]
    newlist = [x.replace("GZIPSTREAMSOURCES =","set(GZIPSTREAMSOURCES ") for x in newlist]
    newlist = [x.replace("PYTHIASOURCES =","set(PYTHIASOURCES ") for x in newlist]
    newlist = [x.replace("GZIPEXTRADIST =","set(GZIPEXTRADIST ") for x in newlist]
    newlist = [x.replace("$(GZIPSTREAMSOURCES)","${GZIPSTREAMSOURCES}") for x in newlist]
    newlist = [x.replace("$(GZIPSTREAMHEADERS)","${GZIPSTREAMHEADERS}") for x in newlist]
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
    
    return newlist

def transform_imake_source(argv, dbg):
    path_to_file=argv
    with open(path_to_file) as f:
      file_contents = f.readlines()

    #filtered = ''.join(file_contents).split('\n')
    filtered = [""]
    for a in file_contents:
#      if len(a)>2:
#       print(a,"->",a[-2],"<-")
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
    incfiles = [ x  for x in newlist if re.match(r'#include.*',x)]
    incfiles =  [ x.replace("#include","") for x in incfiles]

    incfiles =  [ x.replace("\"","") for x in incfiles]
    incfiles =  [ x.replace(" ","") for x in incfiles]
    newlist.append(temp)
    newlist = filter(lambda st: st != '' , newlist)

    newlist = [x.replace("if USING__Analysis","if (SHERPA_ENABLE_ANALYSIS)") for x in newlist]	
    newlist = [x.replace("if USING__EWSud","if (SHERPA_ENABLE_EWSUD)") for x in newlist]	
    newlist = [x.replace("if USING__LHOLE","if (SHERPA_ENABLE_LHOLE)") for x in newlist]	
    newlist = [x.replace("if ENABLE_UFO","if (SHERAP_ENABLE_UFO)") for x in newlist]	
    newlist = [x.replace("if PYTHIA8_SUPPORT","if (SHERPA_ENABLE_PYTHIA8)") for x in newlist]
    newlist = [x.replace("if HZTOOL_SUPPORT","if (SHERPA_ENABLE_HZTOOL)") for x in newlist]
    newlist = [x.replace("if CERNLIB_SUPPORT","if (SHERPA_ENABLE_CERNLIB)") for x in newlist]
    newlist = [x.replace("if GOSAM_SUPPORT","if (SHERPA_ENABLE_GOSAM)") for x in newlist]
    newlist = [x.replace("if PYTHIA_SUPPORT","if (SHERPA_ENABLE_PYTHIA)") for x in newlist]
    newlist = [x.replace("if GZIP_SUPPORT","if (SHERPA_ENABLE_GZIP)") for x in newlist]
    newlist = [x.replace("if BLACKHAT_SUPPORT","if (SHERPA_ENABLE_BLACKHAT)") for x in newlist]
    newlist = [x.replace("NNPDF3archive =","set(NNPDF3archive  ") for x in newlist]
    newlist = [x.replace("EWSUD_ADDS =","set(EWSUD_ADDS  ") for x in newlist]
    newlist = [x.replace("LHOLE_ADDS =","set(LHOLE_ADDS  ") for x in newlist]
    newlist = [x.replace("BLACKHAT_ADDS =","set(BLACKHAT_ADDS  ") for x in newlist]
    newlist = [x.replace("GOSAM_ADDS =","set(GOSAM_ADDS  ") for x in newlist]
    newlist = [x.replace("HZTOOL_ADDS =","set(HZTOOL_ADDS  ") for x in newlist]
    newlist = [x.replace("SOBOL_EXT =","set(SOBOL_EXT  ") for x in newlist]
    newlist = [x.replace("PYTHIAHEADERS =","set(PYTHIAHEADERS ") for x in newlist]
    newlist = [x.replace("GZIPSTREAMHEADERS =","set(GZIPSTREAMHEADERS ") for x in newlist]
    newlist = [x.replace("GZIPSTREAMSOURCES =","set(GZIPSTREAMSOURCES ") for x in newlist]
    newlist = [x.replace("PYTHIASOURCES =","set(PYTHIASOURCES ") for x in newlist]
    newlist = [x.replace("GZIPEXTRADIST =","set(GZIPEXTRADIST ") for x in newlist]
    newlist = [x.replace("ANAANA_ADDS =","set(ANAANA_ADDS ") for x in newlist]
    newlist = [x.replace("ANATRIGGER_ADDS =","set(ANATRIGGER_ADDS ") for x in newlist]
    newlist = [x.replace("ANAMAIN_ADDS =","set(ANAMAIN_ADDS ") for x in newlist]
    newlist = [x.replace("ANAOBS_ADDS =","set(ANAOBS_ADDS ") for x in newlist]
    newlist = [x.replace("ANATRIGGER_ADDS =","set(ANATRIGGER_ADDS ") for x in newlist]
    newlist = [x.replace("ANATOOLS_ADDS =","set(ANATOOL_ADDS ") for x in newlist]
    newlist = [x.replace("DIHIGGS_ADDS =","set(DIHIGGS_ADDS ") for x in newlist]
    newlist = [x.replace("$(GOSAM_SOURCES)","${GoSam_SOURCES}") for x in newlist]
    newlist = [x.replace("$(BLACKHAT_SOURCES)","${BLACKHAT_SOURCES}") for x in newlist]
    newlist = [x.replace("$(EWSUD_SOURCES)","${EWSud_SOURCES}") for x in newlist]
    newlist = [x.replace("$(LHOLE_SOURCES)","${LH_OLE_SOURCES}") for x in newlist]
    newlist = [x.replace("$(HZTOOL_SOURCES)","${HZTOOL_SOURCES}") for x in newlist]
    newlist = [x.replace("$(EWSUD_ADDS)","${EWSUD_ADDS}") for x in newlist]
    newlist = [x.replace("$(GZIPSTREAMSOURCES)","${GZIPSTREAMSOURCES}") for x in newlist]
    newlist = [x.replace("$(GZIPSTREAMHEADERS)","${GZIPSTREAMHEADERS}") for x in newlist]    
    newlist = [x.replace("$(PYTHIASOURCES)","${PYTHIASOURCES}") for x in newlist]
    newlist = [x.replace("$(DIHIGGS_SOURCES)","${DIHIGGS_SOURCES}") for x in newlist]
    newlist = [x.replace("$(ANAMAIN_SOURCES)","${ANAMAIN_SOURCES}") for x in newlist]
    newlist = [x.replace("$(ANAMAIN_ADDS)","${ANAMAIN_ADDS}") for x in newlist]

    newlist = [x.replace("$(ANAOBS_SOURCES)","${ANAOBS_SOURCES}") for x in newlist]
    newlist = [x.replace("$(ANAOBS_ADDS)","${ANAOBS_ADDS}") for x in newlist]

    newlist = [x.replace("$(ANAANA_SOURCES)","${ANAANA_SOURCES}") for x in newlist]
    newlist = [x.replace("$(ANAANA_ADDS)","${ANAANA_ADDS}") for x in newlist]
    newlist = [x.replace("$(DIHIGGS_ADDS)","${DIHIGGS_ADDS}") for x in newlist]


    newlist = [x.replace("$(ANATOOLS_SOURCES)","${ANATOOLS_SOURCES}") for x in newlist]
    newlist = [x.replace("$(ANATRIGGER_SOURCES)","${ANATRIGGER_SOURCES}") for x in newlist]
    newlist = [x.replace("$(ANATOOLS_ADDS)","${ANATOOLS_ADDS}") for x in newlist]
    newlist = [x.replace("$(ANATRIGGER_ADDS)","${ANATRIGGER_ADDS}") for x in newlist]
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
    newlist = [x.replace("GITTAG =","set(GITTAG ")  for x in newlist]
    newlist = [x.replace("ANAMAIN_EXTRA_DIST","#ANAMAIN_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("EWSUD_EXTRA_DIST","#EWSUD_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("DIHIGGS_EXTRA_DIST","#DIHIGGS_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("HIGGS_EXTRA_DIST","#HIGGS_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("LHOLE_EXTRA_DIST","#LHOLE_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("BLACKHAT_EXTRA_DIST","#BLACKHAT_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("GOSAM_EXTRA_DIST","#GOSAM_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("HZTOOL_EXTRA_DIST","#HZTOOL_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("NNLOqT_EXTRA_DIST","#NNLOqT_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("ANATOOLS_EXTRA_DIST","#ANATOOLS_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("ANAOBS_EXTRA_DIST","#ANAOBS_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("ANAANA_EXTRA_DIST","#ANAANA_EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("ANATRIGGER_EXTRA_DIST","#ANATRIGGER_EXTRA_DIST")  for x in newlist]

    newlist = [x.replace("EXTRA_DIST","#EXTRA_DIST")  for x in newlist]
    newlist = [x.replace("rm -f","#rm -f")  for x in newlist]
    newlist = [x.replace("MAKE =","#MAKE =")  for x in newlist]
    newlist = [x.replace("AUTOMAKE_OPTIONS","#AUTOMAKE_OPTIONS")  for x in newlist]
   #newlist = [x.replace("nobase_","#nobase_")  for x in newlist]
    newlist = [x.replace("nobase_dist_pkgdata","#nobase_dist_pkgdata")  for x in newlist]
    newlist = [x.replace("dist_bin_SCRIPTS","#dist_bin_SCRIPTS")  for x in newlist]
    newlist = [x.replace("dist_pkgdata","#dist_pkgdata")  for x in newlist]
    newlist = [x.replace("DIST_SUBDIRS","#DIST_SUBDIRS")  for x in newlist]
    newlist = [x.replace("SUBDIRS","#SUBDIRS")  for x in newlist]
    newlist = [x.replace("SYSLIBS","#SYSLIBS")  for x in newlist]
    newlist = [x.replace("include","#include")  for x in newlist]
    newlist = [x.replace("MD5_EXCLUDE =","set(MD5_EXCLUDE ")  for x in newlist]
    newlist = [x.replace("uninstall-hook","#uninstall-hook")  for x in newlist]
    newlist = [x.replace("install-data-hook","#install-data-hook")  for x in newlist]
    newlist = [x.replace("tar xjf","#tar xjf")  for x in newlist]
    newlist = [x.replace("-rm -rf","#-rm -rf")  for x in newlist]
    newlist = [x.replace("Sobol/%.gz","#Sobol/%.gz")  for x in newlist]
    newlist = [x.replace("mkdir -p Sobol","#mkdir -p Sobol")  for x in newlist]
    newlist = [x.replace("gzip < $< > $@","#gzip < $< > $@")  for x in newlist]
    newlist = [x.replace("dist_dihiggs","#dist_dihiggs")  for x in newlist]
    newlist = [x.replace("dihiggs","#dihiggs")  for x in newlist]


    newlist=[ x.replace(" \n","\n") for x in newlist]
    newlist=[ x.replace(" \n","\n") for x in newlist]

    newlistx = ''.join(newlist)
    newlistx=[ x.replace("\\\\\n","\n") for x in newlistx]
    newlistx=[ x.replace("\\\n","\n") for x in newlistx]
  #  print(newlistx) 
    newlist =''.join(newlistx).split('\n')
    newlist = [x.replace("Git_Info.C","") for x in newlist]
    l=[]
    for x in newlist:
      if re.match(r'.*defined.*',x):
        l.append("#ORIGINAL "+x)
      l.append(x)

    newlist=l    
  #  newlist = [x.replace("_SOURCES =","set("+LISTNAME+"_SOURCES") for x in newlist]
    newlist = ["set("+LISTNAME+"_SOURCES" + x.split("=")[1] if re.match(r'.*_SOURCES =.*',x) else x  for x in newlist]
    newlist = ["set("+LISTNAME+"_HEADERS" + x.split("=")[1] if re.match(r'.*_HEADERS =.*',x) else x  for x in newlist]
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
   # print(fin,FSRC)
    
    diractual="../"+dirf
    diractual=diractual.replace("//","/")
    if FSRC==1:
      fin.append("creategitinfo("+LISTNAME+" "+diractual+")")
      fin.append("list(TRANSFORM "+LISTNAME+"_SOURCES PREPEND \"${CMAKE_CURRENT_SOURCE_DIR}/"+diractual+"\")")
      if ( LISTNAME != "Reconnections" and LISTNAME != "ToolsYaml" and LISTNAME != "ShrimpsLadders"): 
        fin.append("list(APPEND "+LISTNAME+"_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)")
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
   headerx=""
   try:
    ff = open(lname+"/CMakeLists.txt", "r")
    headerx=ff.readline()
    ff.close()
    #print(lname,"->",headerx,"<-")
   except:
    pass
   if (headerx=="#Manually edited\n"):
     return
   
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
   suffixes=get_suffixes(lname)
   lname=subdir+lname
   lname=output_name_static(lname)
   #f.write("set("+lname+"_esources )\n")
   for suff in suffixes:
     if suff=="": 
       #f.write("if (SHERPA_BUILD_SHARED)\n")
       f.write("add_library("+lname+" SHARED \n")
       #f.write("add_library("+lname+" SHARED ${"+lname+"_esources}\n")
       for z in t: f.write("                             "+z + " \n")
       f.write(")\n")
     if suff=="_static": 
       f.write("if (SHERPA_BUILD_STATIC)\n")
       f.write("add_library("+lname+"_static STATIC ${"+lname+"_esources}\n")
       for z in t: f.write("                             "+z + " \n")
       f.write(")\n")
     #f.write("target_include_directories("+lname+suff+" PRIVATE ${PROJECT_SOURCE_DIR}/include)\n")
     f.write("target_include_directories("+lname+suff+" PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})\n")
     for inc in includes:
       if os.path.exists(lname+"/"+inc):
         u="/"+inc
         u=u.replace("//","/")
         f.write("target_include_directories("+lname+suff+" PRIVATE \"${CMAKE_CURRENT_SOURCE_DIR}"+u+"\")\n")
#     f.write("target_include_directories("+lname+suff+" PRIVATE ${FREETYPE_INCLUDE_DIRS})\n")
     for ll in linklibs:
       f.write("target_link_libraries("+lname+suff+" PRIVATE "+ll+")\n")
     f.write("sherpa_mpi_link_libraries("+lname+suff+")\n")
     if if_install_library(lname):
       f.write("install(TARGETS "+lname+suff+" DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)\n")   
     if suff=="_static": 
       f.write("set_target_properties("+lname+"_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME "+output_name_static(lname)+")\n")
       f.write("set_target_properties("+lname+"_static PROPERTIES DEFINE_SYMBOL \"\")\n")
     if suff=="":
       f.write("set_target_properties("+lname+"        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME "+output_name(lname)+" SOVERSION "+get_full_so_version(lname)+")\n")
       f.write("set_target_properties("+lname+"        PROPERTIES DEFINE_SYMBOL \"\")\n")
     #f.write("endif()\n")
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

   os.chdir("AMEGIC++")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Amplitude",  "DipoleSubtraction", "Main", "Phasespace", "String"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Amegic")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")
   os.chdir("AMEGIC++/Amplitude")
   includes=" ".split(" ") 
   installincludes=" ".split(" ")
   create_library(["Makefile.am"],"Zfunctions",includes,installincludes,[],[],[],"")
   os.chdir("../")
   f = open("CMakeLists.txt", "a") 
   f.write("add_subdirectory(Amplitude/Zfunctions)\n")
   f.close()
   os.chdir("../")
   
########################################################################
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

   os.chdir("EXTRA_XS")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Main", "NLO", "One2Three", "One2Two", "Two2Two", "Special"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"ExtraXS")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")  

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
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Model")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   os.chdir("../")  

         ########################################################################

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
   os.chdir("ATOOLS")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am yaml-cpp/Makefile.am yaml-cpp/node/Makefile.am yaml-cpp/node/detail/Makefile.am".split(" ")
            ] 
   MClname= ["Math", "Org", "Phys", "YAML"]
   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Tools")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   ff=open("Org/CMakeLists.txt","a")
   ff.write("""
#if (SHERPA_ENABLE_GZIP)
target_link_libraries(ToolsOrg PRIVATE LibZip::LibZip)
#endif()
""")
   ff.close()   
   ff=open("Phys/CMakeLists.txt","a")
   ff.write("""
if (SHERPA_ENABLE_LHAPDF)
target_link_libraries(ToolsPhys PRIVATE ${LHAPDF_LIBRARIES})
target_include_directories(ToolsPhys PRIVATE ${LHAPDF_INCLUDE_DIRS})
endif()
""")
   ff.close() 
   os.chdir("../")    
   
########################################################################
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
   f=open("LHAPDF/CMakeFiles.txt","a")
   f.write("""
target_link_libraries(LHAPDFSherpa PRIVATE ${LHAPDF_LIBRARIES})
target_include_directories(LHAPDFSherpa PRIVATE ${LHAPDF_INCLUDE_DIRS})
""")
   f.close()
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


   os.chdir("SHERPA")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Main", "Initialization", "LundTools",  "PerturbativePhysics", "Single_Events", "SoftPhysics", "Tools"]

   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Sherpa")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   
   ff=open("Tools/CMakeLists.txt","a")
   ff.write("""
if (SHERPA_ENABLE_HEPMC2)
target_link_libraries(SherpaTools PRIVATE ${HEPMC2_LIBRARIES})
target_include_directories(SherpaTools PRIVATE ${HEPMC2_INCLUDE_DIRS})
endif()
""")
   ff.close() 
   os.chdir("../")

############################

   os.chdir("AddOns")
   MCldirs =["Makefile.am".split(" "),
#             "Makefile.am".split(" "),
#             "Makefile.am".split(" "),
#             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["EWSud","Higgs","LH_OLE","Weights"]

#   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Sherpa")
#     f.write("add_subdirectory("+lname+")\n")
#   f.close()   

   os.chdir("../")



   os.chdir("AddOns")
   MCldirs =["Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Rivet","HZTool","Pythia","NNLO", "OpenLoops", "HepMC", "BlackHat", "GoSam", "Root", "DiHiggsNLO", "PGS"]

#   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"Sherpa")
#     f.write("add_subdirectory("+lname+")\n")
#   f.close()   
   ff=open("Pythia/CMakeLists.txt","a")
   ff.write("""
target_link_libraries(SherpaPythia PRIVATE ${PYTHIA8_LIBRARIES})
target_include_directories(SherpaPythia PRIVATE ${PYTHIA8_INCLUDE_DIRS})
""")
   ff.close() 

   os.chdir("../")


   os.chdir("AddOns/Analysis")
   MCldirs =["Makefile.am".split(" "),
#             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" "),
             "Makefile.am".split(" ")
            ] 
   MClname= ["Tools", "Triggers", "Observables", "Analyses", "Main", "Scripts"]

   f = open("CMakeLists.txt", "w")  
   for x in range(0,len(MClname) ):
     lname=MClname[x]
     ldirs=MCldirs[x]
     includes="/".split(" ")+[lname]
     installincludes=[]
     create_library(ldirs,lname,includes,installincludes,[],[],[],"SherpaAnalysis")
     f.write("add_subdirectory("+lname+")\n")
   f.close()   

   os.chdir("../../")


# 'AddOns/Rivet/Git_Info.C', 
# 'AddOns/HZTool/Git_Info.C', 
# 'AddOns/Pythia/Git_Info.C', 
# 'AddOns/NNLO/Git_Info.C', 
# 'AddOns/OpenLoops/Git_Info.C', 
# 'AddOns/HepMC/Git_Info.C', 
# 'AddOns/BlackHat/Git_Info.C', 
# 'AddOns/GoSam/Git_Info.C', 




