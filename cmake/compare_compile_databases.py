
#This script is designed to analyse the compilation databases created with cmake and autotools build systems.
import json,re,sys,os
def get_compilation_DB(fname):
 f = open(fname)
 data = json.load(f)
 L = {}
 for i in data:
  ar = i['arguments']
  fil = i['file']
  fil = fil.replace("../../../../","")
  fil = fil.replace("../../../","")
  fil = fil.replace("../../","")
  fil = fil.replace("../","")
  current_directory = os.getcwd()
  comp = ar[0]
  comparg = ar[1:]
  comparg = [ x.replace(current_directory+"build","") for x in comparg]
  comparg = [ x.replace(current_directory+"","") for x in comparg]
  comparg = [ x for x in comparg if not re.match(r'^-I/usr/include$',x)]
  comparg = [ x for x in comparg if x.startswith('-')]
  comparg = list( dict.fromkeys(comparg) )
  comparg.sort()
  if  not   re.match(r'.*Sherpa_wrap.cxx$',fil):
   if  not   re.match(r'.*pythia-6.4.18.f$',fil):
     L[fil] = comparg
 f.close()
 return L

def get_list_difference(li1, li2):
  return list(set(li1) - set(li2)) + list(set(li2) - set(li1))

def get_list(dict):
  list = []
  for key in dict.keys():
    list.append(key)
  return list

cmakeDB = get_compilation_DB("CM.json")
imakeDB = get_compilation_DB("AT.json")

cmakeList = get_list(cmakeDB)
cmakeList.sort()

imakeList = get_list(imakeDB)
imakeList.sort()


differenceList = get_list_difference(imakeList,cmakeList)

#Bear does not catch some complex cases in imake or cmake builds just more tests
differenceList = [ x for x in differenceList if not re.match(r'.*kernbit/test.*',x)]  
differenceList = [ x for x in differenceList if not re.match(r'.*kerngen/test.*',x)]  

filteredCmakeList = cmakeList
filteredCmakeList = [ x for x in filteredCmakeList if not re.match(r'.*kernbit/test.*',x)]  

print(differenceList)

#Everything
packages =["allfiles"]

comdefines ={}

for pack in packages:
 allDefinescmake=set()

 allDefinesimake=set()
 allWDefinescmake=set()

 allWDefinesimake=set()
 DiffDefines=[]
 for a in filteredCmakeList:
  if ( pack == "allfiles" ):
    compileOptions = get_list_difference(cmakeDB[a],imakeDB[a])
    Includes = [ x for x in compileOptions if re.match(r'^-I.*',x)]  
    compileOptions = [ x for x in compileOptions if not re.match(r'^-I.*',x)]  
    Definescmake = [ x for x in cmakeDB[a] if re.match(r'^-D.*',x)]  
    Definesimake = [ x for x in imakeDB[a] if re.match(r'^-D.*',x)]  
    comdefines[pack] = Definesimake
    Definescmake.sort()
    Definesimake.sort()
    DiffDefines=get_list_difference(Definesimake,Definescmake)
    #Filter some known differences
    DiffDefines=[ x for x in DiffDefines if not re.match(r'^-Wl,--no-as-needed.*',x)]  

    allDefinesimake.update(Definesimake)
    allDefinescmake.update(Definescmake)

    WDefinescmake = [ x for x in cmakeDB[a] if re.match(r'^-W.*',x)]  
    WDefinesimake = [ x for x in imakeDB[a] if re.match(r'^-W.*',x)]  
    allWDefinesimake.update(WDefinesimake)
    allWDefinescmake.update(WDefinescmake)

    #Filter non-flags and known differences
#    compileOptions = [ x for x in compileOptions if not re.match(r'^-D.*',x)]  
    compileOptions = [ x for x in compileOptions if not re.match(r'^-Wno-.*',x)]  
#    compileOptions = [ x for x in compileOptions if not re.match(r'^-no-pie.*',x)]  
    compileOptions = [ x for x in compileOptions if not re.match(r'^-pipe.*',x)]  
    if (len(compileOptions)>0 or len(DiffDefines)>0):
     print("Note file ",a)
     print("Defines: diff, cmake, imake")
     print(DiffDefines)
     print(Definescmake)
     print(Definesimake)
     print("Compile options: diff")
     print(compileOptions)
 allDiffDefines=get_list_difference(list(allDefinescmake),list(allDefinesimake))
 allDiffDefines=[ x for x in allDiffDefines if not re.match(r'^-Wl,--no-as-needed.*',x)]

 print(pack)
 if len(allDiffDefines)>0:
  print(allDiffDefines)
  cm=list(allDefinescmake)
  cm.sort()
  print(cm)
  im=list(allDefinesimake)
  im.sort()
  print(im)


 Wcm=list(allWDefinescmake)
 Wcm.sort()
 #print(Wcm)
 Wim=list(allWDefinesimake)
 Wim.sort()
 #print(Wim)

alld = set()
comm = set()
for x in comdefines.keys(): alld.update(comdefines[x])

comm = alld 
print("All defines")
print (comm)
for x in comdefines.keys():
    comm= comm.intersection(comdefines[x])
print("Common defines")
print (comm)

