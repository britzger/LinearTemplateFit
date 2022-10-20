############################################################################
# CMakeLists.txt file for building ROOT LTF package
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(LTF
  HEADERS
    LTF/LTF.h
    LTF/LTF_ROOTTools.h
    LTF/LTF_Tools.h
  SOURCES
    src/LTF.cxx
    src/LTF_ROOTTools.cxx
    src/LTF_Tools.cxx
  LINKDEF
    LTF_LinkDef.h
  DEPENDENCIES
  ROOT::MathCore ROOT::Gpad ROOT::Hist ROOT::Graf ROOT::Minuit ROOT::Matrix
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)
