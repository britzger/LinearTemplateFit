###########################################################
#
#  Simple GNU makefile for LTF
#
###########################################################


# check if ROOT exists
WITH_ROOT := $(shell command -v root-config 2> /dev/null)
ROOTLIBS=
ROOTINCLUDES=
ifdef WITH_ROOT
WITH_ROOT=-D__WITH_ROOT__
ROOTLIBS=$(shell root-config --libs)
ROOTINCLUDES=-I$(shell root-config --incdir) $(WITH_ROOT)
endif

# take eigen from local copy in 'Eigen_copy' or from
# variable "EIGEN__HOME"
EIGEN_INCLUDE = ./Eigen_copy
ifdef EIGEN__HOME
EIGEN_INCLUDE = $(EIGEN__HOME)/include/eigen3
endif

SRC=src/LTF.cxx
SRC+=src/LTF_Tools.cxx
ifdef WITH_ROOT
SRC+=src/LTF_ROOTTools.cxx
endif

ifeq (,$(wildcard $(EIGEN_INCLUDE)/Eigen/Eigen))
$(error Eigen not found! Please provide EIGEN__HOME, provide directory 'Eigen' or change Makefile accordingly)
endif

###########################################################
# all
all:
	$(MAKE) slib
	$(MAKE) bin


###########################################################
# libraries 
slib:
	@mkdir -p $(@D)/lib 
	$(CXX) -shared -o lib/libLTF.so $(SRC) -O0 -fPIC -g -Wall -Wno-unused-variable -std=gnu++17 -I. -I$(EIGEN_INCLUDE) $(ROOTINCLUDES) $(ROOTLIBS) 

root-slib:
	root -b -q src/LTF.cxx+


###########################################################
# programs
programs =  example_CMSinclusivejets_NN30_BRSSWpaper 
programs += example_CMSinclusivejets_MSTW_CMSpaper
programs += example1_LTF_gaus
programs += example1_LTF_gaus_NoROOT
programs += example2_LTF_gaus2D
programs += example2_LTF_gaus2D_NoROOT
programs += example3_LTF_gaus_sigma
programs += mjet_mtop
bin: $(programs)

$(programs): %: %.cxx
	@mkdir -p bin 
	$(CXX) $< -o bin/$(notdir $@) -O0 -g -Wall -Wno-sign-compare  -Wno-unused-variable -Wl,-rpath,$(PWD)/lib -std=gnu++17 -I. -I$(EIGEN_INCLUDE) $(ROOTINCLUDES) -Llib/ -lLTF -lstdc++ $(ROOTLIBS)


###########################################################
# cleanup
clean:
	rm lib/libLTF.so
	rm -rf lib/libLTF.so.dSYM
	rm -rf bin/*.dSYM
	rm bin/*


.PHONY: all root-slib clean


###########################################################
#bin: 
#	$(CXX) example_CMSinclusivejets_NN30_BRSSWpaper.cxx -o example_CMSinclusivejets_NN30_BRSSWpaper -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -I. -L. -lLTF  -lstdc++ $(WITH_ROOT) 
#	$(CXX) bb4l_mtop.cxx -o bb4l_mtop -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -lstdc++ 
#example:
#	$(CXX) example_LTF.cxx -o example_LTF -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -I$(shell root-config --incdir) -lstdc++ 
