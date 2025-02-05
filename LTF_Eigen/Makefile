###########################################################
#
#  Simple GNU makefile for LTF
#
###########################################################

# Source directories
EXAMPLE_DIR=examples
LTF_DIR=LTF
BUILD_DIR=build

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

SRC=$(LTF_DIR)/src/LTF.cxx
SRC+=$(LTF_DIR)/src/LTF_Tools.cxx
ifdef WITH_ROOT
SRC+=$(LTF_DIR)/src/LTF_ROOTTools.cxx
endif

ifeq (,$(wildcard $(EIGEN_INCLUDE)/Eigen/Eigen))
$(error Eigen not found! Please provide EIGEN__HOME, provide directory 'Eigen' or change Makefile accordingly)
endif

###########################################################
# all
all:
	$(MAKE) slib
	$(MAKE) bin
	$(MAKE) include


###########################################################
# libraries 
slib:
	@mkdir -p $(@D)/$(BUILD_DIR)/lib 
	$(CXX) -shared -o $(BUILD_DIR)/lib/libLTF.so $(SRC) -O0 -fPIC -g -Wall -Wno-unused-variable -std=gnu++17 -I$(LTF_DIR)/include -I$(EIGEN_INCLUDE) $(ROOTINCLUDES) $(ROOTLIBS) 

root-slib:
	export ROOT_INCLUDE_PATH=$(LTF_DIR)/include:$(EIGEN_INCLUDE) && root -b -q $(LTF_DIR)/src/LTF.cxx+


###########################################################
# programs
#programs =  $(EXAMPLE_DIR)/example_CMSinclusivejets_NN30_BRSSWpaper 
#programs += $(EXAMPLE_DIR)/example_CMSinclusivejets_MSTW_CMSpaper
#programs += $(EXAMPLE_DIR)/example1_LTF_gaus
#programs += $(EXAMPLE_DIR)/example1_LTF_gaus_NoROOT
#programs += $(EXAMPLE_DIR)/example2_LTF_gaus2D
#programs += $(EXAMPLE_DIR)/example2_LTF_gaus2D_NoROOT
#programs += $(EXAMPLE_DIR)/example3_LTF_gaus_sigma
programs += $(EXAMPLE_DIR)/example_ATLAS_topmass
bin: $(programs)

$(programs): %: %.cxx
	@mkdir -p $(BUILD_DIR)/bin 
	$(CXX) $< -o $(BUILD_DIR)/bin/$(notdir $@) -O0 -g -Wall -Wno-sign-compare  -Wno-unused-variable -Wl,-rpath,$(BUILD_DIR)/lib -std=gnu++17 -I$(LTF_DIR)/include -I$(EIGEN_INCLUDE) $(ROOTINCLUDES) -L$(BUILD_DIR)/lib -lLTF -lstdc++ $(ROOTLIBS)


include:
	@cp -rf $(LTF_DIR)/include $(BUILD_DIR)/.

###########################################################
# cleanup
clean:
	rm $(BUILD_DIR)/lib/libLTF.so
	rm -rf $(BUILD_DIR)/lib/libLTF.so.dSYM
	rm -rf $(BUILD_DIR)/bin/*.dSYM
	rm $(BUILD_DIR)/bin/*
	rm -rf $(BUILD_DIR)/include
	rm -rf $(BUILD_DIR)

.PHONY: all root-slib clean


###########################################################
#bin: 
#	$(CXX) example_CMSinclusivejets_NN30_BRSSWpaper.cxx -o example_CMSinclusivejets_NN30_BRSSWpaper -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -I. -L. -lLTF  -lstdc++ $(WITH_ROOT) 
#	$(CXX) bb4l_mtop.cxx -o bb4l_mtop -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -lstdc++ 
#example:
#	$(CXX) example_LTF.cxx -o example_LTF -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -I$(shell root-config --incdir) -lstdc++ 
