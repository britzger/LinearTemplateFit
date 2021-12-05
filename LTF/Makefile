###########################################################
#
#  Simple GNU makefile for LTF
#
###########################################################


# check if ROOT exists
WITH_ROOT := $(shell command -v root-config 2> /dev/null)
ifdef WITH_ROOT
WITH_ROOT=-D__WITH_ROOT__ $(shell root-config --libs) -I$(shell root-config --incdir)
endif


###########################################################
# all
all: slib bin


###########################################################
# libraries 
slib:
	$(CXX) -shared -o libLTF.so src/LTF.cxx -O0 -fPIC -g -Wall -Wno-unused-variable -std=gnu++17 -I.

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
bin: $(programs)

$(programs): %: %.cxx
	$(CXX) -o $@ $< -O0 -g -Wall -Wno-sign-compare  -Wl,-rpath,. -std=gnu++17 -I. -L. -lLTF -lstdc++ $(WITH_ROOT) 


###########################################################
# cleanup
clean:
	rm libLTF.so
	rm $(programs)

.PHONY: all root-slib clean


###########################################################
#bin: 
#	$(CXX) example_CMSinclusivejets_NN30_BRSSWpaper.cxx -o example_CMSinclusivejets_NN30_BRSSWpaper -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -I. -L. -lLTF  -lstdc++ $(WITH_ROOT) 
#	$(CXX) bb4l_mtop.cxx -o bb4l_mtop -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -lstdc++ 
#example:
#	$(CXX) example_LTF.cxx -o example_LTF -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -I$(shell root-config --incdir) -lstdc++ 
