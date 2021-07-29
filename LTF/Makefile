

# check if ROOT exists
WITH_ROOT := $(shell command -v root-config 2> /dev/null)
ifdef WITH_ROOT
WITH_ROOT=-D__WITH_ROOT__ $(shell root-config --libs) -I$(shell root-config --incdir)
endif

all: slib bin

example:
	$(CXX) example_LTF.cxx -o example_LTF -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -I$(shell root-config --incdir) -lstdc++ 

slib:
	$(CXX) -shared -o libLTF.so src/LTF.cxx -O2 -fPIC -g -Wall -std=gnu++17

root-slib:
	root -b -q src/LTF.cxx+

bin: 
	$(CXX) example_CMSinclusivejets_NN30_BRSSWpaper.cxx -o example_CMSinclusivejets_NN30_BRSSWpaper -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF  -lstdc++ $(WITH_ROOT) 
#	$(CXX) bb4l_mtop.cxx -o bb4l_mtop -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -lstdc++ 


clean:
	rm libLTF.so
