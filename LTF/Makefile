

all: slib bin

bin: 
	$(CXX) bb4l_mtop.cxx -o bb4l_mtop -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -lstdc++ 

example:
	$(CXX) example_LTF.cxx -o example_LTF -O2 -g -Wall -Wno-sign-compare -std=gnu++17 -L. -lLTF $(shell root-config --libs) -I$(shell root-config --incdir) -lstdc++ 

slib:
	$(CXX) -shared -o libLTF.so LTF.cxx -O2 -fPIC -g -Wall -std=gnu++17

root-slib:
	root -b -q LTF.cxx+

clean:
	rm libLTF.so
