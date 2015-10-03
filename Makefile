
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIB    = $(shell root-config --libs) -lMinuit

CXXFLAGS   = -g -I.
CXXFLAGS  += -Wno-long-long -fPIC
CXXFLAGS  += $(shell root-config --cflags)

LDFLAGS    = $(ROOTLIB)
LDFLAGS   += -lHistFactory -lRooStats -lRooFit -lRooFitCore

OBJS := $(wildcard Root/*.C)
OBJS      += AtlasStyle.C AtlasUtils.C AtlasLabels.C

# external stuff
# can use RootCore, but keep like this to make it portable
# OBJS += ../TopDataPreparation/Root/SampleXsection.o
# CXXFLAGS += -I../RootCoreBin/include
CXXFLAGS   += -I${ROOTSYS}/include -L${ROOTSYS}/lib

%.o: %.C
	g++ -c $(CXXFLAGS) -o $@ $<

all: myFit.exe     # this is the default executable

%.exe: util/%.o $(OBJS) 
	echo $@
	g++ $(CXXFLAGS) -o $@ $(patsubst %.exe,%.o,util/$@) $(OBJS) $(LDFLAGS)
clean:
	rm -rf *.exe *.o Root/*.o util/*.o
