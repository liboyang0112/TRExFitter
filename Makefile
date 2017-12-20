
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIB    = $(shell root-config --libs) -lMinuit

CXXFLAGS   = -g -I.
CXXFLAGS  += -Wno-long-long -fPIC
CXXFLAGS  += $(shell root-config --cflags)

LDFLAGS    = $(ROOTLIB)
LDFLAGS   += -lHistFactory -lRooStats -lRooFit -lRooFitCore

SOURCE := util/myFit.C
SOURCE += $(wildcard Root/*.C)
SOURCE += AtlasStyle.C AtlasUtils.C AtlasLabels.C

OBJS := $(SOURCE:.C=.o)

%.o: %.C
	g++ -c $(CXXFLAGS) -o $@ $<

all: myFit.exe     # this is the default executable

myFit.exe: $(OBJS)
	@echo "Linking ..."
	g++ $^ -o myFit.exe $(LDFLAGS)

clean:
	rm -rf *.exe *.o Root/*.o util/*.o

hupdate:
	g++ $(CXXFLAGS) -o hupdate hupdate.C $(ROOTCFLAGS) $(LDFLAGS)
