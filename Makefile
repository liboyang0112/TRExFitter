COMMONSYSTSMOOTHINGTOOLDIR =./CommonSystSmoothingTool
COMMONSTATTOOLS = ./CommonStatTools

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIB    = $(shell root-config --libs) -lMinuit

CXXFLAGS   = -g -I.
CXXFLAGS  += -Wno-long-long -fPIC
CXXFLAGS  += $(shell root-config --cflags) -Wall -pedantic

CXXFLAGSSMOOTH   = -g -I${COMMONSYSTSMOOTHINGTOOLDIR}
CXXFLAGSSMOOTH  += -Wno-long-long -fPIC
CXXFLAGSSMOOTH  += $(shell root-config --cflags) -Wall -pedantic

LDFLAGS    = $(ROOTLIB)
LDFLAGS   += -lHistFactory -lRooStats -lRooFit -lRooFitCore
LDFLAGS   += -L${COMMONSTATTOOLS}/build -lExoStats

SOURCE := util/myFit.C
SOURCE += $(wildcard Root/*.C)
SOURCE += AtlasStyle.C AtlasUtils.C AtlasLabels.C

SOURCESMOOTH := $(wildcard ${COMMONSYSTSMOOTHINGTOOLDIR}/Root/*.cxx)

OBJS := $(SOURCE:.C=.o)

OBJSSMOOTH := $(SOURCESMOOTH:.cxx=.o)


%.o: %.C
	g++ -c $(CXXFLAGS) -o $@ $<

%.o: %.cxx
	g++ -c $(CXXFLAGSSMOOTH) -o $@ $<

all: myFit.exe     # this is the default executable

myFit.exe: $(OBJS) ${OBJSSMOOTH}
	@echo "Linking ..."
	g++ $^ -o myFit.exe $(LDFLAGS)

clean:
	rm -rf *.exe *.o Root/*.o util/*.o

hupdate:
	g++ $(CXXFLAGS) -o hupdate hupdate.C $(ROOTCFLAGS) $(LDFLAGS)
