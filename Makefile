
ROOTCFLAGS    = $(shell root-config --cflags)

MYCF=""

ifneq ($(findstring m32,$(ROOTCFLAGS)),)
        MYCF= CXXFLAGS=-m32 CFLAGS=-m32
endif

ROOTLIB    = $(shell root-config --libs) -lMinuit

CXXFLAGS   = -g -I.
CXXFLAGS  += -Wno-long-long -fPIC
CXXFLAGS  += $(shell root-config --cflags)

LDFLAGS    = $(ROOTLIB)
# LDFLAGS   += -lCintex -lHistFactory -lXMLParser -lRooStats -lRooFit -lRooFitCore -lThread -lMinuit -lFoam -lHtml -lMathMore 
LDFLAGS   += -lHistFactory -lXMLParser -lRooStats -lRooFit -lRooFitCore -lThread -lMinuit -lFoam -lHtml -lMathMore 

# OBJ := $(wildcard Root/*.o)
# OBJS       = Root/Common.C Root/FitResults.C Root/NuisParameter.C Root/Sample.C Root/Systematic.C Root/TtHFit.C Root/CorrelationMatrix.C Root/NormFactor.C Root/Region.C Root/SampleHist.C Root/SystematicHist.C Root/TthPlot.C Root/HistoTools.C Root/ConfigParser.C
OBJS := $(wildcard Root/*.C)

# OBJS      += util/%.o

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
