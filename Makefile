#ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
#ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
#ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

#PG da qui per macosx
#PG -----------------

ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
MACOSX_MINOR := $(shell sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2)
CXX           = c++ -lm
#CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual $(ROOTCFLAGS) $(CLHEPCFLAGS)
CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE  
LD           = c++
LDFLAGS       = -O2 -g
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
ifeq ($(MACOSX_MINOR),4)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.4 c++
else
ifeq ($(MACOSX_MINOR),3)
UNDEFOPT      = dynamic_lookup
LD            = MACOSX_DEPLOYMENT_TARGET=10.3 c++
else
UNDEFOPT      = suppress
LD            = c++
endif
endif

endif

#PG fin qui
#PG -------


NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR    = ./
CXX	      += -I$(INCLUDEDIR) -I.
OUTLIB	      = ./lib/

.SUFFIXES: .cc,.C
.PREFIXES: ./lib/


$(OUTLIB)HiggsBase.o: $(INCLUDEDIR)HiggsBase.C $(INCLUDEDIR)HiggsSelection.cc $(INCLUDEDIR)ZSelection.cc $(INCLUDEDIR)WSelection.cc $(INCLUDEDIR)ElectronID.cc $(INCLUDEDIR)plotsEleID.cc $(INCLUDEDIR)ClassEfficiencyStudy.cc $(INCLUDEDIR)WplusJets.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsBase.o $<
$(OUTLIB)HiggsSelection.o: $(INCLUDEDIR)HiggsSelection.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsSelection.o $<
$(OUTLIB)ElectronID.o: $(INCLUDEDIR)ElectronID.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ElectronID.o $<
$(OUTLIB)Counters.o: $(INCLUDEDIR)Counters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Counters.o $<
$(OUTLIB)Selection.o: $(INCLUDEDIR)Selection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Selection.o $<
$(OUTLIB)EfficiencyEvaluator.o: $(INCLUDEDIR)EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)EfficiencyEvaluator.o $<
$(OUTLIB)Monitor.o: $(INCLUDEDIR)Monitor.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Monitor.o $<
$(OUTLIB)SprDataFiller.o: $(INCLUDEDIR)SprDataFiller.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SprDataFiller.o $<
$(OUTLIB)kFactorEvaluator.o: $(INCLUDEDIR)kFactorEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)kFactorEvaluator.o $<
$(OUTLIB)RedHiggsTree.o: $(INCLUDEDIR)RedHiggsTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedHiggsTree.o $<
$(OUTLIB)RedEWKTree.o: $(INCLUDEDIR)RedEWKTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEWKTree.o $<

$(OUTLIB)RedEleIDTree.o: $(INCLUDEDIR)RedEleIDTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDTree.o $<

$(OUTLIB)eleID_Higgs_Studies.o: $(INCLUDEDIR)eleID_Higgs_Studies.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)eleID_Higgs_Studies.o $<

#----------------------------------------------------#

# ==================== HiggsApp =============================================
HiggsApp:  $(INCLUDEDIR)HiggsApp.C $(OUTLIB)HiggsBase.o $(OUTLIB)Selection.o $(OUTLIB)EfficiencyEvaluator.o $(OUTLIB)Counters.o $(OUTLIB)Monitor.o $(OUTLIB)SprDataFiller.o $(OUTLIB)kFactorEvaluator.o $(OUTLIB)RedHiggsTree.o $(OUTLIB)RedEWKTree.o
	$(CXX) $(CXXFLAGS) -o HiggsApp $(OUTLIB)/*.o $(GLIBS) $ $<
HiggsApp.clean:
	rm -f HiggsApp

# ==================== eleID =============================================
eleID_Higgs_Studies:  $(INCLUDEDIR)eleID_Higgs_Studies.cpp $(OUTLIB)RedEleIDTree.o 
	$(CXX) $(CXXFLAGS) -o eleID_Higgs_Studies $(OUTLIB)/*.o $(GLIBS) $ $<
eleID_Higgs_Studies.clean:
	rm -f eleID_Higgs_Studies

# ==================== reduced trees =============================================
ReducedTree_HwwEleId:  $(INCLUDEDIR)ReducedTree_HwwEleId.cpp $(OUTLIB)RedEleIDTree.o 
	$(CXX) $(CXXFLAGS) -o ReducedTree_HwwEleId $(OUTLIB)/*.o $(GLIBS) $ $<
ReducedTree_HwwEleId.clean:
	rm -f ReducedTree_HwwEleId

clean:
	rm -f $(OUTLIB)*.o
	rm -f HiggsApp
	rm -f ReducedTree_HwwEleId
	rm -f eleID_Higgs_Studies

all:  HiggsApp eleID_Higgs_Studies ReducedTree_HwwEleId
