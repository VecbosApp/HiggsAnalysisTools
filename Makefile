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



NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
gGLIBS          = $(filter-out -lNew, $(NGLIBS))

CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
INCLUDEDIRCOMMON = ../
CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/
OUTLIBCOMMON     = $(INCLUDEDIRCOMMON)/CommonTools/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)HiggsBase.o: $(INCLUDEDIR)/src/HiggsBase.C $(INCLUDEDIR)/src/HiggsSelection.cc $(INCLUDEDIR)/src/ZSelection.cc $(INCLUDEDIR)/src/WSelection.cc $(INCLUDEDIR)/src/ElectronID.cc $(INCLUDEDIR)/src/plotsEleID.cc $(INCLUDEDIR)/src/ClassEfficiencyStudy.cc $(INCLUDEDIR)/src/WplusJets.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsBase.o $<
$(OUTLIBCOMMON)Conditions.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Conditions.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Conditions.o $<
$(OUTLIBCOMMON)Utils.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Utils.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIBCOMMON)Utils.o $<
$(OUTLIBCOMMON)Counters.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Counters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Counters.o $<
$(OUTLIBCOMMON)Selection.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Selection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Selection.o $<
$(OUTLIBCOMMON)EfficiencyEvaluator.o: $(INCLUDEDIRCOMMON)/CommonTools/src/EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)EfficiencyEvaluator.o $<
$(OUTLIBCOMMON)Monitor.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Monitor.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Monitor.o $<
$(OUTLIBCOMMON)SprDataFiller.o: $(INCLUDEDIRCOMMON)/CommonTools/src/SprDataFiller.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)SprDataFiller.o $<
$(OUTLIBCOMMON)TriggerMask.o: $(INCLUDEDIRCOMMON)/CommonTools/src/TriggerMask.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)TriggerMask.o $<
$(OUTLIB)HiggsSelection.o: $(INCLUDEDIR)/src/HiggsSelection.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsSelection.o $<
$(OUTLIB)ElectronID.o: $(INCLUDEDIR)/src/ElectronID.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ElectronID.o $<
$(OUTLIB)kFactorEvaluator.o: $(INCLUDEDIR)/src/kFactorEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)kFactorEvaluator.o $<
$(OUTLIB)RedHiggsTree.o: $(INCLUDEDIR)/src/RedHiggsTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedHiggsTree.o $<
$(OUTLIB)RedEWKTree.o: $(INCLUDEDIR)/src/RedEWKTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEWKTree.o $<
$(OUTLIB)RedEleIDTree.o: $(INCLUDEDIR)/src/RedEleIDTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDTree.o $<
$(OUTLIB)eleID_Higgs_Studies.o: $(INCLUDEDIR)/src/eleID_Higgs_Studies.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)eleID_Higgs_Studies.o $<
$(OUTLIB)CutBasedEleIDSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/CutBasedEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)CutBasedEleIDSelector.o $<
$(OUTLIB)CommonHiggsPreselector.o: $(INCLUDEDIRCOMMON)/HiggsAnalysisTools/src/CommonHiggsPreselector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)CommonHiggsPreselector.o $<
$(OUTLIB)CutBasedHiggsSelector.o: $(INCLUDEDIRCOMMON)/HiggsAnalysisTools/src/CutBasedHiggsSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)CutBasedHiggsSelector.o $<

#----------------------------------------------------#

# ==================== HiggsApp =============================================
HiggsApp:  $(INCLUDEDIR)/src/HiggsApp.C $(OUTLIB)HiggsBase.o $(OUTLIBCOMMON)Conditions.o $(OUTLIBCOMMON)Selection.o $(OUTLIBCOMMON)EfficiencyEvaluator.o $(OUTLIBCOMMON)Counters.o $(OUTLIBCOMMON)Monitor.o $(OUTLIBCOMMON)SprDataFiller.o $(OUTLIBCOMMON)TriggerMask.o $(OUTLIBCOMMON)Utils.o $(OUTLIB)kFactorEvaluator.o $(OUTLIB)RedHiggsTree.o $(OUTLIB)RedEWKTree.o $(OUTLIB)CutBasedEleIDSelector.o $(OUTLIB)CommonHiggsPreselector.o $(OUTLIB)CutBasedHiggsSelector.o 
	$(CXX) $(CXXFLAGS) -o HiggsApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(GLIBS) $ $<
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
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f HiggsApp
	rm -f ReducedTree_HwwEleId
	rm -f eleID_Higgs_Studies

all:  HiggsApp eleID_Higgs_Studies ReducedTree_HwwEleId
