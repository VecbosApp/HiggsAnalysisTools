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


$(OUTLIB)HiggsBase.o: $(INCLUDEDIR)/src/HiggsBase.C $(INCLUDEDIR)/src/HiggsSelection.cc $(INCLUDEDIR)/src/HiggsEleIdOptimToyMC.cc $(INCLUDEDIR)/src/RedEleIDOptimTree.cc $(INCLUDEDIR)/src/RedLikeOptimTree.cc $(INCLUDEDIR)/src/HiggsIsolationOptimToyMC.cc $(INCLUDEDIR)/src/RedIsolationOptimTree.cc $(INCLUDEDIR)/src/ZplusJetsSelection.cc $(INCLUDEDIR)/src/LeptonPlusFakeSelection.cc $(INCLUDEDIR)/src/HiggsVertexing.cpp $(INCLUDEDIR)/src/VertexTree.cc
# no more used ------------------------------------
#$(INCLUDEDIR)/src/ZSelection.cc $(INCLUDEDIR)/src/WSelection.cc 
#$(INCLUDEDIR)/src/ElectronID.cc $(INCLUDEDIR)/src/plotsEleID.cc 
#$(INCLUDEDIR)/src/ClassEfficiencyStudy.cc $(INCLUDEDIR)/src/WplusJets.cc 
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
$(OUTLIB)CutBasedEleIDSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/CutBasedEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)CutBasedEleIDSelector.o $<
$(OUTLIB)HiggsSelection.o: $(INCLUDEDIR)/src/HiggsSelection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsSelection.o $<
$(OUTLIB)HiggsEleIdOptimToyMC.o: $(INCLUDEDIR)/src/HiggsEleIdOptimToyMC.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsEleIdOptimToyMC.o $<
$(OUTLIB)RedEleIDOptimTree.o: $(INCLUDEDIR)/src/RedEleIDOptimTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDOptimTree.o $<
$(OUTLIB)RedLikeOptimTree.o: $(INCLUDEDIR)/src/RedLikeOptimTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedLikeOptimTree.o $<
$(OUTLIB)HiggsIsolationOptimToyMC.o: $(INCLUDEDIR)/src/HiggsIsolationOptimToyMC.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsIsolationOptimToyMC.o $<
$(OUTLIB)HiggsKinematicsOptimToyMC.o: $(INCLUDEDIR)/src/HiggsKinematicsOptimToyMC.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsKinematicsOptimToyMC.o $<
$(OUTLIB)RedIsolationOptimTree.o: $(INCLUDEDIR)/src/RedIsolationOptimTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedIsolationOptimTree.o $<
$(OUTLIB)kFactorEvaluator.o: $(INCLUDEDIR)/src/kFactorEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)kFactorEvaluator.o $<
$(OUTLIB)RedHiggsTree.o: $(INCLUDEDIR)/src/RedHiggsTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedHiggsTree.o $<
$(OUTLIB)RedTriggerTree.o: $(INCLUDEDIR)/src/RedTriggerTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedTriggerTree.o $<
$(OUTLIB)CommonHiggsPreselector.o: $(INCLUDEDIR)/src/CommonHiggsPreselector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CommonHiggsPreselector.o $<
$(OUTLIB)CutBasedHiggsSelector.o: $(INCLUDEDIR)/src/CutBasedHiggsSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CutBasedHiggsSelector.o $<
$(OUTLIB)ZplusJetsSelection.o: $(INCLUDEDIR)/src/ZplusJetsSelection.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ZplusJetsSelection.o $<
$(OUTLIB)LeptonPlusFakeSelection.o: $(INCLUDEDIR)/src/LeptonPlusFakeSelection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LeptonPlusFakeSelection.o $<
$(OUTLIB)HiggsVertexing.o: $(INCLUDEDIR)/src/HiggsVertexing.cpp
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsVertexing.o $<
$(OUTLIB)VertexTree.o: $(INCLUDEDIR)/src/VertexTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VertexTree.o $<

# no more used ------------------------------------
#$(OUTLIB)ElectronID.o: $(INCLUDEDIR)/src/ElectronID.C
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ElectronID.o $<
#$(OUTLIB)RedEWKTree.o: $(INCLUDEDIR)/src/RedEWKTree.cc
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEWKTree.o $<
#$(OUTLIB)RedEleIDTree.o: $(INCLUDEDIR)/src/RedEleIDTree.cc
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDTree.o $<
#$(OUTLIB)eleID_Higgs_Studies.o: $(INCLUDEDIR)/src/eleID_Higgs_Studies.cc
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)eleID_Higgs_Studies.o $<

#----------------------------------------------------#

# ==================== HiggsApp =============================================
HiggsApp:  $(INCLUDEDIR)/src/HiggsApp.C $(OUTLIB)HiggsBase.o $(OUTLIBCOMMON)Conditions.o $(OUTLIBCOMMON)Selection.o $(OUTLIBCOMMON)EfficiencyEvaluator.o $(OUTLIBCOMMON)Counters.o $(OUTLIBCOMMON)Monitor.o $(OUTLIBCOMMON)SprDataFiller.o $(OUTLIBCOMMON)TriggerMask.o $(OUTLIBCOMMON)Utils.o $(OUTLIB)kFactorEvaluator.o $(OUTLIB)RedHiggsTree.o $(OUTLIB)RedTriggerTree.o $(OUTLIB)RedEleIDOptimTree.o $(OUTLIB)RedLikeOptimTree.o $(OUTLIB)RedIsolationOptimTree.o $(OUTLIB)CutBasedEleIDSelector.o $(OUTLIB)CommonHiggsPreselector.o $(OUTLIB)CutBasedHiggsSelector.o $(OUTLIB)LeptonPlusFakeSelection.o $(OUTLIB)VertexTree.o
# no more used ------------------------------------
#$(OUTLIB)RedEWKTree.o 
	$(CXX) $(CXXFLAGS) -o HiggsApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(GLIBS) $ $<
HiggsApp.clean:
	rm -f HiggsApp

# ==================== eleID =============================================
eleID_Higgs_Studies:  $(INCLUDEDIR)eleID_Higgs_Studies.cpp $(OUTLIB)RedEleIDTree.o 
	$(CXX) $(CXXFLAGS) -o eleID_Higgs_Studies $(OUTLIB)/*.o $(GLIBS) $ $<
#eleID_Higgs_Studies.clean:
#	rm -f eleID_Higgs_Studies

eleIDtableToy:  $(INCLUDEDIR)/src/eleIDtableToy.cpp
	$(CXX) $(CXXFLAGS) -o eleIDtableToy $(GLIBS) $ $<
eleIDtoyPlot_input:  $(INCLUDEDIR)/src/eleIDtoyPlot_input.cpp
	$(CXX) $(CXXFLAGS) -o eleIDtoyPlot_input $(GLIBS) $ $<
eleIDtableToy.clean:
	rm -f eleIDtableToy
eleIDtoyPlot_input.clean:
	rm -f eleIDtoyPlot_input
eleIDtableLike:  $(INCLUDEDIR)/src/eleIDtableLike.cpp
	$(CXX) $(CXXFLAGS) -o eleIDtableLike $(GLIBS) $ $<
eleIDtableLike.clean:
	rm -f eleIDtableLike
isolationTableToy:  $(INCLUDEDIR)/src/isolationTableToy.cpp
	$(CXX) $(CXXFLAGS) -o isolationTableToy $(GLIBS) $ $<
isolationTableToy.clean:
	rm -f isolationTableToy
isolationStudies_input:  $(INCLUDEDIR)/src/isolationStudies_input.cpp
	$(CXX) $(CXXFLAGS) -o isolationStudies_input $(GLIBS) $ $<
isolationStudies_input.clean:
	rm -f isolationStudies_input
kinematicsTableToy:  $(INCLUDEDIR)/src/kinematicsTableToy.cpp
	$(CXX) $(CXXFLAGS) -o kinematicsTableToy $(GLIBS) $ $<
kinematicsTableToy.clean:
	rm -f kinematicsTableToy
vtxAndIsoOptim:  $(INCLUDEDIR)/src/vtxAndIsoOptim.cpp
	$(CXX) $(CXXFLAGS) -o vtxAndIsoOptim $(GLIBS) $ $<
vtxAndIsoOptim.clean:
	rm -f vtxAndIsoOptim

# ==================== reduced trees =============================================
#ReducedTree_HwwEleId:  $(INCLUDEDIR)ReducedTree_HwwEleId.cpp $(OUTLIB)RedEleIDTree.o 
#	$(CXX) $(CXXFLAGS) -o ReducedTree_HwwEleId $(OUTLIB)/*.o $(GLIBS) $ $<
#ReducedTree_HwwEleId.clean:
#	rm -f ReducedTree_HwwEleId

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f HiggsApp
# rm -f ReducedTree_HwwEleId
# rm -f eleID_Higgs_Studies
	rm -f eleIDtableToy
	rm -f eleIDtoyPlot_input
	rm -f isolationStudies_input
	rm -f eleIDtableLike
	rm -f isolationTableToy
	rm -f kinematicsTableToy
	rm -f vtxAndIsoOptim

all:  HiggsApp
