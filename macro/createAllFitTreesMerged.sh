#! /bin/sh

root -l -b <<EOF

.L createFitTree.C++

createFitTree Higgs("recreate");
Higgs.AddFiles("results/H160_WW_2l/H160_WW_2l_*-datasetEE.root");
Higgs.SetDatasetName("Higgs");
Higgs.SaveRootTree(true);
Higgs.Init();
Higgs.Loop();

createFitTree WjetsMADGRAPH("recreate");
WjetsMADGRAPH.AddFiles("results/WjetsMADGRAPH_Fall08/WjetsMADGRAPH_Fall08-datasetEE.root");
WjetsMADGRAPH.SetDatasetName("WjetsMADGRAPH");
WjetsMADGRAPH.SaveRootTree(true);
WjetsMADGRAPH.Init();
WjetsMADGRAPH.Loop();

createFitTree ZjetsMADGRAPH("recreate");
ZjetsMADGRAPH.AddFiles("results/ZjetsMADGRAPH_Fall08/ZjetsMADGRAPH_Fall08-datasetEE.root");
ZjetsMADGRAPH.SetDatasetName("ZjetsMADGRAPH");
ZjetsMADGRAPH.SaveRootTree(true);
ZjetsMADGRAPH.Init();
ZjetsMADGRAPH.Loop();

createFitTree ttjetsMADGRAPH("recreate");
ttjetsMADGRAPH.AddFiles("results/ttjetsMADGRAPH_Fall08/ttjetsMADGRAPH_Fall08-datasetEE.root");
ttjetsMADGRAPH.SetDatasetName("ttjetsMADGRAPH");
ttjetsMADGRAPH.SaveRootTree(true);
ttjetsMADGRAPH.Init();
ttjetsMADGRAPH.Loop();

createFitTree WW_2l("recreate");
WW_2l.AddFiles("results/WW_2l/WW_2l_*-datasetEE.root");
WW_2l.SetDatasetName("WW_2l");
WW_2l.SaveRootTree(true);
WW_2l.Init();
WW_2l.Loop();

createFitTree ZZ_4l("recreate");
ZZ_4l.AddFiles("results/ZZ_4l/ZZ_4l_*-datasetEE.root");
ZZ_4l.SetDatasetName("ZZ_4l");
ZZ_4l.SaveRootTree(true);
ZZ_4l.Init();
ZZ_4l.Loop();

createFitTree WZ_3l("recreate");
WZ_3l.AddFiles("results/WZ_3l/WZ_3l_*-datasetEE.root");
WZ_3l.SetDatasetName("WZ_3l");
WZ_3l.SaveRootTree(true);
WZ_3l.Init();
WZ_3l.Loop();

createFitTree WZ_incl("recreate");
WZ_incl.AddFiles("results/WZ_incl/WZ_incl_*-datasetEE.root");
WZ_incl.SetDatasetName("WZ_incl");
WZ_incl.SaveRootTree(true);
WZ_incl.Init();
WZ_incl.Loop();

TChain other("data");
other.Add("hww_2e_Tree_WjetsMADGRAPH_21X.root");
other.Add("hww_2e_Tree_ZjetsMADGRAPH_21X.root");
other.Add("hww_2e_Tree_ZZ_4l_21X.root");
other.Add("hww_2e_Tree_WZ_incl_21X.root");
other.Add("hww_2e_Tree_WZ_3l_21X.root");
other.Merge("hww_2e_Tree_other_21X.root");

.L createDatasets.C++
createDatasetsEE();

EOF

rm -f hww*Tree*root
