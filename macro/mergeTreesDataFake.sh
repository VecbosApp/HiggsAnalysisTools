mkdir results_data/merged
hadd results_data/merged/dataset_fake_ee.root results_fake_ee/DoubleElectron/DoubleElectron_*-datasetEE.root
hadd results_data/merged/dataset_fake_me.root results_fake_me/MuEG/MuEG_*-datasetME.root results_fake_me/SingleMu/SingleMu_*-datasetME.root 
