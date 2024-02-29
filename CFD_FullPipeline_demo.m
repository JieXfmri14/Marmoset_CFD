clc
clear all

% set code folder path
mypath='E:\Marmoset_CFD-main\'; 
addpath(genpath(strcat(mypath,'functions_matlab')))

%% Step1: compute and decompose Graph Laplacian 
[W,U]= GSP_Laplacian(mypath);

%% project BOLD-fMRI data onto CC eigenmodes, and split low- and high-frequency eigenmodes
[X_RS,zX_RS,Vlow,Vhigh]= find_filter_cut_off(mypath,U);

%% Step2: the accuracy of the reconstruction the marmoset activity
[CC_recon_act_marmoset,CC_recon_FC_marmoset]= CC_recon_marmoset_activity(W,X_RS,U);

%% Step3: analyze graph signal,i.e., filted low- and high-frequency components
[N_low,N_high,mean_low,mean_high] = GSanalysis(zX_RS,Vhigh,Vlow,U);

%% Step4: system permutation test
[net_low_sig,net_high_sig] = System_permutation_test(zX_RS,U,Vlow,Vhigh,mean_low,mean_high);

%% Step5: compute CFD for empirical BOLD-fMRI data and test significance
[CFDlog,CFD_thr,CFD_surr] = Calcalate_CFD(zX_RS,U,Vlow,Vhigh,mean_low,mean_high,N_low,N_high);

%% Step6: generalizing of marmoset-derived eigenmodes to human
[CC_recon_act_human,CC_recon_FC_human]= Generalizing_eigenmodes_human(mypath);


