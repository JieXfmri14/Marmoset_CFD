function [CFD_log,CFD_thr_log,CFD_surr_log] = Calcalate_CFD(zX_RS,U,Vlow,Vhigh,mean_low,mean_high,N_low,N_high)
% Compute CFD for empirical data and test significance 

% =========================================================================
%     1) Compute cellular-functional decoupling (CFD) of empirical fMRI data   
% =========================================================================
mean_CFD=mean_high./mean_low;          %empirical group-level CFD
CFD=N_high./N_low;                     %emipirical individual CFD

% =========================================================================
%     2) Create CC-informed graph signal surrogates by graph spectral randomization (SR)
% =========================================================================
nSurr = 19;    % significance level = 1/(nSurr+1)
data = zX_RS;
for s=1:size(data,3)
    X=data(:,:,s);
    for n=1:nSurr
        %randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(U,1),1));
        PHIdiag(PHIdiag==0)=-1;
        PHI=diag(PHIdiag);                  
        fdata{s}(:,:,n)=U*PHI*(U'*X);   % graph signal surrogates       
    end
end

% =========================================================================
%     3) find filter low-and high-frequency signal in surrogates
% =========================================================================      
n_ROI = size(U,1);
for s=1:size(fdata,2)   
    for i=1:size(fdata{1},3)  
        X_hat_surr{s}(:,:,i)=U'*fdata{s}(:,:,i);
        X_c_surr{s}(:,:,i)=Vlow*X_hat_surr{s}(:,:,i);     
        X_d_surr{s}(:,:,i)=Vhigh*X_hat_surr{s}(:,:,i);
        for r=1:size(fdata{1},1)
            N_c_surr(r,i,s)=norm(X_c_surr{s}(r,:,i));
            N_d_surr(r,i,s)=norm(X_d_surr{s}(r,:,i));
        end
    end
end

% =========================================================================
%     4) test significance of CFD
% =========================================================================      
% consider and test mean across subjects
CFD_surr=N_d_surr./N_c_surr;                      % surrogates CFD for every subject 
CFD_surr_avgsurr=squeeze(mean(CFD_surr,2));       
CFD_surr_avgsurrsubjs=mean(CFD_surr_avgsurr,2);    

% find threshold for max
for s=1:size(CFD_surr,3)
       max_CFD_surr(:,s)=max(CFD_surr(:,:,s)')';
end
% find threshold for min
for s=1:size(CFD_surr,3)
      min_CFD_surr(:,s)=min(CFD_surr(:,:,s)')';
end

% =========================================================================
%     5) select significant CFD for each subject, across  surrogates
% =========================================================================      
% for each subject, I threshold the ratio based on individual ratio's surrogate distribution 
for s=1:size(fdata,2) 
    significant_values_max(s)=size(find(CFD(:,s)>max_CFD_surr(:,s)),1);
    significant_values_min(s)=size(find(CFD(:,s)<min_CFD_surr(:,s)),1);
    CFD_thr_max(:,s)=CFD(:,s)>max_CFD_surr(:,s);
    CFD_thr_min(:,s)=CFD(:,s)<min_CFD_surr(:,s);
    detect_max=sum(CFD_thr_max'); % amounts of detection per region
    detect_min=sum(CFD_thr_min');
end

% for every region, test across subjects 0.05, correcting for the number of tests (regions), 0.05/n_ROI
x=0:1:100;
y=binocdf(x,100,0.05,'upper');            
THRsubjects=x(min(find(y<0.05/n_ROI))); 
THRsubjects=floor(size(fdata,2)/100*THRsubjects)+1;

CFD_sig_higher=detect_max>THRsubjects;
CFD_sig_lower=detect_min>THRsubjects;

CFD_sig_higher_positions=find(CFD_sig_higher==1);
CFD_sig_lower_positions=find(CFD_sig_lower==1);

CFD_sig_tot_positions=[find(CFD_sig_higher==1),find(CFD_sig_lower==1)];
CFD_sig_tot_positions=sort(unique(CFD_sig_tot_positions));

% threshold empirical mean ratios
mean_CFD_thr=ones(n_ROI,1);
mean_CFD_thr(CFD_sig_tot_positions)=mean_CFD(CFD_sig_tot_positions);


%% surrogates CFD pattern
saturate=1;
CC2=log2(CFD_surr_avgsurrsubjs); 
%% adjust CFD values for saturation (to eliminate outliers peaks)
if saturate
    thr=1;
    CC2new=CC2;
    CC2new(find(CC2>thr))=0;
    CC2new(find(CC2>thr))=max(CC2new);
    CC2new(find(CC2<-thr))=0;
    CC2new(find(CC2<-thr))=min(CC2new);
    CC2=CC2new;
end
CFD_surr_log = CC2;

%% emprical mean_CFD
saturate=1;
CC2=log2(mean_CFD); 
if saturate
    thr=1;
    CC2new=CC2;
    CC2new(find(CC2>thr))=0;
    CC2new(find(CC2>thr))=max(CC2new);
    CC2new(find(CC2<-thr))=0;
    CC2new(find(CC2<-thr))=min(CC2new);
    CC2=CC2new;
end
CFD_log = CC2;


%% Significant decoupling and coupling regions
saturate=1;
CC2=log2(mean_CFD_thr); 
if saturate
    thr=1;
    CC2new=CC2;
    CC2new(find(CC2>thr))=0;
    CC2new(find(CC2>thr))=max(CC2new);
    CC2new(find(CC2<-thr))=0;
    CC2new(find(CC2<-thr))=min(CC2new);
    CC2=CC2new;
end
CFD_thr_log =  CC2;

%% Results: the CFD pattern
Output = [CFD_log,CFD_thr_log,CFD_surr_log];
xlswrite('Result/CFD_index/CFD_pattern.xlsx',Output,'B2:D56')

% =========================================================================
%     6) network-level analysis
% =========================================================================      
network_level = xlsread('marmoset_brain_template/marmoset_55Nodes_11Classes.xlsx','sheet1','D2:D56');
net_num = length(unique(network_level));
%%% null distribution of 11 classes
network_mean_CFD = zeros(net_num,1);
network_median_CFD = zeros(net_num,1);

network = ["VC","PPC","PCC","LIT","AU","SS",...
     "MOT","mPFC","OFC","VLPFC","DLPFC"];

for i=1:net_num
        temp = find(network_level==i);
        network_mean_CFD(i)=mean(CFD_log(temp));
        network_median_CFD(i)=median(CFD_log(temp));
end

end









