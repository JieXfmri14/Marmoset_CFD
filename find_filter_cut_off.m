function [X_RS,zX_RS,Vlow,Vhigh]= find_filter_cut_off(mypath,U)

% % =========================================================================
% %  (1) In main text: The cut-off of CC eigenmodes using the graph spectrum dichotomy approach  
% % =========================================================================
% load functional time series, size: n_ROI * timepoints * subjects
X_RS=load(strcat(mypath,'Data/marmoset_BOLD_demo.mat'));
X_RS = X_RS.Mar_RS;
nsubjs_RS=size(X_RS,3);
n_ROI = size(U,2);  

%% the mean, std and energy distribution 
%=====================================================================
% If all regions have a similar amplitude (or at least std dev), then no rescaling is needed and no bias is introduced. 
%====================================================================
for sub =1:nsubjs_RS
    for i =1:n_ROI
        mu(i,sub) = mean(X_RS(i,:,sub),2);
        sigma(i,sub) = std(X_RS(i,:,sub),0,2);
        energy(i,sub) = norm(X_RS(i,:,sub));
    end
end

% Normalized fMRI timecourses (If all regions have similar amplitude, no need for normalization)
zX_RS=zscore(X_RS,0,2);

% rs-fMRI data projected on the CC eigenmode
clear X_hat_L  
for s=1:nsubjs_RS
    X_hat_L(:,:,s)=U'*zX_RS(:,:,s);              
end

% energy spectral density 
pow=abs(X_hat_L).^2;          
pow_group = mean(pow,3);      
ESD=squeeze(mean(pow,2));     

% cutoff frequency C: graph spectrum dichotomy approach
mESD=mean(ESD,2);                        
AUCTOT=trapz(mESD(1:n_ROI));               

i=0;
AUC=0;
while AUC<AUCTOT/2
    AUC=trapz(mESD(1:i));      
    i=i+1;
end

NN=i-1          
NNL=n_ROI-NN;   

% split CC eigenmodes in low/high frequency
Vlow=zeros(size(U));
Vhigh=zeros(size(U));
Vlow(:,1:NN)=U(:,1:NN);                 
Vhigh(:,NN+1:end)=U(:,NN+1:end);       

% % =========================================================================
% %  (2) Stability_analysis: divide the CC eigenmodes into low, medium, and high components
% %      Chose different lowest KL and highest KH CC eigenmodes 
% % =========================================================================
% for example: KL = 10 and KH = 35
% Vlow=zeros(size(U));
% Vhigh=zeros(size(U));
% KL = 10;
% KH = 35;
% Vlow(:,1:KL)=U(:,1:KL);                      % low-frequency eigenmodes
% Vhigh(:,end-KH+1:end)=U(:,end-KH+1:end);     % high-frequency eigenmodes

end
