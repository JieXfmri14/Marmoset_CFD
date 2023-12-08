function [XrandS] = SR_create_CCinformed_surrogates(data,U,nSurr)

%% Create CC-informed graph signal surrogates by randomization of real eigenmodes Fourier coefficients
clear XrandS 
X_hat_rand=zeros(size(data,1),size(data,2),size(data,3),nSurr); 
for s=1:size(data,3)
    X=data(:,:,s);
    for n=1:nSurr
        %randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(U,1),1));
        PHIdiag(PHIdiag==0)=-1;
        PHI=diag(PHIdiag);                  
        XrandS{s}(:,:,n)=U*PHI*(U'*X);         
    end
end

end

