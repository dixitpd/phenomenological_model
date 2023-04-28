clc
clear
%
rval = 1e-5;

%% Chicken
nK =  5;
filen = strcat('Chicken/cecum_',num2str(nK),'.mat')
load(filen)
nSamp = size(ZC,1);
options = statset('MaxIter',5000);


for ii = 1:20
    ii
    for nG = 1:10
        gmdl  = fitgmdist(ZC,nG,'regularizationValue',rval,'options',options);
        a1(ii,nG) = gmdl.AIC;
        a2(ii,nG) = gmdl.BIC;
        mx{ii,nG} = gmdl;
    end
end
save gaussian_models_chicken mx a1 a2 


%% Cows
nK =  5;
filen = strcat('Cows/cows_',num2str(nK),'.mat')
load(filen)
nSamp = size(ZC,1);
options = statset('MaxIter',5000);

for ii = 1:20
    ii
    for nG = 1:10
        gmdl  = fitgmdist(ZC,nG,'regularizationValue',rval,'options',options);
        a1(ii,nG) = gmdl.AIC;
        a2(ii,nG) = gmdl.BIC;
        mx{ii,nG} = gmdl;
    end
end
save gaussian_models_cows mx a1 a2 

%% Human
nK =  5;
filen = strcat('Human/human_',num2str(nK),'.mat')
load(filen)
nSamp = size(ZC,1);
options = statset('MaxIter',5000);

for ii = 1:20
    ii
    for nG = 1:10
        gmdl  = fitgmdist(ZC,nG,'regularizationValue',rval,'options',options);
        a1(ii,nG) = gmdl.AIC;
        a2(ii,nG) = gmdl.BIC;
        mx{ii,nG} = gmdl;
    end
end
save gaussian_models_human mx a1 a2

