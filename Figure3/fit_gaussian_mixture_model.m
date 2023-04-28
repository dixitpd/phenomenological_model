clc
clear

rval = 1e-5;
options = statset('MaxIter',5000);
nn = 20;
aa = [0.01];
for k=1:length(nn)
    nK = nn(k);
    for ax = 1:length(aa)
        alph = aa(ax);
        [nK alph]
        kxx(k,ax) = nK;
        axx(k,ax) = alph;

        filen = strcat('../Data/Cows_Training/cows_',num2str(nK),'_',num2str(alph),'.mat');
        load(filen)

        for ii = 1:50
            ii
            for nG = 1:10
                gmdl  = fitgmdist(ZC,nG,'regularizationValue',rval,'options',options);
                a1(ii,nG) = gmdl.AIC;
                a2(ii,nG) = gmdl.BIC;
                mx{ii,nG} = gmdl;
            end
        end

    end
end
save gaussian_models a1 a2 mx