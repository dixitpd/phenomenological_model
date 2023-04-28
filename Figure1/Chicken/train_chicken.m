clc
clear
load ../../Data/OTUTables_IC
[nSamp nOc] = size(cecumOTU);
xs_b = cecumOTU;
for nC = 1:10
    eta = 0.002;
    thetB = 0.1*randn(nC,nOc);
    ZC  = 0.1*randn(nSamp,nC);
    filen = strcat('cecum_',num2str(nC),'.mat')
    flg   = 1;iter = 1;
    while flg > 0
        % predictions
        QB    = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);
        deltC = xs_b-QB;

        % Gradients
        grthetC     = [ZC]'*deltC;
        grZcom      = deltC*thetB';

        nrg = norm(grthetC)/norm(thetB) + norm(grZcom)/norm(ZC);

        if nrg < 0.005
            flg = 0;
        end

        % update
        thetB    = thetB - eta*grthetC;
        ZC     = ZC  - eta*grZcom;
        if mod(iter,500) == 0
            nrg
        end
        iter = iter + 1;
    end

    save(filen,'ZC','thetB','xs_b')
end
