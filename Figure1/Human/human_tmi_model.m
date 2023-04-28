clc
clear


load combined_met.mat	


% Only abundant OTUs 
ctf_otu      = 0.001;% based on Brian's paper
mn_bact      = mean(motuncbispecies);
goodb        = find(mn_bact > ctf_otu);
restb = 1 - sum(motuncbispecies(:,goodb)')';
xs_b  = [motuncbispecies(:,goodb) restb];


[nSamp nOb] = size(xs_b);

for nCom = 1:10

    ZC     = 0.1*randn(nSamp,nCom);
    thetB  = 0.1*randn(nCom,nOb);

    % hyperparameters
    etaZ = 0.005;etaT = 0.005;
    grdnorm = 1;iter = 1;ctf_grad = 0.005;
    while grdnorm > ctf_grad
        QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);
        deltB = xs_b-QB;

        % gradients
        grthetB = [ZC]'*deltB;
        grz     = deltB*thetB';

        % update the variables
        ZC     = ZC - etaZ*grz;
        thetB  = thetB - etaT*grthetB;

        % errors
        grdnorm = norm(grz)/norm(ZC) + norm(grthetB)/norm(thetB);

        % Output
        if mod(iter,1000) == 0
            grdnorm
        end
        iter = iter + 1;
    end

    filen = strcat('human_',num2str(nCom),'.mat')
    save(filen,'xs_b','ZC','thetB')
end



