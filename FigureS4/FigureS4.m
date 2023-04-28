clc
clear
%

load ../Data/three_kingdoms_cleaned_up

%% Only abundant OTUs in species 1 of cows
sps1 = find(cow_species==1);

ctf_otu      = 0.001;% based on Brian's paper
bacteria_xs  = bacteria_xs(sps1,:);
archea_xs    = archea_xs(sps1,:);
fungi_xs     = fungi_xs(sps1,:);
metadata_cow = metadata_cow(sps1,:);
mn_bact      = mean(bacteria_xs);
goodb        = find(mn_bact > ctf_otu);
mn_arch      = mean(archea_xs);
gooda        = find(mn_arch > ctf_otu);
mn_fung      = mean(fungi_xs);
goodf        = find(mn_fung > ctf_otu);
%
restb = 1 - sum(bacteria_xs(:,goodb)')';
xs_b  = [bacteria_xs(:,goodb) restb];
resta = 1 - sum(archea_xs(:,gooda)')';
xs_a  = [archea_xs(:,gooda) resta];
restf = 1 - sum(fungi_xs(:,goodf)')';
xs_f  = [fungi_xs(:,goodf) restf];
%
mu_met = mean(metadata_cow);st_met = std(metadata_cow);
metz   = (metadata_cow - mu_met)./st_met;

badmets = [2 4 5 6 7 8 37];
metz(:,badmets) = [];metnames(badmets) = [];
nMet = size(metz,2);

%% Filter samples based on large deviations
badsamps         = find(sum(abs(metz > 10)'));
mcow             = metadata_cow;mcow(:,badmets) = [];
mcow(badsamps,:) = [];
xs_a(badsamps,:) = [];
xs_f(badsamps,:) = [];
xs_b(badsamps,:) = [];

load ../Figure2/testing_data
xs_a_test = xs_a(testing_data,:);
xs_b_test = xs_b(testing_data,:);
xs_f_test = xs_f(testing_data,:);
xs_f_test(xs_f_test<0) = 0;
mcow_test = mcow(testing_data,:);
xs_a(testing_data,:) = [];
xs_b(testing_data,:) = [];
xs_f(testing_data,:) = [];
nTest = length(testing_data);
nSamp = size(xs_a,1);

mcow(testing_data,:) = [];
mu_met               = mean(mcow);
st_met               = std(mcow);
metz                 = (mcow-mu_met)./st_met;
metz_test            = (mcow_test-mu_met)./st_met;


%% Load
alph  = 0.01;
nC  = 20;
filen = strcat('../Data/Cows_Training/cows_',num2str(nC),'_',num2str(alph),'.mat')
load(filen)
% 

[u s v] = svd(ZC*C);
Zr  = u(:,1:nC)*s(1:nC,1:nC);v = v(:,1:nC);
Cr  = v';
rotmat1 = Cr*pinv(C);
rotmat2 = pinv(ZC)*Zr;
thetAr = inv(rotmat2)*thetA;
thetBr = inv(rotmat2)*thetB;
thetFr = inv(rotmat2)*thetF;


sK = nMet-1;

goodM = [];
for s=1:sK
    s
    for m=1:nMet
        if ismember(m,goodM)
            continue
        end
        inds = [goodM m];
        clear coeffs
        for k = 1:nC
            [mdl c]    = lasso(metz(:,inds),Zr(:,k),'lambda',0.0);
            coeffs(k,:) = mdl;
            intr(k) = c.Intercept;
        end
        Zpd = metz(:,inds)*coeffs' + intr;
        eN(m) = mean(mean(abs(Zpd - Zr)));
    end
    eN(goodM) = 1e6;
    id = find(eN==min(eN));
    goodM = [goodM id];
    fx(s) = min(eN);
end
% 
plot(fx)



