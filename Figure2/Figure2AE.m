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

load testing_data
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
goodMs = goodM;

goodM = goodM(1:10);
clear coeffs
for k = 1:nC
    [mdl c]    = lasso(metz(:,goodM),ZC(:,k),'lambda',0.0);
    coeffs(k,:) = mdl;
    intr(k) = c.Intercept;
end
Zp = metz_test(:,goodM)*coeffs' + intr;
QB_t = exp(-Zp*thetB);QB_t = normalize(QB_t,2,'norm',1);
QA_t = exp(-Zp*thetA);QA_t = normalize(QA_t,2,'norm',1);
QF_t = exp(-Zp*thetF);QF_t = normalize(QF_t,2,'norm',1);
klA_10 = JSD(xs_a_test,QA_t);klB_10 = JSD(xs_b_test,QB_t);klF_10 = JSD(xs_f_test,QF_t);
bcB_10 = bc_pairs(xs_b_test,QB_t);bcA_10 = bc_pairs(xs_a_test,QA_t);bcF_10 = bc_pairs(xs_f_test,QF_t);
% 
% 

'Best fit model'
Zp0 = 0.1*randn(nTest,nC);
%

grdnorm = 1;eta = 0.005;
iter = 1;ctf_grad = 0.01;
while grdnorm > ctf_grad
    QA = exp(-[Zp0 ]*thetA);QA = normalize(QA,2,'norm',1);
    QB = exp(-[Zp0 ]*thetB);QB = normalize(QB,2,'norm',1);
    QF = exp(-[Zp0 ]*thetF);QF = normalize(QF,2,'norm',1);
    deltA = xs_a_test-QA;
    deltB = xs_b_test-QB;
    deltF = xs_f_test-QF;
   
    grzQA   = deltA*thetA';
    grzQB   = deltB*thetB';
    grzQF   = deltF*thetF';
    grz     = grzQA + grzQB + grzQF;

    % update the variables
    Zp0     = Zp0 - eta*grz;
    grdnorm = norm(grz)/norm(Zp0); 
    % Output
    if mod(iter,100) == 0
        grdnorm
    end
    iter = iter + 1;
end


QB_t = exp(-Zp0*thetB);QB_t = normalize(QB_t,2,'norm',1);
QA_t = exp(-Zp0*thetA);QA_t = normalize(QA_t,2,'norm',1);
QF_t = exp(-Zp0*thetF);QF_t = normalize(QF_t,2,'norm',1);
klA0 = JSD(xs_a_test,QA_t);klB0 = JSD(xs_b_test,QB_t);klF0 = JSD(xs_f_test,QF_t);
bcB0 = bc_pairs(xs_b_test,QB_t);bcA0 = bc_pairs(xs_a_test,QA_t);bcF0 = bc_pairs(xs_f_test,QF_t);


'JSDM'

ZC_met = metz(:,goodM);sK = size(ZC_met,2);
nB = size(xs_b,2);nA = size(xs_a,2);nF = size(xs_f,2);
thxb = randn(sK,nB);thxa = randn(sK,nA);thxf = randn(sK,nF);
grdnorm = 1;etaT = 0.005;
iter = 1;ctf_grad = 0.01;
while grdnorm > ctf_grad
    QA = exp(-[ZC_met]*thxa);QA = normalize(QA,2,'norm',1);
    QB = exp(-[ZC_met]*thxb);QB = normalize(QB,2,'norm',1);
    QF = exp(-[ZC_met]*thxf);QF = normalize(QF,2,'norm',1);
    deltA = xs_a-QA;
    deltB = xs_b-QB;
    deltF = xs_f-QF;

    % gradients
    grthetA = [ZC_met]'*deltA;
    grthetB = [ZC_met]'*deltB;
    grthetF = [ZC_met]'*deltF;


    % update the variables
    thxa  = thxa - etaT*grthetA;
    thxb  = thxb - etaT*grthetB;
    thxf  = thxf - etaT*grthetF;


    % errors
    grdnorm = norm(grthetA)/norm(thxa) + norm(grthetB)/norm(thxb) + norm(grthetF)/norm(thxf);

    % Output
    if mod(iter,100) == 0
        grdnorm
    end
    iter = iter + 1;

end


Zp   = metz_test(:,goodM);
QB_t_jsdm = exp(-Zp*thxb);QB_t_jsdm = normalize(QB_t_jsdm,2,'norm',1);
QA_t_jsdm = exp(-Zp*thxa);QA_t_jsdm = normalize(QA_t_jsdm,2,'norm',1);
QF_t_jsdm = exp(-Zp*thxf);QF_t_jsdm = normalize(QF_t_jsdm,2,'norm',1);
%
klA_JSDM = JSD(xs_a_test,QA_t_jsdm);
klB_JSDM = JSD(xs_b_test,QB_t_jsdm);
klF_JSDM = JSD(xs_f_test,QF_t_jsdm);
bcB_JSDM = bc_pairs(xs_b_test,QB_t_jsdm);
bcA_JSDM = bc_pairs(xs_a_test,QA_t_jsdm);
bcF_JSDM = bc_pairs(xs_f_test,QF_t_jsdm);


% Arcaea KL
subplot(2,9,1)
ylim([-0.05 1])
s = violinplot(klA_10);
s.ViolinColor = [0.5 0.5 0];
set(gca,'FontSize',15)

subplot(2,9,2)
ylim([-0.05 1])
s = violinplot(klA0);
s.ViolinColor = [0.5 0.5 0];

subplot(2,9,3)
s = violinplot(full(klA_JSDM))
s.ViolinColor = [0.5 0.5 0];
ylim([-0.05 1])

% Bacteria KL
subplot(2,9,4)
ylim([-0.05 1])
s = violinplot(klB_10);
s.ViolinColor = [0.5 0 0.5];
set(gca,'FontSize',15)

subplot(2,9,5)
ylim([-0.05 1])
s = violinplot(klB0);
s.ViolinColor = [0.5 0 0.5];

subplot(2,9,6)
s = violinplot(full(klB_JSDM))
s.ViolinColor = [0.5 0 0.5];
ylim([-0.05 1])

% Fungal KL
subplot(2,9,7)
ylim([-0.05 1.5])
s = violinplot(klF_10);
s.ViolinColor = [0 0.5 0.5];
set(gca,'FontSize',15)

subplot(2,9,8)
ylim([-0.05 1.5])
s = violinplot(klF0);
s.ViolinColor = [0 0.5 0.5];

subplot(2,9,9)
s = violinplot(full(klF_JSDM))
s.ViolinColor = [0 0.5 0.5];
ylim([-0.05 1.5])



% Arcaea BC
subplot(2,9,10)
ylim([-0.05 1])
s = violinplot(bcA_10);
s.ViolinColor = [0.5 0.5 0];
set(gca,'FontSize',15)

subplot(2,9,11)
ylim([-0.05 1])
s = violinplot(bcA0);
s.ViolinColor = [0.5 0.5 0];

subplot(2,9,12)
s = violinplot(full(bcA_JSDM))
s.ViolinColor = [0.5 0.5 0];
ylim([-0.05 1])



% Bacteria BC
subplot(2,9,13)
ylim([-0.05 1])
s = violinplot(bcB_10);
s.ViolinColor = [0.5 0 0.5];
set(gca,'FontSize',15)

subplot(2,9,14)
ylim([-0.05 1])
s = violinplot(bcB0);
s.ViolinColor = [0.5 0 0.5];

subplot(2,9,15)
s = violinplot(full(bcB_JSDM))
s.ViolinColor = [0.5 0 0.5];
ylim([-0.05 1])

% Fungal BC
subplot(2,9,16)
ylim([-0.05 1.5])
s = violinplot(bcF_10);
s.ViolinColor = [0 0.5 0.5];
set(gca,'FontSize',15)

subplot(2,9,17)
ylim([-0.05 1.5])
s = violinplot(bcF0);
s.ViolinColor = [0 0.5 0.5];

subplot(2,9,18)
s = violinplot(full(bcF_JSDM))
s.ViolinColor = [0 0.5 0.5];
ylim([-0.05 1.5])



% 
% subplot(2,6,3)
% ylim([-0.05 1])
% s = violinplot(klB_10);
% s.ViolinColor = [0.5 0 0.5];
% set(gca,'FontSize',15)
% subplot(2,6,4)
% ylim([-0.05 1])
% s = violinplot(full(klB_JSDM))
% s.ViolinColor = [0.5 0 0.5];
% 
% subplot(2,6,5)
% ylim([-0.05 1.5])
% set(gca,'FontSize',15)
% s = violinplot(klF_10);
% s.ViolinColor = [0 0.5 0.5];
% subplot(2,6,6)
% ylim([-0.05 1.5])
% s = violinplot(full(klF_JSDM))
% s.ViolinColor = [0 0.5 0.5];
% 
% 
% subplot(2,6,7)
% ylim([-0.05 1])
% set(gca,'FontSize',15)
% s = violinplot(bcA_10);
% s.ViolinColor = [0.5 0.5 0];
% subplot(2,6,8)
% ylim([-0.05 1])
% s = violinplot(full(bcA_JSDM))
% s.ViolinColor = [0.5 0.5 0];
% 
% subplot(2,6,9)
% ylim([-0.05 1])
% set(gca,'FontSize',15)
% s = violinplot(bcB_10);
% s.ViolinColor = [0.5 0 0.5];
% subplot(2,6,10)
% ylim([-0.05 1])
% s = violinplot(full(bcB_JSDM))
% s.ViolinColor = [0.5 0 0.5];
% 
% subplot(2,6,11)
% ylim([-0.05 1])
% set(gca,'FontSize',15)
% s = violinplot(bcF_10);
% s.ViolinColor = [0 0.5 0.5];
% subplot(2,6,12)
% ylim([-0.05 1])
% s = violinplot(full(bcF_JSDM))
% s.ViolinColor = [0 0.5 0.5];