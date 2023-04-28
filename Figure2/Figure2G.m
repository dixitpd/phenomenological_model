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
mcow_test = mcow(testing_data,:);
xs_a(testing_data,:) = [];
xs_b(testing_data,:) = [];
xs_f(testing_data,:) = [];
nTest = length(testing_data);
nSamp = size(xs_a,1);

mcow(testing_data,:) = [];

mu_met           = mean(mcow);
st_met           = std(mcow);
metz             = (mcow-mu_met)./st_met;

metz_test        = (mcow_test-mu_met)./st_met;


%% Load
alph  = 0.01;
nK  = 20;
filen = strcat('../Data/Cows_Training/cows_',num2str(nK),'_',num2str(alph),'.mat')
load(filen)

Zr = ZC;thetAr = thetA;thetBr = thetB;thetFr = thetF;Cr = C;
clear coeffs intr

dt      = full([xs_a]);
dt_test = full([xs_a_test]);
for k = 1:nK
    [mdl c]    = lasso(dt,Zr(:,k),'lambda',0.0);
    coeffs(k,:) = mdl;
    intr(k) = c.Intercept;
end
Zpd1 = dt*coeffs' + intr;
Zp   = dt_test*coeffs' + intr;
MA    = Zp*Cr;
[ca pa] = corr(MA,metz_test);ca = diag(ca);pa = diag(pa);



clear coeffs intr
dt      = full([xs_b]);
dt_test = full([xs_b_test]);
for k = 1:nK
    [mdl c]    = lasso(dt,Zr(:,k),'lambda',0.0);
    coeffs(k,:) = mdl;
    intr(k) = c.Intercept;
end
Zpd1 = dt*coeffs' + intr;
Zp   = dt_test*coeffs' + intr;
MB    = Zp*Cr;
[cb pb] = corr(MB,metz_test);cb = diag(cb);pb = diag(pb);


clear coeffs intr
dt      = full([xs_f]);
dt_test = full([xs_f_test]);
for k = 1:nK
    [mdl c]    = lasso(dt,Zr(:,k),'lambda',0.0);
    coeffs(k,:) = mdl;
    intr(k) = c.Intercept;
end
Zpd1 = dt*coeffs' + intr;
Zp   = dt_test*coeffs' + intr;
MF    = Zp*Cr;
[cf pf] = corr(MF,metz_test);cf = diag(cf);pf = diag(pf);



clear coeffs intr
dt      = full([xs_a xs_b xs_f]);
dt_test = full([xs_a_test xs_b_test xs_f_test]);
for k = 1:nK
    [mdl c]    = lasso(dt,Zr(:,k),'lambda',0.0);
    coeffs(k,:) = mdl;
    intr(k) = c.Intercept;
end
Zpd1 = dt*coeffs' + intr;
Zp   = dt_test*coeffs' + intr;
MAll    = Zp*Cr;
[cll pall] = corr(MAll,metz_test);cll = diag(cll);pall = diag(pall);
[a b] = sort(cll,'descend');

hold on
intr = 1:(46);
bar(intr,cll(b),'FaceColor',[0.9 0.9 0.9],'BarWidth',0.2)
bar(intr+0.2,ca(b),'FaceColor',[0.5 0.5 0],'BarWidth',0.2)
bar(intr+0.4,cb(b),'FaceColor',[0.5 0.0 0.5],'BarWidth',0.2)
bar(intr+0.6,cf(b),'FaceColor',[0.0 0.5 0.5],'BarWidth',0.2)
xticks(intr)
xticklabels(metnames(b))
xtickangle(90)
ytickangle(90)
set(gca,'FontSize',15)
ylim([-0.1 1])


clear coeffs intr
dt = full(1.0*(xs_a>0));
dt_test = full(1.0*(xs_a_test>0));
for m = 1:nMet
    [mdl c]    = lasso(dt,metz(:,m),'lambda',0.0);
    coeffs(m,:) = mdl;
    intr(m) = c.Intercept;
end
MAD = dt_test*coeffs' + intr;
[cad pad] = corr(MAD,metz_test);cad = diag(cad);pad = diag(pad);

clear coeffs intr
dt = full(1.0*(xs_b > 0));
dt_test = full(1.0*(xs_b_test > 0));
for m = 1:nMet
    [mdl c]    = lasso(dt,metz(:,m),'lambda',0.0);
    coeffs(m,:) = mdl;
    intr(m) = c.Intercept;
end
MBD = dt_test*coeffs' + intr;
[cbd pbd] = corr(MBD,metz_test);cbd = diag(cbd);pad = diag(pbd);


clear coeffs intr
dt = full(1.0*(xs_f > 0));
dt_test = full(1.0*(xs_f_test>0));
for m = 1:nMet
    [mdl c]    = lasso(dt,metz(:,m),'lambda',0.0);
    coeffs(m,:) = mdl;
    intr(m) = c.Intercept;
end
MFD = dt_test*coeffs' + intr;
[cfd pfd] = corr(MFD,metz_test);cfd = diag(cfd);pad = diag(pfd);

% 


[ca cad cb cbd cf cfd]
