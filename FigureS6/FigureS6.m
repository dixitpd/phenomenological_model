clc
clear
%

load ../Data/three_kingdoms_cleaned_up
taxonomy = taxonomy(2:end,:);

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
taxb         = taxonomy(goodb,:);

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

%% Filter samples based on large deviations
badsamps         = find(sum(abs(metz > 10)'));
mcow             = metadata_cow;mcow(:,badmets) = [];
mcow(badsamps,:) = [];
xs_a(badsamps,:) = [];
xs_f(badsamps,:) = [];
xs_b(badsamps,:) = [];

load ../Figure2/testing_data
xs_a(testing_data,:) = [];
xs_b(testing_data,:) = [];
xs_f(testing_data,:) = [];
mcow(testing_data,:) = [];

mu_met           = mean(mcow);
st_met           = std(mcow);
metz             = (mcow-mu_met)./st_met;


%% Load model
alph  = 0.01;
nK  = 20;
filen = strcat('../Data/Cows_Training/cows_',num2str(nK),'_',num2str(alph),'.mat');
load(filen)
QA = exp(-[ZC]*thetA);QA = normalize(QA,2,'norm',1);
QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);
QF = exp(-[ZC]*thetF);QF = normalize(QF,2,'norm',1);
M  = ZC*C;
nSamp = size(ZC,1);
nMet = size(M,2);


%% Find clusters

load ../Figure3/gaussian_models
nG = find(a2==min(min(a2)));gmdl = mx{nG};ncmp = gmdl.NumComponents;
props = gmdl.ComponentProportion;
for s=1:nSamp
    for g = 1:ncmp
        mu = gmdl.mu(g,:);
        sg = gmdl.Sigma(:,:,g);
        tmp(g) = mvnpdf(ZC(s,:),mu,sg);
    end
    tmp = tmp.*props;
    id = find(tmp==max(tmp));
    clst(s) = id;

end
%
crctf = 0.4;fd = 0.05;
[crb prb] = corr(metz,xs_b(:,1:end-1),'type','spearman');
t = reshape(prb,46*156,1);t1 = sort(t);t2 = fd*(1:(46*156))/(46*156);

prb_ct = t1(max(find(t2 - t1' > 0)));

for c=1:4
    cl = find(clst==c);
    clst_num(c) = length(cl);
    [t p] = corr(metz(cl,:),xs_b(cl,1:end-1),'type','spearman');
    clr_b(c,:,:) = t;
    plr_b(c,:,:) = p;
    t1 = reshape(p,46*156,1);t2 = fd*(1:(46*156))/(46*156);
    mx = max(find(t2 > t1'));
    prb_cut_cl(c) = t1(mx);
    ttx{c} = p < t1(mx);
end
f0 = mean(prb < prb_ct,2);
g  = 0.25*(ttx{1} + ttx{2} + ttx{3} + ttx{4});gg = mean(g,2);



crb_st = reshape(std(clr_b,1),nMet,156);
crb_mu = reshape(mean(clr_b,1),nMet,156);

colors = [[0.5 0.5 0];[0 0.5 0];[0.8 0 0.5];[0.5 0.5 0.5]];


subplot(1,2,1)
%
xx = reshape(crb,nMet*156,1);
yy = reshape(crb_mu,nMet*156,1);
zz = reshape(crb_st,nMet*156,1);
mdl = fit(xx,yy,'a*x')
hold on

h = histogram2(xx,yy, 'DisplayStyle', 'tile', 'NumBins', 100,'XBinLimits',[-1 1],'YBinLimits',[-1 1]);
h.EdgeColor =  'none';
shading interp
c = gray;
c = flipud(c);c = c*0.75;
colormap(c);


set(gca,'FontSize',15)
set(gca,'FontWeight','bold')


smps = randsample(nSamp,nSamp);
tt = [0 cumsum(clst_num)];
for c=1:4
    ints = (tt(c)+1):tt(c+1);
    cl = smps(ints);
    [t p] = corr(metz(cl,:),xs_b(cl,1:end-1),'type','spearman');
    clr_b(c,:,:) = t;
    plr_b(c,:,:) = p;
    t1 = reshape(p,46*156,1);t2 = fd*(1:(46*156))/(46*156);
    mx = max(find(t2 > t1'));
    prb_cut_cl(c) = t1(mx);
    ttx{c} = p < t1(mx);
end



crb_st = reshape(std(clr_b,1),nMet,156);
crb_mu = reshape(mean(clr_b,1),nMet,156);

subplot(1,2,2)

xx = reshape(crb,nMet*156,1);
yy = reshape(crb_mu,nMet*156,1);
zz = reshape(crb_st,nMet*156,1);
mdl = fit(xx,yy,'a*x')
hold on

h = histogram2(xx,yy, 'DisplayStyle', 'tile', 'NumBins', 100,'XBinLimits',[-1 1],'YBinLimits',[-1 1]);
h.EdgeColor =  'none';
shading interp
c = gray;
c = flipud(c);c = c*0.75;
colormap(c);

set(gca,'FontSize',15)
set(gca,'FontWeight','bold')

