clc
clear
%

%% Chicken
nK =  5;
filen = strcat('Chicken/cecum_',num2str(nK),'.mat')
load(filen)
nSamp = size(ZC,1);
load gaussian_models_chicken
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,10000);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end
% 
subplot(3,4,2)
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,-5],'Upper',[3,5],'StartPoint',[1 1]);
hold on
c1 = cov(xs_b);c2 = cov(QBs);[a b] = size(c1);
c1 = reshape(c1,a*b,1);c2 = reshape(c2,a*b,1);
scatter(abs(c1),abs(c2),20,'k','filled')
intrvl = 10.^(-10:0.1:0);
plot(intrvl,intrvl,'k--')
set(gca,'xscale','log')
set(gca,'yscale','log')
'Correlations Chicken'
corr(c1,c2,'type','spearman')
%% Cows
nK =  5;
filen = strcat('Cows/cows_',num2str(nK),'.mat')
load(filen)
nSamp = size(ZC,1);
load gaussian_models_cows
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,10000);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
clear QBs;
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end
% 
subplot(3,4,3)
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,-5],'Upper',[3,5],'StartPoint',[1 1]);
hold on
c1 = cov(xs_b);c2 = cov(QBs);[a b] = size(c1);
c1 = reshape(c1,a*b,1);c2 = reshape(c2,a*b,1);
scatter(abs(c1),abs(c2),20,'b','filled')
intrvl = 10.^(-10:0.1:0);
plot(intrvl,intrvl,'k--')
set(gca,'xscale','log')
set(gca,'yscale','log')
'Correlations Cows'
corr(c1,c2,'type','spearman')% 
% 
%% Human
nK =  5;
filen = strcat('Human/human_',num2str(nK),'.mat')
load(filen)
nSamp = size(ZC,1);
load gaussian_models_human
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,10000);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
clear QBs;
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end
% 
subplot(3,4,4)
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,-5],'Upper',[3,5],'StartPoint',[1 1]);
hold on
c1 = cov(xs_b);c2 = cov(QBs);[a b] = size(c1);
c1 = reshape(c1,a*b,1);c2 = reshape(c2,a*b,1);
scatter(abs(c1),abs(c2),20,'m','filled')
intrvl = 10.^(-10:0.1:0);
plot(intrvl,intrvl,'m--')
set(gca,'xscale','log')
set(gca,'yscale','log')
'Correlations Human'
corr(c1,c2,'type','spearman')