clc
clear
%


%% Chicken
nK =  5;
filen = strcat('Chicken/cecum_',num2str(nK),'.mat')
load(filen)
[nSamp nB] = size(xs_b);
options = statset('MaxIter',5000);
load gaussian_models_chicken
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,10000);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end
t1 = reshape(xs_b(:,1:end-1),nSamp*(nB-1),1);t1(t1==0) = [];
t2 = reshape(QBs(:,1:end-1),10000*(nB-1),1);t2(t2==0) = [];

subplot(3,4,9)
hold on
intrvl = 10.^(-2.5:0.01:-1);
[b a] = hist(t1,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
b = b(1:end-1);a = a(1:end-1);b1 = log10(b);
[b a] = hist(t2,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
b = b(1:end-1);a = a(1:end-1);b2 = log10(b);
x  = log10(a);
mdl = fit(x',b1','m*x+n');m = mdl.m;n = mdl.n;
mdl
intrvl = 10.^(-5:0.01:0);
plot(intrvl,10^(n)*intrvl.^(m),'r')
mdl = fit(x',b2','m*x+n');m = mdl.m;n = mdl.n;
mdl
plot(intrvl,10^(n)*intrvl.^(m),'r--')
[b a] = hist(t1,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
plot(a,b,'k')
[b a] = hist(t2,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
plot(a,b,'k--')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([1e-4 1])
ylim([1e-3 1])


%% Cows
nK =  5;
filen = strcat('Cows/cows_',num2str(nK),'.mat')
load(filen)
[nSamp nB] = size(xs_b);
options = statset('MaxIter',5000);
load gaussian_models_cows
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,10000);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
clear QBs
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end
t1 = reshape(xs_b(:,1:end-1),nSamp*(nB-1),1);t1(t1==0) = [];
t2 = reshape(QBs(:,1:end-1),10000*(nB-1),1);t2(t2==0) = [];

subplot(3,4,10)
hold on
intrvl = 10.^(-2.5:0.01:-1);
[b a] = hist(t1,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
b = b(1:end-1);a = a(1:end-1);b1 = log10(b);
[b a] = hist(t2,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
b = b(1:end-1);a = a(1:end-1);b2 = log10(b);
x  = log10(a);
mdl = fit(x',b1','m*x+n');m = mdl.m;n = mdl.n;
mdl
intrvl = 10.^(-5:0.01:0);
plot(intrvl,10^(n)*intrvl.^(m),'r')
mdl = fit(x',b2','m*x+n');m = mdl.m;n = mdl.n;
mdl
plot(intrvl,10^(n)*intrvl.^(m),'r--')
[b a] = hist(t1,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
plot(a,b,'b')
[b a] = hist(t2,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
plot(a,b,'b--')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([1e-4 1])
ylim([1e-3 1])


%% Human
nK =  5;
filen = strcat('Human/human_',num2str(nK),'.mat')
load(filen)
[nSamp nB] = size(xs_b);
options = statset('MaxIter',5000);
load gaussian_models_human
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,10000);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
clear QBs
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end
t1 = reshape(xs_b(:,1:end-1),nSamp*(nB-1),1);t1(t1==0) = [];
t2 = reshape(QBs(:,1:end-1),10000*(nB-1),1);t2(t2==0) = [];

subplot(3,4,11)
hold on
intrvl = 10.^(-2.5:0.01:-1);
[b a] = hist(t1,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
b = b(1:end-1);a = a(1:end-1);b1 = log10(b);
[b a] = hist(t2,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
b = b(1:end-1);a = a(1:end-1);b2 = log10(b);
x  = log10(a);
mdl = fit(x',b1','m*x+n');m = mdl.m;n = mdl.n;
mdl
intrvl = 10.^(-5:0.01:0);
plot(intrvl,10^(n)*intrvl.^(m),'r')
mdl = fit(x',b2','m*x+n');m = mdl.m;n = mdl.n;
mdl
plot(intrvl,10^(n)*intrvl.^(m),'r--')
[b a] = hist(t1,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
plot(a,b,'m')
[b a] = hist(t2,intrvl);b = b/sum(b);b = 1-cumsum(b);b = abs(b);
plot(a,b,'m--')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([1e-4 1])
ylim([1e-3 1])