clc
clear
%
% Chicken data:  1.63   0.21
% Chicken model: 2.03   0.32
% Cow data:      1.37   0.16
% Cow model:     1.62   0.23
% Human data:    1.61   0.22
% Human model:   1.98   0.32

epl = 0.99;
%% Chicken
nK =  5;
filen = strcat('Chicken/cecum_',num2str(nK),'.mat');
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
subplot(3,4,6)
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,-5],'Upper',[3,5],'StartPoint',[1 1]);
hold on
intrvl = (-10:0.1:0);m = mean(xs_b);m(end) = [];v = var(xs_b);v(end) = [];
m = log10(m);v = log10(v);mdl = fit(m',v','a*x+b',fo);a = mdl.a;b = mdl.b;
plot(10.^(intrvl),10.^(a*intrvl+b),'k')
plot(10.^(m),10.^(v),'ko')
t = confint(mdl,epl);
'Chicken data'
[a abs(t(:,1)' - a)]

m = mean(QBs);m(end) = [];v = var(QBs);v(end) = [];
m = log10(m);v = log10(v);mdl = fit(m',v','a*x+b',fo);a = mdl.a;b = mdl.b;
plot(10.^(intrvl),10.^(a*intrvl+b),'k--')
plot(10.^(m),10.^(v),'k*')
t = confint(mdl,epl);
'Chicken model'
[a abs(t(:,1)' - a)]
set(gca,'xscale','log')
set(gca,'yscale','log')

%% Cows
nK =  5;
filen = strcat('Cows/cows_',num2str(nK),'.mat');
load(filen)
nSamp = size(ZC,1);
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
% 
subplot(3,4,7)
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,-5],'Upper',[3,5],'StartPoint',[1 1]);
hold on
intrvl = (-10:0.1:0);m = mean(xs_b);m(end) = [];v = var(xs_b);v(end) = [];
m = log10(m);v = log10(v);mdl = fit(m',v','a*x+b',fo);a = mdl.a;b = mdl.b;
plot(10.^(intrvl),10.^(a*intrvl+b),'b')
plot(10.^(m),10.^(v),'bo')
t = confint(mdl,epl);
'Cow data'
[a abs(t(:,1)' - a)]

m = mean(QBs);m(end) = [];v = var(QBs);v(end) = [];
m = log10(m);v = log10(v);mdl = fit(m',v','a*x+b',fo);a = mdl.a;b = mdl.b;
plot(10.^(intrvl),10.^(a*intrvl+b),'b--')
plot(10.^(m),10.^(v),'b*')
t = confint(mdl,epl);
'Cow model'
[a abs(t(:,1)' - a)]
set(gca,'xscale','log')
set(gca,'yscale','log')

%% Human
nK =  5;
filen = strcat('Human/human_',num2str(nK),'.mat');
load(filen)
nSamp = size(ZC,1);
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
% 
subplot(3,4,8)
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,-5],'Upper',[3,5],'StartPoint',[1 1]);
hold on
intrvl = (-10:0.1:0);m = mean(xs_b);m(end) = [];v = var(xs_b);v(end) = [];
m = log10(m);v = log10(v);mdl = fit(m',v','a*x+b',fo);a = mdl.a;b = mdl.b;
plot(10.^(intrvl),10.^(a*intrvl+b),'m')
plot(10.^(m),10.^(v),'mo')
t = confint(mdl,epl);
'Human data'
[a abs(t(:,1)' - a)]

m = mean(QBs);m(end) = [];v = var(QBs);v(end) = [];
m = log10(m);v = log10(v);mdl = fit(m',v','a*x+b',fo);a = mdl.a;b = mdl.b;
plot(10.^(intrvl),10.^(a*intrvl+b),'m--')
plot(10.^(m),10.^(v),'m*')
t = confint(mdl,epl);
'Human model'
[a abs(t(:,1)' - a)]
set(gca,'xscale','log')
set(gca,'yscale','log')