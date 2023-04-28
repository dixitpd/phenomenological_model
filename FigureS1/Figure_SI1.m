clc
clear


%% Chicken
clear
nK =  5;
filen = strcat('../Figure1/Chicken/cecum_',num2str(nK),'.mat')
load(filen)
[nSamp nB] = size(xs_b);
load ../Figure1/gaussian_models_chicken
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,nSamp);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end

g = f_braycurtis([QBs;xs_b]');g = g(1:nSamp,nSamp+1:end);bc_chicken_model = min(g');
g = f_braycurtis(xs_b');g = g + 1e10*eye(nSamp);bc_chicken_neighbor = min(g);

for s=1:nSamp
    f = JSD(xs_b(s,:),xs_b);f(s) = [];
    kl_chicken_neighbor(s) = min(f);
    f = JSD(QBs(s,:),xs_b);
    kl_chicken_model(s) = min(f);
end

subplot(2,3,1)
hold on
[b a] = hist(kl_chicken_neighbor,0:0.025:1);b = b/sum(b);
plot(a,b,'k','linewidth',1.5)
[b a] = hist(kl_chicken_model,0:0.025:1);b = b/sum(b);
plot(a,b,'k--','linewidth',1.5)
xlim([0 0.6])

subplot(2,3,4)
hold on
[b a] = hist(bc_chicken_neighbor,0:0.025:1);b = b/sum(b);
plot(a,b,'k','linewidth',1.5)
[b a] = hist(bc_chicken_model,0:0.025:1);b = b/sum(b);
plot(a,b,'k--','linewidth',1.5)
xlim([0 0.6])


%% Cows
clear
nK =  5;
filen = strcat('../Figure1/Cows/cows_',num2str(nK),'.mat')
load(filen)
[nSamp nB] = size(xs_b);
load ../Figure1/gaussian_models_cows
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,nSamp);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
clear QBs
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end

g = f_braycurtis([QBs;xs_b]');g = g(1:nSamp,nSamp+1:end);bc_cow_model = min(g);
g = f_braycurtis(xs_b');g = g + 1e10*eye(nSamp);bc_cow_neighbor = min(g);

for s=1:nSamp
    f = JSD(xs_b(s,:),xs_b);f(s) = [];
    kl_cow_neighbor(s) = min(f);
    f = JSD(QBs(s,:),xs_b);
    kl_cow_model(s) = min(f);
end

subplot(2,3,2)
hold on
[b a] = hist(kl_cow_neighbor,0:0.025:1);b = b/sum(b);
plot(a,b,'b','linewidth',1.5)
[b a] = hist(kl_cow_model,0:0.025:1);b = b/sum(b);
plot(a,b,'b--','linewidth',1.5)
xlim([0 0.4])
subplot(2,3,5)
hold on
[b a] = hist(bc_cow_neighbor,0:0.025:1);b = b/sum(b);
plot(a,b,'b','linewidth',1.5)
[b a] = hist(bc_cow_model,0:0.025:1);b = b/sum(b);
plot(a,b,'b--','linewidth',1.5)
xlim([0 0.4])

%% Human
clear
nK =  5;
filen = strcat('../Figure1/Human/human_',num2str(nK),'.mat')
load(filen)
[nSamp nB] = size(xs_b);
load ../Figure1/gaussian_models_human
nG = find(a2==min(min(a2)));
gmdl = mx{nG};
Zsamp = random(gmdl,nSamp);
qq   = exp(-Zsamp*thetB);qq = normalize(qq,2,'norm',1);
clear QBs
for iter=1:size(Zsamp,1);
    x = mnrnd(5000,qq(iter,:));x = x/sum(x);
    QBs(iter,:) = x;
end

g = f_braycurtis([QBs;xs_b]');g = g(1:nSamp,nSamp+1:end);bc_human_model = min(g);
g = f_braycurtis(xs_b');g = g + 1e10*eye(nSamp);bc_human_neighbor = min(g);

for s=1:nSamp
    f = JSD(xs_b(s,:),xs_b);f(s) = [];
    kl_human_neighbor(s) = min(f);
    f = JSD(QBs(s,:),xs_b);
    kl_human_model(s) = min(f);
end

subplot(2,3,3)
hold on
[b a] = hist(kl_human_neighbor,0:0.05:1);b = b/sum(b);
plot(a,b,'m','linewidth',1.5)
[b a] = hist(kl_human_model,0:0.05:1);b = b/sum(b);
plot(a,b,'m--','linewidth',1.5)
xlim([0 0.6])

subplot(2,3,6)
hold on
[b a] = hist(bc_human_neighbor,0:0.05:1);b = b/sum(b);
plot(a,b,'m','linewidth',1.5)
[b a] = hist(bc_human_model,0:0.05:1);b = b/sum(b);
plot(a,b,'m--','linewidth',1.5)
xlim([0 0.6])