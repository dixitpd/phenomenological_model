clc
clear
%
% 

% Cows
filen = strcat('Cows/cows_1.mat')
load(filen)
nSamp = size(xs_b,1);

klCow = [];
for s = 1:nSamp
    f = JSD(xs_b(s,:),xs_b);f(s) = [];klCow = [klCow f];
end
g = f_braycurtis(xs_b');g = triu(g);g(g==0) = [];bc_cow = g;


nCMax = 10;
nSamp = size(xs_b,1);
for nC = 1:nCMax
    filen = strcat('Cows/cows_',num2str(nC),'.mat')
    load(filen)
    Q = exp(-ZC*thetB);Q = normalize(Q,2,'norm',1);
    kl_cow_model(nC,:) = JSD(Q,xs_b);

    for s=1:nSamp
        bc_cow_model(nC,s) = bc_pairs(xs_b(s,:),Q(s,:));
    end
end

clear xs_b
filen = strcat('Chicken/cecum_1.mat')
load(filen)
nSamp = size(xs_b,1);
klChicken = [];
for s = 1:nSamp
    f = JSD(xs_b(s,:),xs_b);f(s) = [];klChicken = [klChicken f];
end
g = f_braycurtis(xs_b');g = triu(g);g(g==0) = [];bc_chicken = g;

for nC = 1:nCMax
    filen = strcat('Chicken/cecum_',num2str(nC),'.mat')
    load(filen)
    Q = exp(-ZC*thetB);Q = normalize(Q,2,'norm',1);
    kl_chicken_model(nC,:) = JSD(Q,xs_b);

    for s=1:nSamp
        bc_chicken_model(nC,s) = bc_pairs(xs_b(s,:),Q(s,:));
    end
end

clear xs_b
filen = strcat('Human/human_1.mat')
load(filen)
nSamp = size(xs_b,1);
klHuman = [];
for s = 1:nSamp
    f = JSD(xs_b(s,:),xs_b);f(s) = [];klHuman = [klHuman f];
end
g = f_braycurtis(xs_b');g = triu(g);g(g==0) = [];bc_human = g;
for nC = 1:nCMax
    filen = strcat('Human/human_',num2str(nC),'.mat')
    load(filen)
    Q = exp(-ZC*thetB);Q = normalize(Q,2,'norm',1);
    kl_human_model(nC,:) = JSD(Q,xs_b);

    for s=1:nSamp
        bc_human_model(nC,s) = bc_pairs(xs_b(s,:),Q(s,:));
    end
end
% 
% 

prct = 1;
ints = 1:nCMax;
subplot(3,4,1)
hold on
plot(ints,ones(nCMax,1)*prctile(klCow,prct),'b--')
plot(ints,mean(kl_cow_model'),'b')
%
plot(ints,ones(nCMax,1)*prctile(klChicken,prct),'k--')
plot(ints,mean(kl_chicken_model'),'k')
%
plot(ints,ones(nCMax,1)*prctile(klHuman,prct),'m--')
plot(ints,mean(kl_human_model'),'m')
ylim([0 0.5])

subplot(3,4,5)
hold on
plot(ints,ones(nCMax,1)*prctile(bc_cow,prct),'b--')
plot(ints,mean(bc_cow_model'),'b')
%
plot(ints,ones(nCMax,1)*prctile(bc_chicken,prct),'k--')
plot(ints,mean(bc_chicken_model'),'k')
%
plot(ints,ones(nCMax,1)*prctile(bc_human,prct),'m--')
plot(ints,mean(bc_human_model'),'m')
ylim([0 0.5])
