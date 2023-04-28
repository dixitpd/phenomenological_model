clc
clear

nK = 16;alph = 0.001; prcs = 1;

filen = strcat('../Data/Cows_Training/cows_',num2str(nK),'_',num2str(alph),'.mat');
load(filen)
[nSamp] = size(xs_a,1);
xs_a(xs_a < 0) = 0;xs_f(xs_f < 0) = 0;

fa = [];fb = [];ff = [];
for s = 1:nSamp
    f = JSD(xs_a(s,:),xs_a);f(s) = [];fa = [fa f];
    f = JSD(xs_b(s,:),xs_b);f(s) = [];fb = [fb f];
    f = JSD(xs_f(s,:),xs_f);f(s) = [];ff = [ff f];
end
klA_pr = prctile(fa,prcs);
klB_pr = prctile(fb,prcs);
klF_pr = prctile(ff,prcs);

% 
% % Bray curtis to closest data point
g = f_braycurtis(xs_a');g = triu(g);g(g==0) = [];bcA_pr = prctile(g,prcs);
g = f_braycurtis(xs_b');g = triu(g);g(g==0) = [];bcB_pr = prctile(g,prcs);
g = f_braycurtis(xs_f');g = triu(g);g(g==0) = [];bcF_pr = prctile(g,prcs);


nn = 4:4:28;
aa = [0.0001 0.001 0.01 0.02 0.05];

for k=1:length(nn)
    nK = nn(k);
    for ax = 1:length(aa)
        alph = aa(ax);
        [nK alph]
        kxx(k,ax) = nK;
        axx(k,ax) = alph;

        filen = strcat('../Data/Cows_Training/cows_',num2str(nK),'_',num2str(alph),'.mat');
        load(filen)
        QA = exp(-[ZC]*thetA);QA = normalize(QA,2,'norm',1);
        QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);
        QF = exp(-[ZC]*thetF);QF = normalize(QF,2,'norm',1);
        M  = ZC*C;
        xs_a(xs_a < 0) = 0;xs_f(xs_f < 0) = 0;

        t1 = sum(sum(xs_a.*log(QA)));
        t2 = sum(sum(xs_b.*log(QB)));
        t3 = sum(sum(xs_f.*log(QF)));
        t4 = (M-metz);t4 = sum(sum(t4.*t4));
        cst_q(k,ax) = -(1-alph)*(t1+t2+t3);
        cst_m(k,ax) = alph*t4;

        % KL divergence with the model
        klA = JSD(QA,xs_a);
        klB = JSD(QB,xs_b);
        klF = JSD(QF,xs_f);

        mcrM(k,ax)  = mean(diag(corr(M,metz)));
        micrM(k,ax) = min(diag(corr(M,metz)));
        mxcrM(k,ax) = max(diag(corr(M,metz)));

        mklA(k,ax) = mean(klA);
        mklB(k,ax) = mean(klB);
        mklF(k,ax) = mean(klF);
        fklA(k,ax) = mean(klA < klA_pr);
        fklB(k,ax) = mean(klB < klB_pr);
        fklF(k,ax) = mean(klF < klF_pr);


        % Bray curtis with the model
        for s=1:nSamp
            e1 = xs_a(s,:);e2 = QA(s,:);
            f  =  f_braycurtis([e1;e2]');f = f(1,2);
            bcA(s) = f;

            e1 = xs_b(s,:);e2 = QB(s,:);
            f  =  f_braycurtis([e1;e2]');f = f(1,2);
            bcB(s) = f;

            e1 = xs_f(s,:);e2 = QF(s,:);
            f  =  f_braycurtis([e1;e2]');f = f(1,2);
            bcF(s) = f;
        end

        mbcA(k,ax) = mean(bcA);
        mbcB(k,ax) = mean(bcB);
        mbcF(k,ax) = mean(bcF);

        fbcA(k,ax) = mean(bcA < bcA_pr);
        fbcB(k,ax) = mean(bcB < bcA_pr);
        fbcF(k,ax) = mean(bcF < bcF_pr);

    end
end

color = [[0.5 0.5 0];[0.5 0 0.5];[0 0.5 0.5]];


subplot(2,5,1)
hold on
plot(nn,mklA,'o-','color',color(1,:))
plot(0:30,ones(31,1)*mean(klA_pr),'--','color',color(1,:))
ylim([0 0.2])
subplot(2,5,6)
hold on
plot(nn,mbcA,'o-','color',color(1,:))
plot(0:30,ones(31,1)*mean(bcA_pr),'--','color',color(1,:))
ylim([0 0.25])

subplot(2,5,2)
hold on
plot(nn,mklB,'o-','color',color(2,:))
plot(0:30,ones(31,1)*mean(klB_pr),'--','color',color(2,:))
ylim([0 0.2])
subplot(2,5,7)
hold on
plot(nn,mbcB,'o-','color',color(2,:))
plot(0:30,ones(31,1)*mean(bcB_pr),'--','color',color(2,:))
ylim([0 0.25])

subplot(2,5,3)
hold on
plot(nn,mklF,'o-','color',color(3,:))
plot(0:30,ones(31,1)*mean(klF_pr),'--','color',color(3,:))
ylim([0 0.2])
subplot(2,5,8)
hold on
plot(nn,mbcB,'o-','color',color(3,:))
plot(0:30,ones(31,1)*mean(bcF_pr),'--','color',color(3,:))
ylim([0 0.25])

subplot(2,5,4)
hold on
plot(nn,fklA(:,3),'o-','color',color(1,:))
plot(nn,fklB(:,3),'o-','color',color(2,:))
plot(nn,fklF(:,3),'o-','color',color(3,:))
ylim([0 1])

plot([nn ],0.4*ones(1,length(nn)),'r--')
plot([nn ],0.5*ones(1,length(nn)),'r--')
plot([nn ],0.6*ones(1,length(nn)),'r--')

subplot(2,5,9)
hold on
plot(nn,fbcA(:,3),'o-','color',color(1,:))
plot(nn,fbcB(:,3),'o-','color',color(2,:))
plot(nn,fbcF(:,3),'o-','color',color(3,:))
ylim([0 1])

plot([nn ],0.4*ones(1,length(nn)),'r--')
plot([nn ],0.5*ones(1,length(nn)),'r--')
plot([nn ],0.6*ones(1,length(nn)),'r--')


t2 = mklA - (klA_pr)
t3 = mklB - (klB_pr)
t4 = mklF - (klF_pr)
t5 = mbcA - (bcA_pr)
t6 = mbcB - (bcB_pr)
t7 = mbcF - (bcF_pr)
% 


