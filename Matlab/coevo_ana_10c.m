function [x cA cB cC w maxGen pA pB pC]=coevo_ana_08c(...
        xInit,...                   % initial conditions
        tmax,...                    % periods per generation
        gen,...                     % generations
        selCoef,...                 % maximum selection coefficient
        nStrats,...                 % total possible number of strategies
        vertTransm,...              % whether there is vertical transmission
        indLearn,...                % indiv. learing perf.; 0-> win-stay lose-shift
        mutationRate,...            % probability of mutation occurring
        regime,...                  % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
        pA0,...
        pB0,...
        pC0,...                     % initial values of pA, pB, pC
        incr,...                    % increment at which the environment becomes better or worse
        pincr,...                   % probability that environmental quality changes at all after each period
        dpA,...                     % shift of pA, pB, pC
        dpB,...
        dpC);


% INITIALIZATION
x=zeros(nStrats,gen);       % frequency matrix
cA=zeros(nStrats,tmax);     % percent choices for A
cB=zeros(nStrats,tmax);     % percent choices for B
cC=zeros(nStrats,tmax);     % percent choices for C
w=zeros(nStrats,gen);       % fitness matrix

x(:,1)=xInit;

extinctVec=zeros(nStrats,1);

for g=1:gen
    
    % ENVIRONMENT
    % routine to determine pA, pB, and pC
    [pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);
    [pC,pD] = randomenvironment4(tmax,regime,incr,pincr,pC0,pC0);
    % possible shift
    pA=min(1-incr,max(incr,pA+dpA));pB=min(1-incr,max(pB+dpB,incr));pC=min(1-incr,max(pC+dpC,incr));
    % we define C is the better and F as the worse choice
    % remember last pA and pB
    pA0=pA(tmax);pB0=pB(tmax);pC0=pC(tmax);
    
    
    % FIRST PERIOD'S CHOICE OF THE STRATEGIES
    if g==1||vertTransm==0     % in very first generation
        cA(:,1)=1/3;
        cB(:,1)=1/3;
        cC(:,1)=1/3;
    elseif vertTransm==1
        cA(:,1)=cAlast;
        cB(:,1)=cBlast;
        cC(:,1)=cClast;
    end
    
    
    % CHOICES OF THE STRATEGIES
    for t=2:tmax
        
        cAmean=sum(x(:,g).*cA(:,t-1));
        cBmean=sum(x(:,g).*cB(:,t-1));
        cCmean=sum(x(:,g).*cC(:,t-1));
        
        % individual learners
        if indLearn==0
            cAind=(1-pB(t))*(1-pC(t))/(3+pB(t)*(pC(t)-2)-2*pC(t)+pA(t)*(pB(t)+pC(t)-2));
            cBind=(1-pA(t))*(1-pC(t))/(3+pB(t)*(pC(t)-2)-2*pC(t)+pA(t)*(pB(t)+pC(t)-2));
            cCind=1-cAind-cBind;
        else
            cind=((1-2*(1-indLearn))*sign(pC(t)-pF(t))+1)/2;
            'error: not defined yet'
            break
        end
        cA(1,t)=cAind;
        cB(1,t)=cBind;
        cC(1,t)=cCind;
        

        % 02 conformists
        cA(2,t)=cAmean.*(-2.*cAmean.^2+cAmean.*(3-2.*cBmean)-2.*(-1+cBmean).*cBmean);
        cB(2,t)=cBmean.*(3.*cBmean-2.*(cAmean.^2+cAmean.*(-1+cBmean)+cBmean.^2));
        cC(2,t)=1-cA(2,t)-cB(2,t);

        % 03 contrarian
        cA(3,t)=(-(-1+cAmean).^3-(-1+cAmean).*(3+cAmean).*cBmean-(3+cAmean).*cBmean.^2)/2;
        cB(3,t)=(-(-1+cBmean).^3-cAmean.^2.*(3+cBmean)-cAmean.*(-1+cBmean).*(3+cBmean))/2;
        cC(3,t)=1-cA(3,t)-cB(3,t);

        %  4 PBSL [1/-1]    equal weights
        cA(4,t)=(1-3.*cAmean+3.*cAmean.^2-cAmean.^3+3.*cBmean-2.*cAmean.*cBmean-cAmean.^2.*cBmean-3.*cBmean.^2-cAmean.*cBmean.^2+6.*cAmean.*pA(t)-6.*cAmean.^2.*pA(t)-4.*cAmean.*cBmean.*pA(t)+4.*cAmean.^2.*cBmean.*pA(t)+4.*cAmean.*cBmean.^2.*pA(t)+6.*cAmean.^3.*pA(t).^2-4.*cAmean.^3.*pA(t).^3-6.*cBmean.*pB(t)+8.*cAmean.*cBmean.*pB(t)-2.*cAmean.^2.*cBmean.*pB(t)+6.*cBmean.^2.*pB(t)-2.*cAmean.*cBmean.^2.*pB(t)-2.*cAmean.*cBmean.*pA(t).*pB(t)-4.*cAmean.^2.*cBmean.*pA(t).*pB(t)+2.*cAmean.*cBmean.^2.*pA(t).*pB(t)+6.*cAmean.^2.*cBmean.*pA(t).^2.*pB(t)-3.*cBmean.^3.*pB(t).^2-6.*cAmean.*cBmean.^2.*pA(t).*pB(t).^2+2.*cBmean.^3.*pB(t).^3+2.*(-1+cAmean+cBmean).*(-3.*cAmean.^2.*(-1+pA(t)).*pA(t)+cBmean.*(3-3.*pB(t)+cAmean.*(-1+pA(t)+pB(t)))).*pC(t)-3.*(-1+cAmean+cBmean).^2.*(1-cAmean-cBmean+2.*cAmean.*pA(t)).*pC(t).^2-2.*(-1+cAmean+cBmean).^3.*pC(t).^3)/2;
        cB(4,t)=(1+3.*cAmean-3.*cAmean.^2-3.*cBmean-2.*cAmean.*cBmean-cAmean.^2.*cBmean+3.*cBmean.^2-cAmean.*cBmean.^2-cBmean.^3-6.*cAmean.*pA(t)+6.*cAmean.^2.*pA(t)+8.*cAmean.*cBmean.*pA(t)-2.*cAmean.^2.*cBmean.*pA(t)-2.*cAmean.*cBmean.^2.*pA(t)-3.*cAmean.^3.*pA(t).^2+2.*cAmean.^3.*pA(t).^3+6.*cBmean.*pB(t)-4.*cAmean.*cBmean.*pB(t)+4.*cAmean.^2.*cBmean.*pB(t)-6.*cBmean.^2.*pB(t)+4.*cAmean.*cBmean.^2.*pB(t)-2.*cAmean.*cBmean.*pA(t).*pB(t)+2.*cAmean.^2.*cBmean.*pA(t).*pB(t)-4.*cAmean.*cBmean.^2.*pA(t).*pB(t)-6.*cAmean.^2.*cBmean.*pA(t).^2.*pB(t)+6.*cBmean.^3.*pB(t).^2+6.*cAmean.*cBmean.^2.*pA(t).*pB(t).^2-4.*cBmean.^3.*pB(t).^3+2.*(-1+cAmean+cBmean).*(cAmean.*(-3+cBmean).*(-1+pA(t))+cAmean.*cBmean.*pB(t)-3.*cBmean.^2.*(-1+pB(t)).*pB(t)).*pC(t)+3.*(-1+cAmean+cBmean).^2.*(-1+cAmean+cBmean-2.*cBmean.*pB(t)).*pC(t).^2-2.*(-1+cAmean+cBmean).^3.*pC(t).^3)/2;
        cC(4,t)=1-cA(4,t)-cB(4,t);
        
        %  5 PBSL [3/-1]    gains-biased
        cA(5,t)=(cAmean.^3.*(2.*pA(t).^3+(-1+pC(t)).^3-12.*pA(t).^2.*pC(t)-6.*pA(t).*(-2+pC(t)).*pC(t))-(-1+cBmean.*(pB(t)-pC(t))+pC(t)).*(cBmean.*(-3+5.*pB(t)-2.*pC(t)).*(-1+pC(t))+(-1+pC(t)).^2+cBmean.^2.*(-3+pB(t).^2+pB(t).*(3-5.*pC(t))+pC(t).*(3+pC(t))))+cAmean.^2.*(3.*(-(-1+pC(t)).^3+4.*pA(t).*(-1+pC(t)).*pC(t)+pA(t).^2.*(-2+4.*pC(t)))+cBmean.*(12.*pA(t).^2.*(pB(t)-pC(t))+2.*pA(t).*(2-5.*pB(t)+(7-6.*pC(t)).*pC(t))+(-1+pC(t)).*(1+2.*pB(t)-6.*pB(t).*pC(t)+3.*pC(t).^2)))+cAmean.*(3.*(-1+pC(t)).*((-1+pC(t)).^2-2.*pA(t).*(1+pC(t)))+2.*cBmean.*(-(pA(t).*(2+pB(t)+pC(t)-6.*pC(t).^2))+(-1+pC(t)).*(1-4.*pB(t)+6.*pB(t).*pC(t)-3.*pC(t).^2))+cBmean.^2.*(2.*pA(t).*(2+pB(t)-3.*pB(t).^2+pC(t)-3.*pC(t).^2)+(-1+pC(t)).*(1+6.*pB(t).^2+3.*pC(t).*(2+pC(t))-4.*pB(t).*(1+3.*pC(t))))))/2;
        cB(5,t)=(cAmean.^3.*(-1+pA(t)).^3)/2+3.*cAmean.^2.*cBmean.*(-1+pA(t)).^2.*pB(t)-6.*cAmean.^2.*cBmean.*(-1+pA(t)).*pA(t).*pB(t)+6.*cAmean.*cBmean.^2.*(-1+pA(t)).*(-1+pB(t)).*pB(t)+3.*cBmean.^3.*(-1+pB(t)).^2.*pB(t)-3.*cAmean.*cBmean.^2.*(-1+pA(t)).*pB(t).^2+3.*cAmean.*cBmean.^2.*pA(t).*pB(t).^2-3.*cBmean.^3.*(-1+pB(t)).*pB(t).^2+cBmean.^3.*pB(t).^3+3.*cAmean.^2.*(-1+cAmean+cBmean).*(-1+pA(t)).^2.*(-1+pC(t))+2.*cAmean.*cBmean.*(-1+cAmean+cBmean).*(-1+pA(t)).*(-1+pB(t)).*(-1+pC(t))-6.*cAmean.*cBmean.*(-1+cAmean+cBmean).*(-1+pA(t)).*pB(t).*(-1+pC(t))+3.*cAmean.*cBmean.*(-1+cAmean+cBmean).*pA(t).*pB(t).*(-1+pC(t))-6.*cBmean.^2.*(-1+cAmean+cBmean).*(-1+pB(t)).*pB(t).*(-1+pC(t))+3.*cBmean.^2.*(-1+cAmean+cBmean).*pB(t).^2.*(-1+pC(t))-3.*cAmean.*(-1+cAmean+cBmean).^2.*(-1+pA(t)).*(-1+pC(t)).^2+3.*cBmean.*(-1+cAmean+cBmean).^2.*pB(t).*(-1+pC(t)).^2+((-1+cAmean+cBmean).^3.*(-1+pC(t)).^3)/2+3.*cAmean.*cBmean.*(-1+cAmean+cBmean).*(-1+pA(t)).*pB(t).*pC(t)-2.*cAmean.*cBmean.*(-1+cAmean+cBmean).*pA(t).*pB(t).*pC(t)-3.*cBmean.^2.*(-1+cAmean+cBmean).*pB(t).^2.*pC(t)-6.*cBmean.*(-1+cAmean+cBmean).^2.*pB(t).*(-1+pC(t)).*pC(t);
        cC(5,t)=1-cA(5,t)-cB(5,t);
        
        %  6 PBSL [1/-3]    loss-averse
        cA(6,t)=(1-3.*cAmean+3.*cAmean.^2-cAmean.^3+3.*cBmean-2.*cAmean.*cBmean-cAmean.^2.*cBmean-3.*cBmean.^2-cAmean.*cBmean.^2+6.*cAmean.*pA(t)-12.*cAmean.^2.*pA(t)+6.*cAmean.^3.*pA(t)-4.*cAmean.*cBmean.*pA(t)+4.*cAmean.^2.*cBmean.*pA(t)+4.*cAmean.*cBmean.^2.*pA(t)+6.*cAmean.^2.*pA(t).^2-6.*cAmean.^3.*pA(t).^2+2.*cAmean.^3.*pA(t).^3-6.*cBmean.*pB(t)+8.*cAmean.*cBmean.*pB(t)-2.*cAmean.^2.*cBmean.*pB(t)+12.*cBmean.^2.*pB(t)-8.*cAmean.*cBmean.^2.*pB(t)-6.*cBmean.^3.*pB(t)-2.*cAmean.*cBmean.*pA(t).*pB(t)+2.*cAmean.^2.*cBmean.*pA(t).*pB(t)+2.*cAmean.*cBmean.^2.*pA(t).*pB(t)-6.*cBmean.^2.*pB(t).^2+6.*cAmean.*cBmean.^2.*pB(t).^2+6.*cBmean.^3.*pB(t).^2-6.*cAmean.*cBmean.^2.*pA(t).*pB(t).^2-cBmean.^3.*pB(t).^3+2.*cBmean.*(-1+cAmean+cBmean).*(cAmean.*(2+pA(t)-2.*pB(t))-3.*cBmean.*(-1+pB(t).^2)).*pC(t)-6.*(-1+cAmean+cBmean).^2.*(cBmean+cAmean.*pA(t)-cBmean.*pB(t)).*pC(t).^2+(-1+cAmean+cBmean).^3.*pC(t).^3)/2;
        cB(6,t)=(1+cBmean.*(-3+6.*pB(t)+cBmean.*(3+6.*(-2+pB(t)).*pB(t)+cBmean.*(-1+2.*pB(t).*(3+(-3+pB(t)).*pB(t)))))-6.*(-1+cBmean).^2.*cBmean.*pB(t).*pC(t).^2+(-1+cBmean).^3.*pC(t).^3+cAmean.*(-1+cBmean).*(-3+6.*pA(t)+cBmean.*(-1+2.*pA(t).*(-1+pB(t))+4.*pB(t))+2.*cBmean.*(2-2.*pA(t)+pB(t)).*pC(t)+6.*((-1+cBmean).*(-1+pA(t))-2.*cBmean.*pB(t)).*pC(t).^2+3.*(-1+cBmean).*pC(t).^3)-cAmean.^3.*(pA(t)-pC(t)).*(6+pA(t).^2+(-6+pC(t)).*pC(t)+pA(t).*(-6+7.*pC(t)))+cAmean.^2.*(-3.*(-1+pC(t)).*(-1-2.*pA(t).^2+(-3+pC(t)).*pC(t)+4.*pA(t).*(1+pC(t)))+cBmean.*(-6.*pA(t).^2.*(-1+pB(t)+pC(t))+(-1+pC(t)).*(1+3.*(-3+pC(t)).*pC(t)-2.*pB(t).*(2+3.*pC(t)))+2.*pA(t).*(pB(t)+2.*(-1+pC(t)).*(2+3.*pC(t))))))/2;
        cC(6,t)=1-cA(6,t)-cB(6,t);
        
        %  7 PBSL [1/0]     only gains
        cA(7,t)=(1+6.*cAmean.*pA(t)-6.*cAmean.^2.*pA(t).^2+2.*cAmean.^3.*pA(t).^3-3.*cBmean.*pB(t)-3.*cAmean.*cBmean.*pA(t).*pB(t)+6.*cAmean.^2.*cBmean.*pA(t).^2.*pB(t)+3.*cBmean.^2.*pB(t).^2-3.*cAmean.*cBmean.^2.*pA(t).*pB(t).^2-cBmean.^3.*pB(t).^3+3.*(-1+cAmean+cBmean).*(cAmean.*pA(t).*(1-2.*cAmean.*pA(t))+(-1+cBmean.*pB(t)).^2).*pC(t)-3.*(-1+cAmean+cBmean).^2.*(-1+cAmean.*pA(t)+cBmean.*pB(t)).*pC(t).^2+(-1+cAmean+cBmean).^3.*pC(t).^3)/3;
        cB(7,t)=(1-3.*cAmean.*pA(t)+3.*cAmean.^2.*pA(t).^2-cAmean.^3.*pA(t).^3+6.*cBmean.*pB(t)-3.*cAmean.*cBmean.*pA(t).*pB(t)-3.*cAmean.^2.*cBmean.*pA(t).^2.*pB(t)-6.*cBmean.^2.*pB(t).^2+6.*cAmean.*cBmean.^2.*pA(t).*pB(t).^2+2.*cBmean.^3.*pB(t).^3+3.*(-1+cAmean+cBmean).*((-1+cAmean.*pA(t)).^2+cBmean.*pB(t)-2.*cBmean.^2.*pB(t).^2).*pC(t)-3.*(-1+cAmean+cBmean).^2.*(-1+cAmean.*pA(t)+cBmean.*pB(t)).*pC(t).^2+(-1+cAmean+cBmean).^3.*pC(t).^3)/3;
        cC(7,t)=1-cA(7,t)-cB(7,t);
        
        %  8 PBSL [0/-1]    only losses
        cA(8,t)=(3-9.*cAmean+9.*cAmean.^2-3.*cAmean.^3+9.*cBmean-6.*cAmean.*cBmean-3.*cAmean.^2.*cBmean-9.*cBmean.^2-3.*cAmean.*cBmean.^2+9.*cAmean.*pA(t)-18.*cAmean.^2.*pA(t)+9.*cAmean.^3.*pA(t)+6.*cAmean.*cBmean.*pA(t)-6.*cAmean.^2.*cBmean.*pA(t)-6.*cAmean.*cBmean.^2.*pA(t)+9.*cAmean.^2.*pA(t).^2-9.*cAmean.^3.*pA(t).^2+2.*cAmean.^3.*pA(t).^3-9.*cBmean.*pB(t)+6.*cAmean.*cBmean.*pB(t)+3.*cAmean.^2.*cBmean.*pB(t)+18.*cBmean.^2.*pB(t)-6.*cAmean.*cBmean.^2.*pB(t)-9.*cBmean.^3.*pB(t)-6.*cAmean.*cBmean.*pA(t).*pB(t)+6.*cAmean.^2.*cBmean.*pA(t).*pB(t)+6.*cAmean.*cBmean.^2.*pA(t).*pB(t)-3.*cAmean.^2.*cBmean.*pA(t).^2.*pB(t)-9.*cBmean.^2.*pB(t).^2+9.*cAmean.*cBmean.^2.*pB(t).^2+9.*cBmean.^3.*pB(t).^2-3.*cAmean.*cBmean.^2.*pA(t).*pB(t).^2-cBmean.^3.*pB(t).^3+3.*(-1+cAmean+cBmean).*(cAmean.^2.*pA(t).^2+2.*cAmean.*cBmean.*(2+pA(t)-2.*pB(t))+cBmean.^2.*(3-2.*pB(t).^2)).*pC(t)-3.*(-1+cAmean+cBmean).^2.*(cAmean.*pA(t)+cBmean.*(3-2.*pB(t))).*pC(t).^2+(-1+cAmean+cBmean).^3.*pC(t).^3)/6;
        cB(8,t)=(3+9.*cAmean-9.*cAmean.^2-3.*cAmean.^2.*cBmean-9.*cAmean.*pA(t)+18.*cAmean.^2.*pA(t)-9.*cAmean.^3.*pA(t)-6.*cAmean.^2.*cBmean.*pA(t)-9.*cAmean.^2.*pA(t).^2+9.*cAmean.^3.*pA(t).^2+9.*cAmean.^2.*cBmean.*pA(t).^2-cAmean.^3.*pA(t).^3-6.*cAmean.^2.*cBmean.*pB(t)+6.*cAmean.^2.*cBmean.*pA(t).*pB(t)-3.*cAmean.^2.*cBmean.*pA(t).^2.*pB(t)+cBmean.*(9.*(-1+pB(t))+9.*cBmean.*(-1+pB(t)).^2+cBmean.^2.*(-3+(-3+pB(t)).*pB(t).*(-3+2.*pB(t))))-3.*cAmean.*cBmean.*(2.*(-1+pA(t)).*(-1+pB(t))+cBmean.*(1+2.*pB(t)+pA(t).*(-1+(-2+pB(t)).*pB(t))))+3.*(-1+cAmean+cBmean).*(cAmean.^2.*(3-2.*pA(t).^2)+cBmean.^2.*pB(t).^2+2.*cAmean.*cBmean.*(2-2.*pA(t)+pB(t))).*pC(t)-3.*(-1+cAmean+cBmean).^2.*(cAmean.*(3-2.*pA(t))+cBmean.*pB(t)).*pC(t).^2+(-1+cAmean+cBmean).^3.*pC(t).^3)/6;
        cC(8,t)=1-cA(8,t)-cB(8,t);
        
        %  9 PBSL [3/1]     pbsl-conformist
        cA(9,t)=cAmean.*(3.*cAmean-2.*cAmean.^2+2.*cBmean-2.*cAmean.*cBmean-2.*cBmean.^2+3.*pA(t)-6.*cAmean.*pA(t)+3.*cAmean.^2.*pA(t)-2.*cBmean.*pA(t)+2.*cAmean.*cBmean.*pA(t)+2.*cBmean.^2.*pA(t)-2.*cBmean.*pB(t)-cAmean.*cBmean.*pB(t)+2.*cBmean.^2.*pB(t)-cBmean.*pA(t).*pB(t)+7.*cAmean.*cBmean.*pA(t).*pB(t)-5.*cBmean.^2.*pA(t).*pB(t)-3.*cAmean.*cBmean.*pA(t).^2.*pB(t)+3.*cBmean.^2.*pA(t).*pB(t).^2+(-1+cAmean+cBmean).*(6.*pA(t)+3.*cAmean.*(1+(-4+pA(t)).*pA(t))+cBmean.*(2-5.*pA(t)-2.*pB(t))).*pC(t)+3.*(-1+cAmean+cBmean).^2.*pA(t).*pC(t).^2);
        cB(9,t)=cBmean.*(2.*cAmean-2.*cAmean.^2+3.*cBmean-2.*cAmean.*cBmean-2.*cBmean.^2-2.*cAmean.*pA(t)+2.*cAmean.^2.*pA(t)-cAmean.*cBmean.*pA(t)+3.*pB(t)-2.*cAmean.*pB(t)+2.*cAmean.^2.*pB(t)-6.*cBmean.*pB(t)+2.*cAmean.*cBmean.*pB(t)+3.*cBmean.^2.*pB(t)-cAmean.*pA(t).*pB(t)-5.*cAmean.^2.*pA(t).*pB(t)+7.*cAmean.*cBmean.*pA(t).*pB(t)+3.*cAmean.^2.*pA(t).^2.*pB(t)-3.*cAmean.*cBmean.*pA(t).*pB(t).^2+(-1+cAmean+cBmean).*(cAmean.*(2-2.*pA(t)-5.*pB(t))+3.*cBmean.*(-4+pB(t)).*pB(t)+3.*(cBmean+2.*pB(t))).*pC(t)+3.*(-1+cAmean+cBmean).^2.*pB(t).*pC(t).^2);
        cC(9,t)=1-cA(9,t)-cB(9,t);
        
        % 10 PBSL [-1/-3]   pbsl-contrarian, hates wins, hates losses even more
        cA(10,t)=(1-cAmean.^3-3.*(-1+cBmean).*cBmean+cAmean.*(-3+cBmean.*(-2-cBmean+8.*pA(t)-8.*cBmean.*pA(t)-4.*pB(t)+4.*cBmean.*pB(t)-2.*pA(t).*pB(t)+2.*cBmean.*pA(t).*pB(t)+2.*(-1+cBmean).*(2+pA(t)-2.*pB(t)).*pC(t)))+cAmean.^2.*(3+cBmean.*(-1+4.*pB(t)+4.*pC(t)-4.*pB(t).*pC(t)+2.*pA(t).*(-4+pB(t)+pC(t)))))/2;
        cB(10,t)=(-(-1+cBmean).^3+cAmean.^2.*(-3+cBmean.*(-1-8.*pB(t)+2.*pA(t).*(2+pB(t)-2.*pC(t))+2.*(2+pB(t)).*pC(t)))+cAmean.*(-1+cBmean).*(-3+cBmean.*(-1-8.*pB(t)+2.*pA(t).*(2+pB(t)-2.*pC(t))+2.*(2+pB(t)).*pC(t))))/2;
        cC(10,t)=1-cA(10,t)-cB(10,t);
        
        % 11 PBSL a la McElreath
        cA(11,t)=3.*(-1+cBmean).^2.*cBmean.*pB(t).*(-1+pC(t)).^2+cAmean.*(3.*pA(t)+cBmean.*(2.*(-1+cBmean).*(-1+pA(t))+(-1+cBmean).*(8+pA(t)).*pB(t)-3.*cBmean.*pA(t).*pB(t).^2)+(-1+cBmean).*cBmean.*(2+pA(t)-14.*pB(t)).*pC(t)-3.*(-1+cBmean).*((-1+cBmean).*pA(t)-2.*cBmean.*pB(t)).*pC(t).^2)+cAmean.^3.*(-2+3.*pC(t)-3.*pA(t).*(-1+pC(t).*(pA(t)+pC(t))))+cAmean.^2.*(3-3.*pC(t)+3.*pA(t).*(-2+pC(t).*(pA(t)+2.*pC(t)))+cBmean.*(-2+3.*pA(t).^2.*(pB(t)-pC(t))+5.*pC(t)+pA(t).*(2+pB(t)+pC(t)-6.*pC(t).^2)+pB(t).*(2+pC(t).*(-8+3.*pC(t)))));
        cB(11,t)=cBmean.*(2.*cAmean-2.*cAmean.^2+3.*cBmean-2.*cAmean.*cBmean-2.*cBmean.^2-2.*cAmean.*pA(t)+2.*cAmean.^2.*pA(t)-cAmean.*cBmean.*pA(t)+4.*cAmean.*pB(t)-cAmean.^2.*pB(t)-4.*cAmean.*cBmean.*pB(t)-cAmean.*pA(t).*pB(t)+cAmean.^2.*pA(t).*pB(t)+cAmean.*cBmean.*pA(t).*pB(t)-3.*cAmean.^2.*pA(t).^2.*pB(t)+3.*cAmean.*cBmean.*pA(t).*pB(t).^2-(-1+cAmean+cBmean).*(cAmean.*(-2+2.*pA(t)-7.*pB(t))+6.*pB(t)+3.*cBmean.*(-1+(-2+pB(t)).*pB(t))).*pC(t)-6.*(-1+cAmean+cBmean).^2.*pB(t).*pC(t).^2);
        cC(11,t)=1-cA(11,t)-cB(11,t);
        
        % somehow, matlab achieves to reach values that are infinitesimally
        % greater than 1 when certain strategies have become very rare.
        % This error screws everything, so we have to correct it
        cA(:,t)=max(min(cA(:,t),ones(nStrats,1)),zeros(nStrats,1));
        cB(:,t)=max(min(cB(:,t),ones(nStrats,1)),zeros(nStrats,1));
        cC(:,t)=max(min(cC(:,t),ones(nStrats,1)),zeros(nStrats,1));
    
    end
    
    cAlast=cA(:,tmax);
    cBlast=cB(:,tmax);
    cClast=cC(:,tmax);
    
    
    % SELECTION
    
    % fitnesses
    score=cA.*repmat(pA,nStrats,1)+cB.*repmat(pB,nStrats,1)+cC.*repmat(pC,nStrats,1);
    w(:,g)=1+selCoef*(mean(score,2)-mean(score(1,:)));
    % mean fitness
    wmean=sum(w(:,g).*x(:,g));
    % frequencies updated
    x(:,g+1)=x(:,g).*w(:,g)./wmean;
    
    % allow EXTINCTION
    % find mutants too rare
    extinctVec=x(:,g)>(1/10^9);
    % eliminate rare mutants
    x(:,g)=x(:,g).*extinctVec;
    % normalize frequency
    x(:,g)=x(:,g)./sum(x(:,g));
    
    % allow MUTATION
    if mutationRate>0
        mutationRand=rand(1,nStrats);
        % add random mutations
        x(:,g)=x(:,g)+1/(10^6)*(ones(1,nStrats).*(mutationRand<mutationRate))';
        % if a mutation occurred
        if min(mutationRand)<mutationRate
            x(:,g)=x(:,g)/sum(x(:,g));
        end
    end
    
end

maxGen=g;