function [x cA cB cC w maxGen pA pB pC]=coevo_ana_08_3(...
        xInit,...                   % initial conditions
        tmax,...                    % periods per generation
        gen,...                     % generations
        selCoef,...                 % maximum selection coefficient
        nStrats,...                 % total possible number of strategies
        vertTransm,...              % whether there is vertical transmission
        indLearn,...                % indiv. learing perf.; 0-> win-stay lose-shift
        mutationRate,...            % prob. of mutation occurring
        regime,...                  % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
        pA0,...
        pB0,...                     
        pC0,...                     % initial values of pA and pB
        incr,...                    % increment at which the environment becomes better or worse
        pincr,...                   % probability that environmental quality changes at all after each period
        dpA,...                     % shift of pA and pB
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
    pA=min(1,max(0,pA+dpA));pB=min(1,max(pB+dpB,0));pC=min(1,max(pC+dpC,0));
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
    elseif vertTransm==0.5
        adopt this later!!
        cA(:,1)=cAlast;
        cA([2 3 14 15],1)=1/3;
        cB(:,1)=cBlast;
        cB([2 3 14 15],1)=1/3;
        cC(:,1)=cClast;
        cC([2 3 14 15],1)=1/3;
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
        
        % conformists tn=3
        cConf3A=cAmean*(-2*cAmean^2 + cAmean*(3 - 2*cBmean) - 2*(-1 + cBmean)*cBmean);
        cConf3B=cBmean*(3*cBmean - 2*(cAmean^2 + cAmean*(-1 + cBmean) + cBmean^2));
        cConf3C=1-cConf3A-cConf3B;
        cA(2,t)=cConf3A;
        cB(2,t)=cConf3B;
        cC(2,t)=cConf3C;
        
        % conformists tn=7
        cConf7A=cAmean^3*(cAmean*(35 - 2*cAmean*(42 + 5*cAmean*(-7 + 2*cAmean))) - 70*(-1 + cAmean)^3*cBmean + 140*(-1 + cAmean)*cBmean^3 + 70*cBmean^4);
        cConf7B=cBmean^3*(70*cAmean^4 + 140*cAmean^3*(-1 + cBmean) - 70*cAmean*(-1 + cBmean)^3 + cBmean*(35 - 2*cBmean*(42 + 5*cBmean*(-7 + 2*cBmean))));
        cConf7C=1-cConf7A-cConf7B;
        cA(3,t)=cConf7A;
        cB(3,t)=cConf7B;
        cC(3,t)=cConf7C;
        
        % mixed 2/3-1/3 tn=3
        if x(4,g)>0;
            cA(4,t)=2/3*cAind+1/3*cConf3A;
            cB(4,t)=2/3*cBind+1/3*cConf3B;
            cC(4,t)=1-cA(4,t)-cB(4,t);
        end
        
        % mixed 2/3-1/3 tn=7
        if x(5,g)>0;
            cA(5,t)=2/3*cAind+1/3*cConf7A;
            cB(5,t)=2/3*cBind+1/3*cConf7B;
            cC(5,t)=1-cA(5,t)-cB(5,t);
        end
        
        % mixed 1/3-2/3 tn=3
        if x(6,g)>0;
            cA(6,t)=1/3*cAind+2/3*cConf3A;
            cB(6,t)=1/3*cBind+2/3*cConf3B;
            cC(6,t)=1-cA(6,t)-cB(6,t);
        end        
          
        % mixed 1/3-2/3 tn=7
        if x(7,g)>0;
            cA(7,t)=1/3*cAind+2/3*cConf7A;
            cB(7,t)=1/3*cBind+2/3*cConf7B;
            cC(7,t)=1-cA(7,t)-cB(7,t);
        end
                
        % PBSL type 1, equal weight tn=3
        if x(8,g)>0;
        cA(8,t)=(1 - 3*cAmean + 3*cAmean^2 - cAmean^3 + 3*cBmean - 2*cAmean*cBmean - cAmean^2*cBmean - 3*cBmean^2 - cAmean*cBmean^2 + 6*cAmean*pA(t) - 6*cAmean^2*pA(t) - 4*cAmean*cBmean*pA(t) + 4*cAmean^2*cBmean*pA(t) + 4*cAmean*cBmean^2*pA(t) + 6*cAmean^3*pA(t)^2 - 4*cAmean^3*pA(t)^3 - 6*cBmean*pB(t) + 8*cAmean*cBmean*pB(t) - 2*cAmean^2*cBmean*pB(t) + 6*cBmean^2*pB(t) - 2*cAmean*cBmean^2*pB(t) - 2*cAmean*cBmean*pA(t)*pB(t) - 4*cAmean^2*cBmean*pA(t)*pB(t) + 2*cAmean*cBmean^2*pA(t)*pB(t) + 6*cAmean^2*cBmean*pA(t)^2*pB(t) - 3*cBmean^3*pB(t)^2 - 6*cAmean*cBmean^2*pA(t)*pB(t)^2 + 2*cBmean^3*pB(t)^3 + 2*(-1 + cAmean + cBmean)*(-3*cAmean^2*(-1 + pA(t))*pA(t) + cBmean*(3 - 3*pB(t) + cAmean*(-1 + pA(t) + pB(t))))*pC(t) - 3*(-1 + cAmean + cBmean)^2*(1 - cAmean - cBmean + 2*cAmean*pA(t))*pC(t)^2 - 2*(-1 + cAmean + cBmean)^3*pC(t)^3)/2;
        cB(8,t)=(1 + 3*cAmean - 3*cAmean^2 - 3*cBmean - 2*cAmean*cBmean - cAmean^2*cBmean + 3*cBmean^2 - cAmean*cBmean^2 - cBmean^3 - 6*cAmean*pA(t) + 6*cAmean^2*pA(t) + 8*cAmean*cBmean*pA(t) - 2*cAmean^2*cBmean*pA(t) - 2*cAmean*cBmean^2*pA(t) - 3*cAmean^3*pA(t)^2 + 2*cAmean^3*pA(t)^3 + 6*cBmean*pB(t) - 4*cAmean*cBmean*pB(t) + 4*cAmean^2*cBmean*pB(t) - 6*cBmean^2*pB(t) + 4*cAmean*cBmean^2*pB(t) - 2*cAmean*cBmean*pA(t)*pB(t) + 2*cAmean^2*cBmean*pA(t)*pB(t) - 4*cAmean*cBmean^2*pA(t)*pB(t) - 6*cAmean^2*cBmean*pA(t)^2*pB(t) + 6*cBmean^3*pB(t)^2 + 6*cAmean*cBmean^2*pA(t)*pB(t)^2 - 4*cBmean^3*pB(t)^3 + 2*(-1 + cAmean + cBmean)*(cAmean*(-3 + cBmean)*(-1 + pA(t)) + cAmean*cBmean*pB(t) - 3*cBmean^2*(-1 + pB(t))*pB(t))*pC(t) + 3*(-1 + cAmean + cBmean)^2*(-1 + cAmean + cBmean - 2*cBmean*pB(t))*pC(t)^2 - 2*(-1 + cAmean + cBmean)^3*pC(t)^3)/2;
        cC(8,t)=1-cA(8,t)-cB(8,t);
        end
        
        % pbsl equal weights tn=7
        if x(9,g)>0;
        cA(9,t)=0;
        cB(9,t)=0;
        cC(9,t)=0;
        end
        
        % pbsl gains/losses 3/1 tn=3
        if x(10,g)>0;
        cA(10,t)=(cAmean^3*(2*pA(t)^3 + (-1 + pC(t))^3 - 12*pA(t)^2*pC(t) - 6*pA(t)*(-2 + pC(t))*pC(t)) - (-1 + cBmean*(pB(t) - pC(t)) + pC(t))*(cBmean*(-3 + 5*pB(t) - 2*pC(t))*(-1 + pC(t)) + (-1 + pC(t))^2 + cBmean^2*(-3 + pB(t)^2 + pB(t)*(3 - 5*pC(t)) + pC(t)*(3 + pC(t)))) + cAmean^2*(3*(-(-1 + pC(t))^3 + 4*pA(t)*(-1 + pC(t))*pC(t) + pA(t)^2*(-2 + 4*pC(t))) + cBmean*(12*pA(t)^2*(pB(t) - pC(t)) + 2*pA(t)*(2 - 5*pB(t) + (7 - 6*pC(t))*pC(t)) + (-1 + pC(t))*(1 + 2*pB(t) - 6*pB(t)*pC(t) + 3*pC(t)^2))) + cAmean*(3*(-1 + pC(t))*((-1 + pC(t))^2 - 2*pA(t)*(1 + pC(t))) + 2*cBmean*(-(pA(t)*(2 + pB(t) + pC(t) - 6*pC(t)^2)) + (-1 + pC(t))*(1 - 4*pB(t) + 6*pB(t)*pC(t) - 3*pC(t)^2)) + cBmean^2*(2*pA(t)*(2 + pB(t) - 3*pB(t)^2 + pC(t) - 3*pC(t)^2) + (-1 + pC(t))*(1 + 6*pB(t)^2 + 3*pC(t)*(2 + pC(t)) - 4*pB(t)*(1 + 3*pC(t))))))/2;
        cB(10,t)=-(cAmean^3*(-1 + pA(t))^3)/2 + 3*cAmean^2*cBmean*(-1 + pA(t))^2*pB(t) - 6*cAmean^2*cBmean*(-1 + pA(t))*pA(t)*pB(t) + 6*cAmean*cBmean^2*(-1 + pA(t))*(-1 + pB(t))*pB(t) + 3*cBmean^3*(-1 + pB(t))^2*pB(t) - 3*cAmean*cBmean^2*(-1 + pA(t))*pB(t)^2 + 3*cAmean*cBmean^2*pA(t)*pB(t)^2 - 3*cBmean^3*(-1 + pB(t))*pB(t)^2 + cBmean^3*pB(t)^3 + 3*cAmean^2*(-1 + cAmean + cBmean)*(-1 + pA(t))^2*(-1 + pC(t)) + 2*cAmean*cBmean*(-1 + cAmean + cBmean)*(-1 + pA(t))*(-1 + pB(t))*(-1 + pC(t)) - 6*cAmean*cBmean*(-1 + cAmean + cBmean)*(-1 + pA(t))*pB(t)*(-1 + pC(t)) + 3*cAmean*cBmean*(-1 + cAmean + cBmean)*pA(t)*pB(t)*(-1 + pC(t)) - 6*cBmean^2*(-1 + cAmean + cBmean)*(-1 + pB(t))*pB(t)*(-1 + pC(t)) + 3*cBmean^2*(-1 + cAmean + cBmean)*pB(t)^2*(-1 + pC(t)) - 3*cAmean*(-1 + cAmean + cBmean)^2*(-1 + pA(t))*(-1 + pC(t))^2 + 3*cBmean*(-1 + cAmean + cBmean)^2*pB(t)*(-1 + pC(t))^2 + ((-1 + cAmean + cBmean)^3*(-1 + pC(t))^3)/2 + 3*cAmean*cBmean*(-1 + cAmean + cBmean)*(-1 + pA(t))*pB(t)*pC(t) - 2*cAmean*cBmean*(-1 + cAmean + cBmean)*pA(t)*pB(t)*pC(t) - 3*cBmean^2*(-1 + cAmean + cBmean)*pB(t)^2*pC(t) - 6*cBmean*(-1 + cAmean + cBmean)^2*pB(t)*(-1 + pC(t))*pC(t);
        cC(10,t)=1-cA(10,t)-cB(10,t);
        end
        
        % pbsl gains/losses 3/1 tn 7
        if x(11,g)>0
        cA(11,t)=0;
        cB(11,t)=0;
        cC(11,t)=0;
        end
        
        % pbsl gains/losses 1/3 tn=3
        if x(12,g)>0;
        cA(12,t)=(1 - 3*cAmean + 3*cAmean^2 - cAmean^3 + 3*cBmean - 2*cAmean*cBmean - cAmean^2*cBmean - 3*cBmean^2 - cAmean*cBmean^2 + 6*cAmean*pA(t) - 12*cAmean^2*pA(t) + 6*cAmean^3*pA(t) - 4*cAmean*cBmean*pA(t) + 4*cAmean^2*cBmean*pA(t) + 4*cAmean*cBmean^2*pA(t) + 6*cAmean^2*pA(t)^2 - 6*cAmean^3*pA(t)^2 + 2*cAmean^3*pA(t)^3 - 6*cBmean*pB(t) + 8*cAmean*cBmean*pB(t) - 2*cAmean^2*cBmean*pB(t) + 12*cBmean^2*pB(t) - 8*cAmean*cBmean^2*pB(t) - 6*cBmean^3*pB(t) - 2*cAmean*cBmean*pA(t)*pB(t) + 2*cAmean^2*cBmean*pA(t)*pB(t) + 2*cAmean*cBmean^2*pA(t)*pB(t) - 6*cBmean^2*pB(t)^2 + 6*cAmean*cBmean^2*pB(t)^2 + 6*cBmean^3*pB(t)^2 - 6*cAmean*cBmean^2*pA(t)*pB(t)^2 - cBmean^3*pB(t)^3 + 2*cBmean*(-1 + cAmean + cBmean)*(cAmean*(2 + pA(t) - 2*pB(t)) - 3*cBmean*(-1 + pB(t)^2))*pC(t) - 6*(-1 + cAmean + cBmean)^2*(cBmean + cAmean*pA(t) - cBmean*pB(t))*pC(t)^2 + (-1 + cAmean + cBmean)^3*pC(t)^3)/2;
        cB(12,t)=(1 + cBmean*(-3 + 6*pB(t) + cBmean*(3 + 6*(-2 + pB(t))*pB(t) + cBmean*(-1 + 2*pB(t)*(3 + (-3 + pB(t))*pB(t))))) - 6*(-1 + cBmean)^2*cBmean*pB(t)*pC(t)^2 + (-1 + cBmean)^3*pC(t)^3 + cAmean*(-1 + cBmean)*(-3 + 6*pA(t) + cBmean*(-1 + 2*pA(t)*(-1 + pB(t)) + 4*pB(t)) + 2*cBmean*(2 - 2*pA(t) + pB(t))*pC(t) + 6*((-1 + cBmean)*(-1 + pA(t)) - 2*cBmean*pB(t))*pC(t)^2 + 3*(-1 + cBmean)*pC(t)^3) - cAmean^3*(pA(t) - pC(t))*(6 + pA(t)^2 + (-6 + pC(t))*pC(t) + pA(t)*(-6 + 7*pC(t))) + cAmean^2*(-3*(-1 + pC(t))*(-1 - 2*pA(t)^2 + (-3 + pC(t))*pC(t) + 4*pA(t)*(1 + pC(t))) + cBmean*(-6*pA(t)^2*(-1 + pB(t) + pC(t)) + (-1 + pC(t))*(1 + 3*(-3 + pC(t))*pC(t) - 2*pB(t)*(2 + 3*pC(t))) + 2*pA(t)*(pB(t) + 2*(-1 + pC(t))*(2 + 3*pC(t))))))/2;
        cC(12,t)=1-cA(12,t)-cB(12,t);
        end
        
        % pbsl gains/losses 1/3 tn 7
        if x(13,g)>0
        cA(13,t)=0;
        cB(13,t)=0;
        cC(13,t)=0;
        end
        
        % pbsl conf gains/losses 3/1 tn 3
        if x(14,g)>0
        cA(14,t)=cAmean*(3*cAmean - 2*cAmean^2 + 2*cBmean - 2*cAmean*cBmean - 2*cBmean^2 + 3*pA(t) - 6*cAmean*pA(t) + 3*cAmean^2*pA(t) - 2*cBmean*pA(t) + 2*cAmean*cBmean*pA(t) + 2*cBmean^2*pA(t) - 2*cBmean*pB(t) - cAmean*cBmean*pB(t) + 2*cBmean^2*pB(t) - cBmean*pA(t)*pB(t) + 7*cAmean*cBmean*pA(t)*pB(t) - 5*cBmean^2*pA(t)*pB(t) - 3*cAmean*cBmean*pA(t)^2*pB(t) + 3*cBmean^2*pA(t)*pB(t)^2 + (-1 + cAmean + cBmean)*(6*pA(t) + 3*cAmean*(1 + (-4 + pA(t))*pA(t)) + cBmean*(2 - 5*pA(t) - 2*pB(t)))*pC(t) + 3*(-1 + cAmean + cBmean)^2*pA(t)*pC(t)^2);
        cB(14,t)=cBmean*(2*cAmean - 2*cAmean^2 + 3*cBmean - 2*cAmean*cBmean - 2*cBmean^2 - 2*cAmean*pA(t) + 2*cAmean^2*pA(t) - cAmean*cBmean*pA(t) + 3*pB(t) - 2*cAmean*pB(t) + 2*cAmean^2*pB(t) - 6*cBmean*pB(t) + 2*cAmean*cBmean*pB(t) + 3*cBmean^2*pB(t) - cAmean*pA(t)*pB(t) - 5*cAmean^2*pA(t)*pB(t) + 7*cAmean*cBmean*pA(t)*pB(t) + 3*cAmean^2*pA(t)^2*pB(t) - 3*cAmean*cBmean*pA(t)*pB(t)^2 + (-1 + cAmean + cBmean)*(cAmean*(2 - 2*pA(t) - 5*pB(t)) + 3*cBmean*(-4 + pB(t))*pB(t) + 3*(cBmean + 2*pB(t)))*pC(t) + 3*(-1 + cAmean + cBmean)^2*pB(t)*pC(t)^2);
        cC(14,t)=1-cA(14,t)-cB(14,t);
        end
        
        % pbsl conf gains/losses 3/1 tn 7
        if x(15,g)>0
        cA(15,t)=0;
        cB(15,t)=0;
        cC(15,t)=0;
        end
        
        % pbsl only gains tn 3
        if x(16,g)>0
        cA(16,t)=(1 + 6*cAmean*pA(t) - 6*cAmean^2*pA(t)^2 + 2*cAmean^3*pA(t)^3 - 3*cBmean*pB(t) - 3*cAmean*cBmean*pA(t)*pB(t) + 6*cAmean^2*cBmean*pA(t)^2*pB(t) + 3*cBmean^2*pB(t)^2 - 3*cAmean*cBmean^2*pA(t)*pB(t)^2 - cBmean^3*pB(t)^3 + 3*(-1 + cAmean + cBmean)*(cAmean*pA(t)*(1 - 2*cAmean*pA(t)) + (-1 + cBmean*pB(t))^2)*pC(t) - 3*(-1 + cAmean + cBmean)^2*(-1 + cAmean*pA(t) + cBmean*pB(t))*pC(t)^2 + (-1 + cAmean + cBmean)^3*pC(t)^3)/3;
        cB(16,t)=(1 - 3*cAmean*pA(t) + 3*cAmean^2*pA(t)^2 - cAmean^3*pA(t)^3 + 6*cBmean*pB(t) - 3*cAmean*cBmean*pA(t)*pB(t) - 3*cAmean^2*cBmean*pA(t)^2*pB(t) - 6*cBmean^2*pB(t)^2 + 6*cAmean*cBmean^2*pA(t)*pB(t)^2 + 2*cBmean^3*pB(t)^3 + 3*(-1 + cAmean + cBmean)*((-1 + cAmean*pA(t))^2 + cBmean*pB(t) - 2*cBmean^2*pB(t)^2)*pC(t) - 3*(-1 + cAmean + cBmean)^2*(-1 + cAmean*pA(t) + cBmean*pB(t))*pC(t)^2 + (-1 + cAmean + cBmean)^3*pC(t)^3)/3;
        cC(16,t)=1-cA(16,t)-cB(16,t);
        end
        
        % pbsl only gains tn 7
        if x(17,g)>0
        cA(17,t)=0;
        cB(17,t)=0;
        cC(17,t)=0;
        end
        
        % pbsl McElreath
        if x(18,g)>0
        cA(18,t)=3*(-1 + cBmean)^2*cBmean*pB(t)*(-1 + pC(t))^2 + cAmean*(3*pA(t) + cBmean*(2*(-1 + cBmean)*(-1 + pA(t)) + (-1 + cBmean)*(8 + pA(t))*pB(t) - 3*cBmean*pA(t)*pB(t)^2) + (-1 + cBmean)*cBmean*(2 + pA(t) - 14*pB(t))*pC(t) - 3*(-1 + cBmean)*((-1 + cBmean)*pA(t) - 2*cBmean*pB(t))*pC(t)^2) + cAmean^3*(-2 + 3*pC(t) - 3*pA(t)*(-1 + pC(t)*(pA(t) + pC(t)))) + cAmean^2*(3 - 3*pC(t) + 3*pA(t)*(-2 + pC(t)*(pA(t) + 2*pC(t))) + cBmean*(-2 + 3*pA(t)^2*(pB(t) - pC(t)) + 5*pC(t) + pA(t)*(2 + pB(t) + pC(t) - 6*pC(t)^2) + pB(t)*(2 + pC(t)*(-8 + 3*pC(t)))));
        cB(18,t)=0;cBmean*(2*cAmean - 2*cAmean^2 + 3*cBmean - 2*cAmean*cBmean - 2*cBmean^2 - 2*cAmean*pB(t) + 2*cAmean^2*pB(t) - cAmean*cBmean*pB(t) + 4*cAmean*pB(t) - cAmean^2*pB(t) - 4*cAmean*cBmean*pB(t) - cAmean*pB(t)*pB(t) + cAmean^2*pB(t)*pB(t) + cAmean*cBmean*pB(t)*pB(t) - 3*cAmean^2*pB(t)^2*pB(t) + 3*cAmean*cBmean*pB(t)*pB(t)^2 - (-1 + cAmean + cBmean)*(cAmean*(-2 + 2*pB(t) - 7*pB(t)) + 6*pB(t) + 3*cBmean*(-1 + (-2 + pB(t))*pB(t)))*pC(t) - 6*(-1 + cAmean + cBmean)^2*pB(t)*pC(t)^2);
        cC(18,t)=1-cA(18,t)-cB(18,t);
        end
        
        % pbsl hybrid
        if x(19,g)>0
        cA(19,t)=0;
        cB(19,t)=0;
        cC(19,t)=0;
        end
        
        % somehow, matlab achieves to reach values that are infinitesimally
        % greater than 1 when certain strategies have become very rare.
        % This error screws everything, so we have to correct it
        cA(:,t)=min(cA(:,t),ones(nStrats,1));
        cB(:,t)=min(cB(:,t),ones(nStrats,1));
        cC(:,t)=min(cC(:,t),ones(nStrats,1));
        cA(:,t)=max(cA(:,t),zeros(nStrats,1));
        cB(:,t)=max(cB(:,t),zeros(nStrats,1));
        cC(:,t)=max(cC(:,t),zeros(nStrats,1));
        
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