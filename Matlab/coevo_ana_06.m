function [x c w maxGen pA pB]=coevo_ana_06(...
        xInit,...                   % initial conditions
        tmax,...                    % periods per generation
        gen,...                     % generations
        selCoef,...                 % maximum selection coefficient
        nStrats,...                 % total possible number of strategies
        vertTransm,...              % whether there is vertical transmission
        regime,...                  % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
        pA0,...
        pB0,...                     % initial values of pA and pB
        incr,...                    % increment at which the environment becomes better or worse
        pincr);                     % probability that environmental quality changes at all after each period
        
            
% INITIALIZATION
x=zeros(nStrats,gen);       % frequency matrix
c=zeros(nStrats,tmax);      % percent correct choices matrix
w=zeros(nStrats,gen);       % fitness matrix

x(:,1)=xInit;

extinctVec=zeros(nStrats,1);

for g=1:gen
    
    % ENVIRONMENT
    % routine to determine pA and pB
    [pA,pB] = randomenvironment2(tmax,regime,incr,pincr,pA0,pB0);
    % we define C is the better and F as the worse choice
    pC=pA;
    pF=pB;
    % remember last pA and pB
    pA0=pA(tmax);pB0=pB(tmax);
%     % when do switches occur? if switch -> tSwitch=-1
%     pDif=sign(pA-pB);
%     tSwitch=[1 pDif(2:tmax).*pDif(1:tmax-1)];
    
    % FIRST PERIOD'S CHOICE OF THE STRATEGIES
    if g==1||vertTransm==0     % in very first generation
        c(:,1)=1/2;
    else
        c(:,1)=clast;
    end
    
    % CHOICES OF THE STRATEGIES
    for t=2:tmax
                
%         % in case of environmental change
%         if tSwitch(t-1)<0
%             c(:,t-1)=1-c(:,t-1);
%         end
        
        % mean proportion of correct choices
%         cmean=x(1,g)*c(1,t-1)+x(2,g)*c(2,t-1)+x(3,g)*c(3,t-1)+...
%             x(4,g)*c(4,t-1)+x(5,g)*c(5,t-1)+x(6,g)*c(6,t-1)+...
%             x(7,g)*c(7,t-1)+x(8,g)*c(8,t-1)+x(9,g)*c(9,t-1)+...
%             x(10,g)*c(10,t-1)+x(11,g)*c(11,t-1)+x(12,g)*c(12,t-1);
        cmean=sum(x(:,g).*c(:,t-1));
        
        % individual learners
        cind=(1-pF(t))/(2-pC(t)-pF(t));
        c(1,t)=cind;
        
        % conformists
        c(2,t)=(3-2*cmean)*cmean^2;

        % PBSL type 1, equal weight tn=3
        if x(3,t)>0;
        c(3,t)=-((1 - pF(t) + cmean*(-1 + pC(t) + pF(t))).^2.*(-1 - ...
            2*pF(t) + 2*cmean*(-1 + pC(t) + pF(t))));
        end

        
        % PBSL type 2, more weight gains, tn=3
        if x(4,t)>0;
        c(4,t)=1 + cmean^2*(-1 + pC(t))^2*(-3 + cmean*(2 + ...
            pC(t))) - 3*(-1 + cmean)*(-1 + cmean*(cmean + ...
            2*(1 + cmean*(-2 + pC(t)))*pC(t)))*pF(t) - 3*(-1 + ...
            cmean)^2*(-1 + 2*cmean*pC(t))*pF(t)^2 + (-1 + ...
            cmean)^3*pF(t)^3;
        end

        
        % PBSL type 3 a la McElreath et al 2008
        if x(5,t)>0;
        c(5,t)=cmean*(-3*(-1 + cmean)*cmean*pC(t)^2*pF(t) + ...
            cmean*(3 - 2*cmean + 3*(-1 + cmean)*pF(t)) - ...
            3*(-1 + cmean)^2*pC(t)*(-1 + pF(t)^2));
        end

        
        % opportunistic individual learner
        if x(6,t)>0;
        c(6,t)=(cind-cind*pF(t)+(3-2*cmean)*cmean^2*pF(t))/...
            (1+(cind+cmean^2*(-3+2*cmean))*(pC(t)-pF(t)));
        end
        
          
        % opportunistic conformist
        if x(7,t)>0;
            if (pC(t)==1&pF(t)==0)
                % this check must be included to avert dividing by 0
                c(7,t)=1;
            else
                c(7,t)=(-cmean^2*(-3+2*cmean)*(pF(t)-1)-cind*pF(t))/...
                    (-1+(cind+cmean^2*(-3+2*cmean))*(pC(t)-pF(t)));            
            end
        end
        
        % pbsl tn=7, loss/gains = 3/1
        if x(8,t)>0;
        c(8,t)=1 - cmean.^4.*(-1 + pC(t)).^2.*(35.*(-1 + pC(t)).^2 -...
            70*cmean.^2.*(-1 + 4*pC(t) - 4*pC(t).^3 + pC(t).^4) - ...
            42*cmean.*(2 + pC(t).*(-6 + pC(t) + 3*pC(t).^2)) + ...
            cmean.^3.*(-20 + pC(t).*(100 + pC(t).*(10 + 3*pC(t)*...
            (-50 + pC(t).*(25 + 2*pC(t))))))) + 35.*(-1 + cmean)...
            *cmean.^3.*(-1 + pC(t)).^2.*(4 - 12*cmean.*(-1 + pC(t)).^2 -...
            4*pC(t) + 12*cmean.^2.*(1 + (-1 + pC(t))*pC(t).*(3 + pC(t))) +...
            cmean.^3.*(-4 + pC(t).*(16 + 3*pC(t).*(-4 + (-4 + pC(t))*pC(t)))))*pF(t) +...
            105.*(-1 + cmean).^2*cmean.^3.*(-1 + pC(t)).*(-2 + 10*pC(t) -...
            8*pC(t).^2 - 2*cmean.*(-1 + pC(t)).*(2 + 5.*(-2 + pC(t))*pC(t)) +...
            cmean.^2.*(-2 + pC(t).*(14 + pC(t).*(-22 + pC(t).*(8 + 3*pC(t))))))*pF(t).^2 -...
            70.*(-1 + cmean).^3*cmean.^2.*(-1 + pC(t)).*(3 - 3*pC(t) + ...
            cmean.*(-4 + 10*pC(t).*(-1 + 2*pC(t)) + cmean.*(1 + ...
            5*pC(t).*(3 + 2.*(-3 + pC(t))*pC(t)))))*pF(t).^3 - 35.*(-1 + ...
            cmean).^4*cmean.*(3 - 15*cmean.*(-1 + pC(t)).^2 - 3*pC(t) + ...
            cmean.^2.*(11 + 5*pC(t).*(-3 + pC(t).*(-3 + 4*pC(t)))))*pF(t).^4 + ...
            21.*(-1 + cmean).^5*cmean.*(-8 + 10*pC(t) + cmean.*(14 + ...
            15.*(-2 + pC(t))*pC(t)))*pF(t).^5 + 7.*(-1 + cmean).^6.*(-1 + ...
            3*cmean.*(-3 + 5*pC(t)))*pF(t).^6 - 6.*(-1 + cmean).^7*pF(t).^7;
        end  
        
        % pbsl tn=7, equal weights
        if x(9,t)>0;
        c(9,t)=-((1 - pF(t) + cmean.*(-1 + pC(t) + pF(t))).^4.*(-1 + ...
            20*cmean.^3.*(-1 + pC(t) + pF(t)).^3 - 10*cmean.^2.*(-1 + ...
            pC(t) + pF(t)).^2.*(1 + 6*pF(t)) - 2*pF(t).*(2 + 5*pF(t).*(1 + 2*pF(t))) + ...
            4*cmean.*(-1 + pC(t) + pF(t)).*(1 + 5*pF(t).*(1 + 3*pF(t)))));
        end
        
        % pbsl tn=7, loss/gains = 1/3
        if x(10,t)>0;
        c(10,t)=1 - cmean.^4.*(-1 + pC(t)).^4.*(35 - 21*cmean.*(4 + ...
            pC(t)) - 35*cmean.^2.*(-2 + (-2 + pC(t))*pC(t)) + cmean.^3.*(-20 + ...
            3*pC(t).*(-15 + 2*pC(t).*(4 + pC(t))))) + 35.*(-1 + cmean)*cmean.^2.*(-1 + ...
            pC(t)).^2.*(3 + cmean.*(-4.*(2 + pC(t)) + cmean.*(6 - 12.*(-2 + ...
            pC(t))*pC(t) + 6*cmean*pC(t).*(-5 + pC(t).*(2 + pC(t))) + cmean.^2.*(-1 + ...
            pC(t).*(10 + 3*pC(t).*(1 + (-4 + pC(t))*pC(t)))))))*pF(t) + 21.*(-1 + ...
            cmean).^2.*(-1 + 5*cmean*pC(t) + 10*cmean.^2.*(2 + 3.*(-2 + ...
            pC(t))*pC(t)) - 10*cmean.^3.*(4 - 9*pC(t) + 4*pC(t).^3) + 5*cmean.^4.*(5 - ...
            2*pC(t).^2.*(18 + 5.*(-4 + pC(t))*pC(t))) + cmean.^5.*(-4 + ...
            5*pC(t).*(-6 + pC(t).*(24 - 20*pC(t) + 3*pC(t).^3))))*pF(t).^2 - 70.*(-1 + ...
            cmean).^3.*(1 - 6*cmean*pC(t) - 9*cmean.^2.*(1 + 2.*(-2 + pC(t))*pC(t)) + ...
            2*cmean.^4.*(-1 + pC(t)).*(2 + 5.*(-3 + pC(t))*pC(t).^2) + 4*cmean.^3.*(3 - ...
            9*pC(t) + 5*pC(t).^3))*pF(t).^3 - 35.*(-1 + cmean).^4.*(3 + ...
            2*cmean.*(-9*pC(t) + cmean.*(-6 - 15.*(-2 + pC(t))*pC(t) + cmean.*(4 - ...
            15*pC(t) + 10*pC(t).^3))))*pF(t).^4 + 21.*(-1 + cmean).^5.*(-4 + ...
            5*cmean.*(4*pC(t) + cmean.*(1 + 3.*(-2 + pC(t))*pC(t))))*pF(t).^5 + ...
            35.*(-1 + cmean).^6.*(-1 + 3*cmean*pC(t))*pF(t).^6 - 6.*(-1 + ...
            cmean).^7*pF(t).^7;
        end
        
        % payoff conformist, weight A+2,a+1,B-2,b-1
        if x(11,t)>0;
        c(11,t)=(cmean.^3.*(105.*(-1 + cmean).^3*cmean*pC(t).^4.*(-1 + ...
            pF(t))*pF(t).*(-1 + 3*pF(t)) + cmean.*(70 - 4*cmean.*(42 + ...
            5*cmean.*(-7 + 2*cmean)) + 105.*(-1 + cmean).^3*pF(t) - ...
            35.*(-1 + cmean).^3*pF(t).^3) + 105.*(-1 + cmean).^3*pC(t).*(-1 + ...
            pF(t)).*(1 - cmean + pF(t).*(-3 + 7*cmean + (-1 + cmean).*(-3 + ...
            pF(t))*pF(t))) - 210.*(-1 + cmean).^3*pC(t).^2*pF(t).*(-2.*(-1 + pF(t)).^3 + ...
            cmean.*(-5 + 2*pF(t).*(6 + (-4 + pF(t))*pF(t)))) + 35.*(-1 + ...
            cmean).^3*pC(t).^3.*(1 - cmean + pF(t).*(-12*cmean + pF(t).*(-12 + ...
            48*cmean + pF(t).*(20 - 40*cmean + 9.*(-1 + cmean)*pF(t)))))))/2;
        end
        
        % only count wins: weight A+1,a+0,B-1,b-0
        if x(12,t)>0;
        c(12,t)=(1 + cmean*pC(t).*(7 + cmean*pC(t).*(-21 + ...
            cmean*pC(t).*(35 + cmean*pC(t).*(-35 + cmean*pC(t).*(21 + ...
            cmean*pC(t).*(-7 + cmean*pC(t))))))) - 7.*(-1 + cmean).*(-1 + ...
            cmean*pC(t)).^5.*(1 + 5*cmean*pC(t))*pF(t) + 21.*(-1 + cmean).^2.*(-1 + ...
            cmean*pC(t)).^3.*(-1 + cmean*pC(t).*(2 + 9*cmean*pC(t)))*pF(t).^2 - ...
            35.*(-1 + cmean).^3.*(-1 + cmean*pC(t).*(8 - 12*cmean*pC(t) + ...
            5*cmean.^3*pC(t).^3))*pF(t).^3 - 35.*(-1 + cmean).^4.*(-1 + ...
            cmean*pC(t).*(9 + 5*cmean*pC(t).*(-3 + cmean*pC(t))))*pF(t).^4 + ...
            21.*(-1 + cmean).^5.*(1 + cmean*pC(t).*(-8 + 9*cmean*pC(t)))*pF(t).^5 - ...
            7.*(-1 + cmean).^6.*(-1 + 5*cmean*pC(t))*pF(t).^6 + (-1 + ...
            cmean).^7*pF(t).^7)/2;
        end
        
        % somehow, matlab achieves to reach values that are infinitesimally
        % greater than 1 when certain strategies have become very rare.
        % This error screws everything, so we have to correct it
        c(:,t)=min(c(:,t),ones(nStrats,1));
    end
    
    clast=c(:,tmax);
    
    % SELECTION
    
    % fitnesses
    score=c.*repmat(pC,nStrats,1)+(1-c).*repmat(pF,nStrats,1);
    w(:,g)=1+selCoef*(mean(score,2)-mean(score(1,:)));
    % mean fitness
    wmean=sum(w(:,g).*x(:,g));
    % frequencies updated
    x(:,g+1)=x(:,g).*w(:,g)./wmean;
    
    % allow extinction
    % find mutants too rare
    extinctVec=x(:,g)>(1/10^9);
    % eliminate rare mutants
    x(:,g)=x(:,g).*extinctVec;
    % normalize frequency
    x(:,g)=x(:,g)./sum(x(:,g));
    
end

maxGen=g;