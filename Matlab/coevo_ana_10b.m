function [x c w maxGen pA pB]=coevo_ana_08(...
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
        pB0,...                     % initial values of pA and pB
        incr,...                    % increment at which the environment becomes better or worse
        pincr,...                   % probability that environmental quality changes at all after each period
        dpA,...                     % shift of pA and pB
        dpB);
        
            
% INITIALIZATION
x=zeros(nStrats,gen);       % frequency matrix
c=zeros(nStrats,tmax);      % percent correct choices matrix
w=zeros(nStrats,gen);       % fitness matrix

x(:,1)=xInit;

extinctVec=zeros(nStrats,1);

for g=1:gen
    
    % ENVIRONMENT
    % routine to determine pA and pB
    [pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);
    % possible shift
    pA=min(1-incr,max(incr,pA+dpA));pB=min(1-incr,max(pB+dpB,incr));
    % we define C is the better and F as the worse choice
    pC=pA;
    pF=pB;
    % remember last pA and pB
    pA0=pA(tmax);pB0=pB(tmax);
    
    % FIRST PERIOD'S CHOICE OF THE STRATEGIES
    if g==1||vertTransm==0     % in very first generation
        c(:,1)=1/2;
    elseif vertTransm==1
        c(:,1)=clast;
    elseif vertTransm==0.5
        c(:,1)=clast;
        c([2 3 14 15],1)=1/2;
    end
    
    % CHOICES OF THE STRATEGIES
    for t=2:tmax
        cmean=sum(x(:,g).*c(:,t-1));

        % individual learners
        if indLearn==0
            cind=(1-pF(t))/(2-pC(t)-pF(t));
        else
            cind=((1-2*(1-indLearn))*sign(pC(t)-pF(t))+1)/2;
        end
        c(1,t)=cind;

        % 02 conformists
        cconf7=cmean.^4.*(35-2.*cmean.*(42+5.*cmean.*(-7+2.*cmean)));
        c(2,t)=cconf7;

        % 03 contrarian
        c(3,t)=(-1+cmean).^4.*(1+2.*cmean.*(2+5.*cmean.*(1+2.*cmean)));

        %  4 PBSL [1/-1]    equal weights
        c(4,t)=-((1-pF(t)+cmean.*(-1+pC(t)+pF(t))).^4.*(-1+20.*cmean.^3.*(-1+pC(t)+pF(t)).^3-10.*cmean.^2.*(-1+pC(t)+pF(t)).^2.*(1+6.*pF(t))-2.*pF(t).*(2+5.*pF(t).*(1+2.*pF(t)))+4.*cmean.*(-1+pC(t)+pF(t)).*(1+5.*pF(t).*(1+3.*pF(t)))));

        %  5 PBSL [3/-1]    gains-biased
        c(5,t)=1+cmean.^4.*(-1+pC(t)).^4.*(-35+21.*cmean.*(4+pC(t))+35.*cmean.^2.*(-2+(-2+pC(t)).*pC(t))+cmean.^3.*(20-3.*pC(t).*(-15+2.*pC(t).*(4+pC(t)))))+35.*(-1+cmean).*cmean.^2.*(-1+pC(t)).^2.*(3+cmean.*(-4.*(2+pC(t))+cmean.*(6-12.*(-2+pC(t)).*pC(t)+6.*cmean.*pC(t).*(-5+pC(t).*(2+pC(t)))+cmean.^2.*(-1+pC(t).*(10+3.*pC(t).*(1+(-4+pC(t)).*pC(t))))))).*pF(t)+21.*(-1+cmean).^2.*(-1+5.*cmean.*pC(t)+10.*cmean.^2.*(2+3.*(-2+pC(t)).*pC(t))-10.*cmean.^3.*(4-9.*pC(t)+4.*pC(t).^3)+5.*cmean.^4.*(5-2.*pC(t).^2.*(18+5.*(-4+pC(t)).*pC(t)))+cmean.^5.*(-4+5.*pC(t).*(-6+pC(t).*(24-20.*pC(t)+3.*pC(t).^3)))).*pF(t).^2-70.*(-1+cmean).^3.*(1-6.*cmean.*pC(t)-9.*cmean.^2.*(1+2.*(-2+pC(t)).*pC(t))+2.*cmean.^4.*(-1+pC(t)).*(2+5.*(-3+pC(t)).*pC(t).^2)+4.*cmean.^3.*(3-9.*pC(t)+5.*pC(t).^3)).*pF(t).^3-35.*(-1+cmean).^4.*(3+2.*cmean.*(-9.*pC(t)+cmean.*(-6-15.*(-2+pC(t)).*pC(t)+cmean.*(4-15.*pC(t)+10.*pC(t).^3)))).*pF(t).^4+21.*(-1+cmean).^5.*(-4+5.*cmean.*(cmean+(4+3.*cmean.*(-2+pC(t))).*pC(t))).*pF(t).^5+35.*(-1+cmean).^6.*(-1+3.*cmean.*pC(t)).*pF(t).^6-6.*(-1+cmean).^7.*pF(t).^7;
        
        %  6 PBSL [1/-3]    loss-averse
        c(6,t)=1+cmean.^4.*(-1+pC(t)).^2.*(-35.*(-1+pC(t)).^2+70.*cmean.^2.*(-1+4.*pC(t)-4.*pC(t).^3+pC(t).^4)+42.*cmean.*(2+pC(t).*(-6+pC(t)+3.*pC(t).^2))-cmean.^3.*(-20+pC(t).*(100+pC(t).*(10+3.*pC(t).*(-50+pC(t).*(25+2.*pC(t)))))))+35.*(-1+cmean).*cmean.^3.*(-1+pC(t)).^2.*(4-12.*cmean.*(-1+pC(t)).^2-4.*pC(t)+12.*cmean.^2.*(1+(-1+pC(t)).*pC(t).*(3+pC(t)))+cmean.^3.*(-4+pC(t).*(16+3.*pC(t).*(-4+(-4+pC(t)).*pC(t))))).*pF(t)+105.*(-1+cmean).^2.*cmean.^3.*(-1+pC(t)).*(-2+4.*cmean+2.*(5-4.*pC(t)).*pC(t)-2.*cmean.*pC(t).*(12+5.*(-3+pC(t)).*pC(t))+cmean.^2.*(-2+pC(t).*(14+pC(t).*(-22+pC(t).*(8+3.*pC(t)))))).*pF(t).^2-70.*(-1+cmean).^3.*cmean.^2.*(-1+pC(t)).*(3-3.*pC(t)+cmean.*(-4+10.*pC(t).*(-1+2.*pC(t))+cmean.*(1+5.*pC(t).*(3+2.*(-3+pC(t)).*pC(t))))).*pF(t).^3-35.*(-1+cmean).^4.*cmean.*(3-15.*cmean.*(-1+pC(t)).^2-3.*pC(t)+cmean.^2.*(11+5.*pC(t).*(-3+pC(t).*(-3+4.*pC(t))))).*pF(t).^4+21.*(-1+cmean).^5.*cmean.*(-8+14.*cmean+5.*(2+3.*cmean.*(-2+pC(t))).*pC(t)).*pF(t).^5+7.*(-1+cmean).^6.*(-1+3.*cmean.*(-3+5.*pC(t))).*pF(t).^6-6.*(-1+cmean).^7.*pF(t).^7;
        
        %  7 PBSL [1/0]     only gains
        c(7,t)=(1+cmean.*pC(t).*(7+cmean.*pC(t).*(-21+cmean.*pC(t).*(35+cmean.*pC(t).*(-35+cmean.*pC(t).*(21+cmean.*pC(t).*(-7+cmean.*pC(t)))))))-7.*(-1+cmean).*(-1+cmean.*pC(t)).^5.*(1+5.*cmean.*pC(t)).*pF(t)+21.*(-1+cmean).^2.*(-1+cmean.*pC(t)).^3.*(-1+cmean.*pC(t).*(2+9.*cmean.*pC(t))).*pF(t).^2-35.*(-1+cmean).^3.*(-1+cmean.*pC(t).*(8-12.*cmean.*pC(t)+5.*cmean.^3.*pC(t).^3)).*pF(t).^3-35.*(-1+cmean).^4.*(-1+cmean.*pC(t).*(9+5.*cmean.*pC(t).*(-3+cmean.*pC(t)))).*pF(t).^4+21.*(-1+cmean).^5.*(1+cmean.*pC(t).*(-8+9.*cmean.*pC(t))).*pF(t).^5-7.*(-1+cmean).^6.*(-1+5.*cmean.*pC(t)).*pF(t).^6+(-1+cmean).^7.*pF(t).^7)/2;
        
        %  8 PBSL [0/-1]    only losses
        c(8,t)=(2+cmean.^4.*(-70+168.*cmean+20.*cmean.^2.*(-7+2.*cmean)-140.*(-1+cmean).^3.*pC(t)+70.*(-1+cmean).^2.*(-2+5.*cmean).*pC(t).^3-70.*(-1+cmean).^2.*(-1+4.*cmean).*pC(t).^4-42.*(-1+cmean).*cmean.*pC(t).^5+28.*(-1+cmean).*cmean.^2.*pC(t).^6+cmean.^3.*pC(t).^7)+7.*(-1+cmean).*cmean.^3.*(-20.*(-1+pC(t)).^3+30.*cmean.*(-1+pC(t)).^2.*(-2+pC(t).*(4+pC(t)))+30.*cmean.^2.*(2+pC(t).*(-10+12.*pC(t)-5.*pC(t).^3+pC(t).^4))-cmean.^3.*(20+pC(t).*(-120+pC(t).*(150+pC(t).*(40+pC(t).*(-120+pC(t).*(24+5.*pC(t)))))))).*pF(t)+21.*(-1+cmean).^2.*cmean.^3.*pC(t).*(-30.*(-1+pC(t)).^2+20.*cmean.*(-1+pC(t)).*(-4+5.*pC(t))+cmean.^2.*(-50+pC(t).*(120+pC(t).*(-60+pC(t).*(-20+9.*pC(t)))))).*pF(t).^2-35.*(-1+cmean).^3.*cmean.^2.*(-6.*(-1+pC(t)).^2+4.*cmean.*(-1+pC(t)).*(-4+pC(t).*(2+5.*pC(t)))+cmean.^2.*(-10+pC(t).*(8+pC(t).*(36+5.*(-8+pC(t)).*pC(t))))).*pF(t).^3-35.*(-1+cmean).^4.*cmean.^2.*(-6+8.*cmean+pC(t).*(18-12.*pC(t)+cmean.*(-24+pC(t).*(12+5.*pC(t))))).*pF(t).^4+21.*(-1+cmean).^5.*cmean.*(2+pC(t).*(-2+cmean.*(-8+9.*pC(t)))).*pF(t).^5-7.*(-1+cmean).^6.*cmean.*(-4+5.*pC(t)).*pF(t).^6+(-1+cmean).^7.*pF(t).^7)/2;
        
        %  9 PBSL [3/1]     pbsl-conformist
        c(9,t)=cmean.^2.*(21.*(-1+cmean).^2.*cmean.^3.*pC(t).^5.*pF(t).^2+cmean.^2.*(35-2.*cmean.*(42+5.*cmean.*(-7+2.*cmean))+105.*(-1+cmean).^3.*pF(t)-21.*(-1+cmean).^2.*(-5+6.*cmean).*pF(t).^2+35.*(-1+cmean).^3.*pF(t).^3)+105.*(-1+cmean).^2.*cmean.^2.*pC(t).^4.*pF(t).*(-1+cmean+pF(t).*(5-6.*cmean+5.*(-1+cmean).*pF(t)))+105.*(-1+cmean).^2.*cmean.*pC(t).*((-1+pF(t)).^4-2.*cmean.*(-1+pF(t)).^2.*(1+(-4+pF(t)).*pF(t))+cmean.^2.*(1+(-1+pF(t)).*pF(t).*(8+(-7+pF(t)).*pF(t))))+35.*(-1+cmean).^2.*cmean.*pC(t).^3.*((-1+pF(t)).^2.*(1+5.*pF(t).*(-2+3.*pF(t)))-2.*cmean.*(1+3.*pF(t).*(-6+5.*(-2+pF(t)).^2.*pF(t)))+cmean.^2.*(1+pF(t).*(-24+5.*pF(t).*(18+pF(t).*(-16+3.*pF(t))))))+21.*(-1+cmean).^2.*pC(t).^2.*(-(-1+pF(t)).^5+cmean.*(-1+pF(t)).^3.*(8+pF(t).*(-31+3.*pF(t)))-cmean.^2.*(-1+pF(t)).*(13+pF(t).*(-112+pF(t).*(188+pF(t).*(-62+3.*pF(t)))))+cmean.^3.*(-6+pF(t).*(75+pF(t).*(-200+pF(t).*(150+(-30+pF(t)).*pF(t)))))));
        
        % 10 PBSL [-1/-3]   pbsl-contrarian, hates wins, hates losses even more
        c(10,t)=(-1+cmean).^2.*(1+cmean.*(2+cmean.*(3+cmean.*(4+cmean.*(-30+20.*cmean+7.*pC(t).^2.*(-30.*(-1+cmean)+pC(t).*(40.*(-1+cmean)+3.*(5+cmean.*(-5+pC(t))).*pC(t)))))-42.*cmean.^2.*pC(t).^2.*(15.*(-1+pC(t)).^2+cmean.*(-15+pC(t).*(30+(-15+pC(t)).*pC(t)))).*pF(t)+21.*cmean.*(10.*(-1+pC(t)).^3+10.*cmean.*(-1+pC(t)).*(-2+pC(t).*(4+5.*(-1+pC(t)).*pC(t)))+cmean.^2.*(-10+pC(t).*(30+pC(t).*(-60+pC(t).*(90+(-50+pC(t)).*pC(t)))))).*pF(t).^2+35.*(-1+cmean).*cmean.*(4.*(-1+pC(t)).^2.*(-2+5.*pC(t))+cmean.*(8+pC(t).*(-36+pC(t).*(54+5.*pC(t).*(-8+3.*pC(t)))))).*pF(t).^3+105.*(-1+cmean).^2.*cmean.*(-1+pC(t)).*(1+5.*(-1+pC(t)).*pC(t)).*pF(t).^4+21.*(-1+cmean).^3.*(-1+pC(t)).^2.*pF(t).^5)));
        
        % 11 PBSL a la McElreath
        c(11,t)=cmean.*(-3.*(-1+cmean).*cmean.*pC(t).^2.*pF(t)+cmean.*(3-2.*cmean+3.*(-1+cmean).*pF(t))-3.*(-1+cmean).^2.*pC(t).*(-1+pF(t).^2));
        
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