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
        cconf3=(3-2*cmean)*cmean^2;
        c(2,t)=cconf3;

        % 03 contrarian
        c(3,t)=(-1+cmean).^2.*(1+2.*cmean);

        %  4 PBSL [1/-1]    equal weights
        c(4,t)=-((1-pF(t)+cmean.*(-1+pC(t)+pF(t))).^2.*(-1-2.*pF(t)+2.*cmean.*(-1+pC(t)+pF(t))));   

        %  5 PBSL [3/-1]    gains-biased
        c(5,t)=1+cmean.^2.*(-1+pC(t)).^2.*(-3+cmean.*(2+pC(t)))-3.*(-1+cmean).*(-1+cmean.*(cmean+2.*(1+cmean.*(-2+pC(t))).*pC(t))).*pF(t)-3.*(-1+cmean).^2.*(-1+2.*cmean.*pC(t)).*pF(t).^2+(-1+cmean).^3.*pF(t).^3;        
        
        %  6 PBSL [1/-3]    loss-averse
        c(6,t)=1+cmean.^2.*(-1+pC(t)).*(3-3.*pC(t)+cmean.*(-2+pC(t).*(4+pC(t))))-6.*(-1+cmean).*cmean.*(1+cmean.*(-1+pC(t))).*(-1+pC(t)).*pF(t)-3.*(-1+cmean).^2.*cmean.*(-1+2.*pC(t)).*pF(t).^2+(-1+cmean).^3.*pF(t).^3;

        %  7 PBSL [1/0]     only gains
        c(7,t)=(1+cmean.*pC(t).*(3+cmean.*pC(t).*(-3+cmean.*pC(t)))-3.*(-1+cmean).*(-1+cmean.^2.*pC(t).^2).*pF(t)-3.*(-1+cmean).^2.*(-1+cmean.*pC(t)).*pF(t).^2+(-1+cmean).^3.*pF(t).^3)/2;

        %  8 PBSL [0/-1]    only losses
        c(8,t)=(2+6.*cmean.^2.*(-1+pC(t))+cmean.^3.*(4-6.*pC(t)+pC(t).^3)-3.*(-1+cmean).*cmean.*(2.*(-1+pC(t))+cmean.*(2+(-4+pC(t)).*pC(t))).*pF(t)-3.*(-1+cmean).^2.*cmean.*pC(t).*pF(t).^2+(-1+cmean).^3.*pF(t).^3)/2;

        %  9 PBSL [3/1]     pbsl-conformist
        c(9,t)=cmean.*(3.*(-1+cmean).*cmean.*pC(t).^2.*pF(t)+cmean.*(3-2.*cmean+3.*(-1+cmean).*pF(t))+3.*(-1+cmean).*pC(t).*(-1+cmean+(2+cmean.*(-4+pF(t))-pF(t)).*pF(t)));

        % 10 PBSL [-1/-3]   pbsl-contrarian, hates wins, hates losses even more
        c(10,t)=(-1+cmean).*(-1+cmean.*(-1+cmean.*(2-3.*pC(t).^2)+3.*cmean.*pC(t).^2.*pF(t)+3.*(-1+cmean).*(-1+pC(t)).*pF(t).^2));

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