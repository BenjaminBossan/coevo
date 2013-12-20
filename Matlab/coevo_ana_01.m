function [x c w maxGen]=coevo_ana_01(...
    xInit,...                   % initial conditions
    cind,...                    % performance individual learners
    tmax,...                    % periods per generation
    gen,...                     % generations
    selCoef,...                 % maximum selection coefficient
    wind,...                    % fitness individual learners
    probEnvChange,...           % probability of environmental change
    pC,pF,...                   % probabilities that correct/false options yield success
    nStrats,...                 % total possible number of strategies
    windowSize,...              % window size for moving average filter
    stabCrit,...                % stability criterion to estimate equilibria
    genBonus);                  % # of generations that simulation continues to
                                % run after stability criterion is reached

% INITIALIZATION
x=zeros(nStrats,gen);       % frequency matrix
c=zeros(nStrats,tmax);      % percent correct choices matrix
w=zeros(nStrats,gen);       % fitness matrix

% certain fixed probabilities that are used repeatedly:
p31=(pC^3+3*pC^2*(1-pC));
p32=(pC^2+pC*(1-pC)*(1-pF));
p33=(pC*(1-pF)+(1-pC)*(1-pF)^2);
p34=((1-pF)^3+3*pF*(1-pF)^2);

p41=(pC^3+3*pC*(1-pC));
p42=(pC^2+pC*(1-pC)*(1-pF));
p43=(pC*(1-pF)+(1-pC)*(1-pF)^2);
p44=(1-pF)^3;
        
p52=(pC^2+(1-pC)*(1-pF));
p53=pC*(1-pF);

x(:,1)=xInit;

for g=1:gen
    
    % FIRST PERIOD'S CHOICE OF THE STRATEGIES
    
    % 1st period's choice of conformists is right or wrong depending on the
    % probability of environmental change and the degree of conformism
    if g==1     % in very first generation
        c(:,1)=1/2;
    else
        % if no environmental change
        if rand>probEnvChange;
            c(:,1)=clast;
        % if environmental change
        else
            c(:,1)=1-clast;
        end
    end
    
    % individual learners always choose the correct option with prob. cind
    c(1,2:tmax)=cind;
    
    for t=2:tmax
        
        % mean proportion of correct choices
        
%         cmean=sum(x(:,g).*c(:,t-1));
        cmean=x(1,g)*c(1,t-1)+x(2,g)*c(2,t-1)+x(3,g)*c(3,t-1)+...
            x(4,g)*c(4,t-1)+x(5,g)*c(5,t-1);
        
        % conformists
        c(2,t)=(3-2*cmean)*cmean^2;

        % PBSL equal weight
        c(3,t)=cmean^3*         p31+...
            cmean^2*(1-cmean)*3*p32+...
            cmean*(1-cmean)^2*3*p33+...
            (1-cmean)^3*        p34;
        
        % PBSL more weight gains
        c(4,t)=cmean^3*         p41+...
            cmean^2*(1-cmean)*3*p42+...
            cmean*(1-cmean)^2*3*p43+...
            (1-cmean)^3*        p44^3;
        
        % PBSL ŕ la McElreath et al 2008
        c(5,t)=cmean^3+...
            cmean^2*(1-cmean)*3*p52+...
            cmean*(1-cmean)^2*3*p53;
            
        
    end
    
    clast=c(:,tmax);
    
    % SELECTION
    
    % fitnesses
    w(:,g)=1+selCoef*(mean(c,2)-cind);
    % mean fitness
    wmean=sum(w(:,g).*x(:,g));
    % frequencies updated
    x(:,g+1)=x(:,g).*w(:,g)./wmean;
    
    % EQUILIBRIUM CONDITION
    % check every 1000th generation only (filter is so slow)
    if mod(g,2500)==0
        % never break during first 1 hundreth of simulation rounds
        if g>max([windowSize*3 ceil(gen/100)])
%             % moving average of frequency over generations
%             xflat=filter(ones(1,windowSize)/windowSize,1,x(:,1:g)')';
%             % first difference
%             dxflat=diff(xflat')';
%             % second difference
%             ddxflat=abs(diff(dxflat')');
%             % in case a strategy becomes fixed
            if max(x(:,g+1))>.999;
                gen=min(g+genBonus,gen);
            % in case frequencies move slowly and there is decelration
%             elseif (max(abs(dxflat(:,g-1)))<stabCrit)&(max(abs(ddxflat(:,g-2)))<stabCrit)
%                 gen=min(g+genBonus,gen);
            end
        end
    end
    
end

maxGen=g;