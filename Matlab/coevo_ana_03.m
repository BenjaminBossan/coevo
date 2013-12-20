% environments can change within generations

function [x c w maxGen]=coevo_ana_03(...
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

% PBSLs 1
p31=(pC^3+3*pC^2*(1-pC));
p32=(3*pC^2*pF+3*pC^2*(1-pF)+6*pC*(1-pC)*(1-pF));
p33=(6*pC*pF*(1-pF)+3*pC*(1-pF)^2+3*(1-pC)*(1-pF)^2);
p34=((1-pF)^3+3*pF*(1-pF)^2);

% PBSLs 2
p41=(pC^3+3*pC^2*(1-pC)+3*pC*(1-pC)^2);
p42=p32;
p43=p33;
p44=(1-pF)^3;

% PBSLs 3
p51=1;
p52=(3*pC^2*pF+3*pC^2*(1-pF)+6*pC*(1-pC)*(1-pF)+3*(1-pC)^2*(1-pF));
p53=(6*pC*pF*(1-pF)+3*pC*(1-pF)^2);
p54=0;

x(:,1)=xInit;

for g=1:gen
    
    % FIRST PERIOD'S CHOICE OF THE STRATEGIES
    
    % 1st period's choice of conformists is right or wrong depending on the
    % probability of environmental change and the degree of conformism
    if g==1     % in very first generation
        c(:,1)=1/2;
    else
        c(:,1)=clast;
    end
    
    % individual learners always choose the correct option with prob. cind
    c(1,2:tmax)=cind;
    
    for t=2:tmax
        
        % in case of environmental change
        if rand<probEnvChange
            c(:,t-1)=1-c(:,t-1);
        end
        
        % mean proportion of correct choices
        cmean=x(1,g)*c(1,t-1)+x(2,g)*c(2,t-1)+x(3,g)*c(3,t-1)+...
            x(4,g)*c(4,t-1)+x(5,g)*c(5,t-1)+x(6,g)*c(6,t-1)+x(7,g)*c(7,t-1);
        
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
            (1-cmean)^3*        p44;
        
        % PBSL a la McElreath et al 2008
        c(5,t)=cmean^3*         p51+...
            cmean^2*(1-cmean)*3*p52+...
            cmean*(1-cmean)^2*3*p53+...
            (1-cmean)^3*        p54;
            
        % opportunistic individual learners
        c(6,t)=(cind-cind*pF+(3-2*cmean)*cmean^2*pF)/...
            (1+(cind+cmean^2*(-3+2*cmean))*(pC-pF));
        
        % opportunistic conformists
        if (pC==1&&pF==0)
            % this check must be included to avert dividing by 0
            c(7,t)=1;
        else
            c(7,t)=(-cmean^2*(-3+2*cmean)*(pF-1)-cind*pF)/...
                (-1+(cind+cmean^2*(-3+2*cmean))*(pC-pF));
        end
        
        % somehow matlab achieves to reach values that are slightly
        % above 1, which screws everything
        c(:,t)=min(c(:,t),ones(nStrats,1));
    end
    
    clast=c(:,tmax);
    
%     if length(find(c>1))|length(find(c<0))
%         'c error'
%         w(:,g)
%         x(:,g)
%         c
%         break
%     end
    
    % SELECTION
    
    % fitnesses
    w(:,g)=1+selCoef*(mean(c,2)-cind);
    % mean fitness
    wmean=sum(w(:,g).*x(:,g));
    % frequencies updated
    x(:,g+1)=x(:,g).*w(:,g)./wmean;
    
%     % eliminate rare strategies
%     x(:,g+1)=x(:,g+1).*(x(:,g+1)>1/1000000);
%     % normalize
%     x(:,g+1)=x(:,g+1)/sum(x(:,g+1));
    
    if sum(w(:,g))>=0
        1;
    else
        w(:,g)
        x(:,g)
        c
        break
    end
    
%     % EQUILIBRIUM CONDITION
%     % check every 1000th generation only (filter is so slow)
%     if mod(g,2500)==0
%         % never break during first 1 hundreth of simulation rounds
%         if g>max([windowSize*3 ceil(gen/100)])
% %             % moving average of frequency over generations
% %             xflat=filter(ones(1,windowSize)/windowSize,1,x(:,1:g)')';
% %             % first difference
% %             dxflat=diff(xflat')';
% %             % second difference
% %             ddxflat=abs(diff(dxflat')');
% %             % in case a strategy becomes fixed
%             if max(x(:,g+1))>.999;
%                 gen=min(g+genBonus,gen);
%             % in case frequencies move slowly and there is decelration
% %             elseif (max(abs(dxflat(:,g-1)))<stabCrit)&(max(abs(ddxflat(:,g-2)))<stabCrit)
% %                 gen=min(g+genBonus,gen);
%             end
%         end
%     end
    
end

maxGen=g;