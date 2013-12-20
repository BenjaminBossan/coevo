function [n,pA,pB,relyOnConf] = coevo_IDC_50(...
        tmax,...    %tmax time of simulation
        nindi,...   %nindi # of individuals
        nstrat,...  %# of strategy types
        incr,...    %the increment
        pincr,...   %probability that patch quality changes
        w0,...      %the base fitness
        b,...       %the benefit from succeeding
        ninitial,...%initial state of the population
        pA0,...     %remember last pA
        pB0,...     %remember last pB
        dpA,...     %change in mean pA
        dpB,...     %change in mean pB
        pfix,...    %a fixed environment for testing purposes
        param,...   %# of parameter values
        q,...       %oblivousness, discount for older memory
        lambda,...  %sensitivity factor for reinforcement learners; higher->steeper
        genetics,...%whether pure (0) or mixed strategies (1)
        kdoubt,...  %threshold value for ODCs
        compare_self,...    %ITW looks at own wealth?
        regime);    %regime, how the environment changes. 1->low variance, 2->high variance


% INDICES
ichoice=1;isucc=2;ifit=3;itally=4;imemo=5;ibias=6;
itrack=7;irecent=8;iskill=9;istrat=10;ibest=11;isource=12;
% STRATS:
ipind=13;   % 1  individual learner (threshold reinforcement learning)
ipcon=14;   % 2  conformist, sample size 3
ipoil=15;   % 3  opportunistic individual learners, sample size 3
ipoc= 16;   % 4  opportunistic conformists, sample size 3
ipidc=17;   % 5  in doubt, conform, sample size 3
ipitw=18;   % 6  imitate the wealthiest, sample size 7
ip4m1=19;   % 7  scoring-type PBSL weights [4/-1], sample size 7
ip10= 20;   % 8  scoring-type PBSL weights [1/0], sample size 3
ipMcE=21;   % 9  PBSLs McElreath, sample size 3
ippct=22;   % 10 PBSLs payoff-conformist trade-off, sample size 6

% INITIALIZE POPULATION TENSOR
n=zeros(nindi,param*tmax);
n=reshape(n,nindi,param,tmax);  %population tensor
n(:,:,1)=ninitial;              %the initial state of the population

n(:,iskill,:)=repmat(ninitial(:,iskill),1,tmax);    % skill stays same

for s=1:nstrat
    n(:,ipind+s-1,:)=repmat(ninitial(:,ipind+s-1),1,tmax);
end

%TRACK whom wealth imitators copy and recentness of this information
n(:,itrack,:) =8*ones(nindi,tmax);  %8 codes for the non-wealth tallier
n(:,itrack,1) =9*ones(nindi,1);     %9 codes for initial guesses
n(:,irecent,1)=zeros(nindi,1);      %0 codes for non-wealth tallying

% change in mean pA and pB:
% first subtract the change to get normal behavior of pA and pB
pA0=pA0-dpA;pB0=pB0-dpB;

%create a RANDOM ENVIRONMENT according to specification
if regime==12345
    [pA pB] = randomenvironment_AR_01(tmax,pA0,pB0,0,0.9925,0.03);
else
    [pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);
end
%use this fixed environment only for testing purposes:
% pA=pfix(1,:);
% pB=pfix(2,:);
% 'achtung pfix'

% apply changes in mean pA and pB, consider the boundaries
pA=min(1,max(0,pA+dpA));pB=min(1,max(pB+dpB,0));

% vector containing the SKILL
skillvec=ninitial(:,iskill);

%random matrix used for determining success
randsucmatA=rand(nindi,tmax);
randsucmatB=rand(nindi,tmax);

% determine which STRATEGY is used in EACH PERIOD,
% according to the PROBABILITY distribution
if sum(mod(unique(n(:,ipind:ipind+nstrat-1,1)),1))==0
    %if there are only pure strategies
    for jj=1:nindi
        n(jj,istrat,1:tmax)=find(n(jj,ipind:ipind+nstrat-1,1)==1);
    end
else
    %if there are mixed strategies
    n=coevo_determine_strat2(n,nindi,nstrat,istrat,ipind,tmax);
end

% how often a strat relies on conformism
relyOnConf = zeros(1,10);

for t=1:tmax
    
    % S U C C E S S
    % update success
    n(:,isucc,t)=n(:,ichoice,t).*(randsucmatA(:,t)-skillvec<pA(t))+...
        mod(n(:,ichoice,t)+1,2).*(randsucmatB(:,t)-skillvec<pB(t));
    
    % CHOICE FREQUENCIES
    % frequency of A choices in last period:
    x=sum(n(:,ichoice,t))/nindi;
    
    
    % update BIAS towards A or B for strategies relying on individual
    % learning
    n(:,ibias,t+1)=q*n(:,ibias,t)+...   %discount older bias
        ((n(:,ichoice,t)==n(:,isucc,t))*2-1);
    
    % bias+1 if A succeeded or B failed, -1 vice versa
    
    % F I T N E S S
    % the benefit vector
    % all successful individuals receive +1
    if t==1
        n(:,ifit,t)=w0+b*n(:,isucc,t);
    else
        n(:,ifit,t)=n(:,ifit,t-1)...    % add fitness of last round
            +b*n(:,isucc,t);            % add benefit if successful
    end
    
    % C H O I C E
    
    % the SCORE vector
    % this vector will determine which option an individual
    % chooses in the next round. A positive score leads to A
    % choice, a negative score to B choice.
    scvec=zeros(nindi,1);
    
    
    % KDOUBT VECTOR
    kdoubt = 100*ones(1,10);
    kdoubt(6:end) = 0;
%     kdoubt(4:7) = 2;
    
    % random vector with size of population
    randIndi=rand(nindi,1);
    
    % IN DOUBT, CONFORM
    
    % choice made by a conformist
    bConf = 2*(randIndi<((3-2*x)*x^2))-1;
    
    % choice made by an IL
    bIndi = n(:,ibias,t+1);
    
    % if there is doubt, use conformism
    for kd = 1:10
        
        % index of individuals that are uncertain
        indDoubt = (abs(n(:,ibias,t+1))<kdoubt(kd));
        
        % strategies that are updated now
        indexIDC = (n(:,istrat,t)==kd);
        
        % strategies that use conformism
        index1 = indexIDC.*indDoubt;
        
        % if uncertain, use conformism
        scvec=scvec+index1.*bConf;
        
        % if certain, use individual learning
        scvec=scvec+(1-index1).*indexIDC.*bIndi;
        
        % remember how often a strat relies on conformism
        if sum(indexIDC) > 0
            relyOnConf(kd) = relyOnConf(kd) + sum(index1)/sum(indexIDC)/tmax;
        else
            relyOnConf(kd) = 0;
        end
        
    end
    
    % choose according to the SCORE VECTOR
    
    % TIE BREAKING rule
    % in case of draw, random choice
    
    scvec(find(scvec==0))=((rand(1,length(find(scvec==0))))>.5)-.5;
    
    %score>0 -> choose A, else choose B
    n(:,ichoice,t+1)=max(sign(scvec),0);
    
    % has the individual chosen the BETTER option?
    % if yes->1, no->-1, pA=pB->0
    dpApB=sign(pA(t)-pB(t));
    n(:,ibest,t)=((2*n(:,ichoice,t)-1)*dpApB);
    
end
    
% in the last round, a choice is made for round tmax+1
% thus truncate very last round of n
n=n(:,:,1:tmax);

% the last entry of isource, how often an individual was source for
% a wealth imitator, is the sum of all entries
n(:,isource,tmax)=sum(n(:,isource,:),3)/tmax;

% omniscient strategy has additional learning costs in terms of fitness,
% which will result in non-integer fitness values. Therefore, round the
% final fitness values.
n(:,ifit,tmax)=round(n(:,ifit,tmax));