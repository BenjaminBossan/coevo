function [n,pA,pB] = coevo_pbsl_02(...
        tmax,...    %tmax time of simulation
        nindi,...   %nindi # of individuals
        nstrat,...  %# of strategy types
        incr,...    %the increment
        pincr,...   %probability that patch quality changes
        w0,...      %the base fitness
        b,...       %the benefit from succeeding
        c,...       %the tallying cost as a function of the benefit
        cm,...      %the cost per memory slot
        ninitial,...%initial state of the population
        pA0,...     %remember last pA
        pB0,...     %remember last pB
        dp,...      %shift in mean pA and pB
        pfix,...    %a fixed environment for testing purposes
        param,...   %# of parameter values
        q,...       %oblivousness, discount for older memory
        lambda,...  %sensitivity factor for reinforcement learners; higher->steeper
        genetics,...%whether pure (0) or mixed strategies (1)
        kdoubt,...  %threshold value for ODCs
        compare_self,...    %ITW looks at own wealth?
        aversion,...%whether positive/negative info is taken into account
        regime);    %regime, how the environment changes. 1->low variance, 2->high variance


% INDICES
ichoice=1;isucc=2;ifit=3;itally=4;imemo=5;ibias=6;
itrack=7;irecent=8;iskill=9;istrat=10;ibest=11;isource=12;
igain=27;   % weight for gains for scoring type PBSLs
iloss=28;   % weight for losses for scoring type PBSLs
% STRATS:
ipind=13;   % 1  individual learner (threshold reinforcement learning)
ipmaj=14;   % 2  conformist
ipsuc=15;   % 3  payoff-biased social learner (PBSL) with equal weight gains/losses
ipwea=16;   % 4  imitate the wealthiest (ITW)
ipomt=17;   % 5  opportunistic conformist (OC)
ipoil=18;   % 6  opportunistic individual learner (OIL)
iprei=19;   % 7  probabilistic reinforcement learning
ipidc=20;   % 8  in doubt, conform (IDC)
ipomn=21;   % 9  omniscient
ipp31=22;   % 10 PBSL with weight on gains/losses of 3/1
ipp13=23;   % 11 PBSL with weight on gains/losses of 1/3
ippog=24;   % 12 PBSL who only factors in gains
ippco=25;   % 13 PBSL-conformist
ippmc=26;   % 14 PBSL à la McElreath

% INITIALIZE POPULATION TENSOR
n=zeros(nindi,param*tmax);
n=reshape(n,nindi,param,tmax);  %population tensor
n(:,:,1)=ninitial;              %the initial state of the population

n(:,itally,:)=repmat(ninitial(:,itally),1,tmax);    % the tally # stays the same within 1 gen
n(:,imemo,1)=0;                                   % at beginning, both equally likely
n(:,iskill,:)=repmat(ninitial(:,iskill),1,tmax);    % skill stays same
for s=1:nstrat
    n(:,ipind+s-1,:)=repmat(ninitial(:,ipind+s-1),1,tmax);
end
n(:,igain,:)=repmat(ninitial(:,igain),1,tmax);  % weight of gains stays the same throughout
n(:,iloss,:)=repmat(ninitial(:,iloss),1,tmax);  % weight of losses stays the same throughout

%TRACK whom wealth imitators copy and recentness of this information
n(:,itrack,:) =8*ones(nindi,tmax);  %8 codes for the non-wealth tallier
n(:,itrack,1) =9*ones(nindi,1);     %9 codes for initial guesses
n(:,irecent,1)=zeros(nindi,1);      %0 codes for non-wealth tallying


%create a RANDOM ENVIRONMENT according to specification
[pA,pB] = randomenvironment2(tmax,regime,incr,pincr,pA0,pB0);
% possible shift
pA=min(1-incr,max(incr,pA+dp));pB=min(1-incr,max(pB+dp,incr));
% %use this fixed environment only for testing purposes:
% pA=pfix(1,:);
% pB=pfix(2,:);
% 'achtung pfix'

% vector containing the SKILL
skillvec=ninitial(:,iskill);

%random matrix used for determining success
randsucmatA=rand(nindi,tmax);
randsucmatB=rand(nindi,tmax);

% determine which STRATEGY is used in EACH PERIOD,
% according to the PROBABILITY distribution
% n=coevo_determine_strat2(n,nindi,nstrat,istrat,ipind,tmax);
n(:,istrat,:)=3;

for t=1:tmax
    
    % S U C C E S S
    % update success
    n(:,isucc,t)=n(:,ichoice,t).*(randsucmatA(:,t)-skillvec<pA(t))+...
        mod(n(:,ichoice,t)+1,2).*(randsucmatB(:,t)-skillvec<pB(t));
    
    % F I T N E S S
    % the benefit vector
    % all individuals receive +1, except omniscients, who pay a cost
    % so as to equalize their fitness and the fitness of individual
    % learners (under standard conditions at least)
    bvec=b*ones(nindi,1)-0.055*(n(:,istrat,1)==9);
    if t==1     %in 1st round, add base fitness
        n(:,ifit,1)=w0*ones(nindi,1)+bvec.*n(:,isucc,t);
%         the first choice is random, therefore no costs
    else
        n(:,ifit,t)=n(:,ifit,t-1)...        % add fitness of last round
            +bvec.*n(:,isucc,t);                % add benefit if successful
    end
    
    % C H O I C E
    
    % the SCORE vector
    % this vector will determine which option an individual
    % chooses in the next round. A positive score leads to A
    % choice, a negative score to B choice.
    scvec=zeros(nindi,1);
    
    % matrix containing random individuals who are sampled
    randtmat=ceil(nindi*rand(nindi,max(n(:,itally,1))));
    
    % PAYOFF-BIASED SOCIAL LEARNING
    
    %maximum tally number
    tallymaxpbsl=max(n(n(:,istrat,t)==3,itally));
    
    % PBSLs with variable weights
    for ii=1:tallymaxpbsl
        scvec=scvec+(n(:,istrat,t)==3).*(...
            n(:,igain).*((n(randtmat(:,ii),ichoice,t)==1)&(n(randtmat(:,ii),isucc,t)==1))-... %choice A + suc
            n(:,igain).*((n(randtmat(:,ii),ichoice,t)==0)&(n(randtmat(:,ii),isucc,t)==1))+... %choice B + suc
            n(:,iloss).*((n(randtmat(:,ii),ichoice,t)==1)&(n(randtmat(:,ii),isucc,t)==0))-... %choice A + fail
            n(:,iloss).*((n(randtmat(:,ii),ichoice,t)==0)&(n(randtmat(:,ii),isucc,t)==0))); %choice B + fail
    end
        
    
    % choose according to the SCORE VECTOR
    
    % TIE BREAKING rule
    % (ties should not occur often if q=/=1 and if tally numbers of
    % success and majority talliers are odd)
    % in case of draw, stick to choice of last round
    repchoices=n(:,ichoice,t);
    scvec(find(scvec==0))=repchoices(find(scvec==0));
    
    %score>0 -> choose A, else choose B
    n(:,ichoice,t+1)=max(sign(scvec),0);
    
    % has the individual chosen the BETTER option?
    % if yes->1, no->-1, pA=pB->0
    n(:,ibest,t)=repmat(sign(pA(t)-pB(t)),nindi,1).*(2*n(:,ichoice,t)-1);
    
end
    
% in the last round, a choice is made for round tmax+1
% thus truncate very last round of n
n=n(:,:,1:tmax);

% the last entry of ibest, whether the better option was chosen
% is the sum of all entries, giving the aggregate performance
% on this level
n(:,ibest,tmax)=sum(n(:,ibest,:),3)/tmax;

% the last entry of isource, how often an individual was source for
% a wealth imitator, is the sum of all entries
n(:,isource,tmax)=sum(n(:,isource,:),3)/tmax;

% omniscient strategy has additional learning costs in terms of fitness,
% which will result in non-integer fitness values. Therefore, round the
% final fitness values.
n(:,ifit,tmax)=round(n(:,ifit,tmax));