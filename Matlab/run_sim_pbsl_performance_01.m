% in this specific instance, different weights are tested for their
% performance, using many iterations and testing all the different weight
% combinations. It is assumed that populations are homogeneous.
% No evolution occurs

clear all



% P A R A M E T E R S

iterations=1000;    %how often the simulation is repeated
gen=1;              %number of simulated generations
w0=10;              %base fitness
tmax=50;            %time per generation
pA0=1/2;            %initial success rate of option A
pB0=1/2;            %initial success rate of option B
dp=.25               %shift in mean pA and pB
pskill=0;           %max var of (pos. or neg.) influence of skill
nindi=1000        %total population
regime=1;           %how the environment changes. 0->no regression to the mean, 1->medium, 2->high
incr=2/100;         %increment at which the environment becomes better or worse
pincr=1;            %probability that environmental quality changes at all after each period
param=28;           %number of parameters per individual
q=.9;               %discount rate for older memory
b=1;                %benefit from choosing correctly
c=0;cm=0;           %cost of learning individually or socially; not functional
mutationrate=0.01;  %the rate at which traits mutate
mutationincr=0;     %the factor by which the probability is modified
monitorCrosscorr=1; %if this option is 1, the crosscorrelation betw. envir. and behavior is monitored
nstrat=14;          %number of strategies
lambda=1;           %sensitiviy factor for reinforcement learners
kdoubt=2;           %threshold value for use of maj. tal. by IDCs
compare_self=0;     % Usually, ITW only screens others and then adopts the
                    % choice of the best of them. However, they could
                    % include themselves in the sample, so that they don't
                    % switch if they are best themselves. To introduce this
                    % effect, set compare_self=1; else =0.
aversion=0;         %gives the weight of losses/gains, ie, =1->equal, >1 loss aversion, <1 win affinity
                    %aversion=0 is the same as =1; does not work anymore
                    %for PBSLs
tallyn=7           %tally number, the number of "screened" individuals
weightscope=5;     %absolute number of maximum weight for gains and losses of PBSLs

% 'random seed'
% rand('seed',17415)

% STRATEGIES:
% each individual agent is characterized by the following parameters:
%   1) history of past choices (A->1, B->0)
%   2) history of past successes (1->success, 0->failure)
%   3) the fitness
%   4) the tally number (how many individuals are screened)
%   5) memory (for probabilistic reinforcement learning)
%   6) bias towards choosing A (for non-probabilistic reinforcement learning)
%   7) whom wealth imitators track
%   8) recentness of information of wealth imitators (age in periods)
%   9) individual differences in skill
%  10) which strategy is chosen this round (important for mixed strategies)
%  11) whether in this period, the individual has chosen the better of the
%      two options (yes->1, no->0, draw->0.5)
%  12) whenever an individual is imitated by ITW, isource gets +1
%  --- an individual has a certain probability to use each strategy
%      usually, each individual has only one strategy that she uses
%      each round, but mixed strategies are possible. In that case
%      the following parameters determine the probabilities. In case
%      of pure strategies, the probability for one strategy is 1 and
%      for the rest 0.
%  13) probability to use individual learning (deterministic reinf. learning)
%  14) probability to use majority tallying (conformity)
%  15) probability to use success tallying (payoff-biased)
%  16) probability to use wealth imitation (ITW, prestige)
%  17) probability to use hierarchical strategy 1: opportunistic majority
%      tallying (OMT)
%      this strat uses ind. learning if successful, else maj. tallying
%  18) probability to use hierarchical strategy 2: opportunistic individual
%      learning (OIL)
%      this strat uses maj. tallying if successful, else ind. learning
%  19) probability to use probabilistic reinforcement learning
%  20) probability to use hierarchical strategy 3: In Doubt, Conform (IDC)
%      this strategy uses individual learning, but if the bias towards A or
%      B is too small (<kdoubt), it uses majority tallying
%  21) probability to use omniscient strategy, a strategy that immediately
%      knows what to choose
% EXPLANATION: For the population tensor, the dimensions are:
%   1) individuals
%   2) parameters (choice, success...)
%   3) time
% INDICES
% for convenience, these indices are not passed to the coevo function.
% Changes here have thus to be updated in coevo as well
ichoice=1;isucc=2;ifit=3;itally=4;imemo=5;ibias=6;
itrack=7;irecent=8;iskill=9;istrat=10;ibest=11;isource=12;
igain=27;   % weight for gains for scoring type PBSLs
iloss=28;   % weight for losses for scoring type PBSLs
% STRATS:
ipind=13;   % 1  individual learner (threshold reinforcement learning)
ipmaj=14;   % 2  conformist
ipsuc=15;   % 3  payoff-biased social learner (PBSL) with variable weight gains/losses
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


% I N I T I A L I Z A T I O N
ninitial=zeros(nindi,param);        %the initial state of the population
ntsize=nindi*param*gen;
if ntsize<=50000000
    nt=zeros(nindi,param*gen);          %tensor for just the final state of the population
    nt=reshape(nt,nindi,param,gen);
    fig5=zeros(1,gen);                  %proportion mean correct
else % if the tensor becomes too big, matlab cannot handle it -> track only every 10th generation
    nt=reshape(zeros(nindi,param*gen/10),nindi,param,gen/10);
    fig5=zeros(1,gen/10);
end
ntn=zeros(nindi,param);             %nt now

pfix=.5*ones(2,tmax);               %fixed environment, for testing purposes only
pfix(1,1:tmax)=45/100;
pfix(1,1:20)=55/100;pfix(1,41:60)=55/100;pfix(1,81:100)=55/100;
% pfix(1,1:10)=55/100;pfix(1,21:30)=55/100;pfix(1,41:50)=55/100;
% pfix(1,61:70)=55/100;pfix(1,81:90)=55/100;pfix(1,101:110)=55/100;
% pfix(1,1:30)=55/100;pfix(1,61:90)=55/100;
% pfix(1,1:60)=55/100;
pfix(2,1:tmax)=50/100;
% pfix(1,151:200)=4.5/10;
% pfix(1,201:300)=3/10;
% pfix(1,251:300)=3.5/10;
% pfix(1,:)=5/10;
% pfix=[pA;pB];


% INITIAL POPULATION
ninitial(:,ichoice)=rand(nindi,1)>1/2;                %initialize random choice in 1st round
ninitial(:,iskill)=(-.5+rand(nindi,1))*pskill;        %skill is uniformly distributed with mean 0
ninitial(:,itally)=tallyn;                            %tally#
ninitial(:,[igain iloss])=repmat([1 -1],nindi,1);     %weight gains and losses

% ~~~~~~~~~~~~ change initial population here~~~~~~~~~~~~~~~
% use this procedure to define the initial population for
% pure strategies
ind=[...
    0;          % INDividual learners
    0;          % MAJority talliers (conformists)
    1000;          % SUCcess talliers (PBSL 1/1)
    0;          % WEAlth imitators (ITW)
    0;          % OCs
    0;          % OILs
    0;          % prob. REInforcement learners
    0;          % IDC
    0;          % OMNiscient strategy
    0;          % PBSL gains/losses 3/1
    0;          % PBSL gains/losses 1/3
    0;          % PBSL only gains
    0;          % PBSL-conformist
    0];         % PBSL McElreath

cc=1;cind=cumsum(ind);pop=zeros(nindi,nstrat);
if cind(length(cind))~=nindi;'error in population length'
    break;end
for i=1:nstrat
    pop([cc:cind(i)],i)=1;
    cc=cc+ind(i);
end
ninitial(:,ipind:ipind+nstrat-1)=pop;
clear cc;clear pop;clear ind;clear cind;

% 'achtung2'
% ninitial(:,[ipind,ipmaj])=repmat([.7 .3],nindi,1);
% ninitial(001:10,itally)=6;

% use this procedure for mixed strategies
% stratmat=zeros(nindi,nstrat);
% for i=1:nindi;stratmat(i,:)=randperm(nstrat);end
% stratmat=stratmat/10;
% ninitial(:,[ipind,ipmaj,ipsuc,ipwea])=stratmat;


% ninitial(:,[ipind,ipwea])=repmat([.01 .99],nindi,1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % manipulate these values to influence the initial choice
% ninitial(ninitial(:,ipind)==1,ichoice)=...
%     rand(sum(ninitial(:,ipind)==1),1)<3/5;
% ninitial(ninitial(:,ipsuc)==1,ichoice)=...
%     rand(sum(ninitial(:,ipsuc)==1),1)<.1;
% 'achtung'

% ERROR MESSAGES
% costs not implemented yet
if c>0||cm>0
    'error: costs not implemented in fitness calculation'
    break
end
if iterations>1&&gen>1
    'error: statistics only available for 1 generation right now'
    break
end

% check whether there are only pure strategies
% if yes -> genetics = 0
% if not -> genetics = 1
% as soon as mutations are possible, genetics are not pure
genetics=coevo_check_genetics2(...
    mutationrate,nstrat,ninitial,ipind,nindi,gen);

% measure duration of simulation
tic

% progress bar
wb=waitbar(0,'progress');

% monitored performance
% performance according to weight
perfmat=zeros(2*weightscope+1,(2*weightscope+1)*iterations);
perfmat=reshape(perfmat,2*weightscope+1,2*weightscope+1,iterations);
% best choices according to weight
bestmat=perfmat;
% fitness according to weight
fitmat=perfmat;

% S T A R T   O F   T H E   S I M U L A T I O N

choicesvec=rand(nindi,1)>1/2; 

for x=1:iterations
    
    % weight of gains
    for w1=1:2*weightscope+1
        
        %weight of losses
        for w2=1:2*weightscope+1
        
        % weight of gains and losses (some combinations will be redundant)
        ninitial(:,[igain iloss])=repmat([w1-weightscope-1 w2-weightscope-1],nindi,1);
        
        % take over choices
        ninitial(:,ichoice)=choicesvec;
        
        % initialize simulation
        [n,pA,pB] = coevo_pbsl_02(...
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
        lambda,...  %sensitivity factor for reinforcement learners, higher->steeper
        genetics,...%whether pure (0) or mixed strategies (1)
        kdoubt,...  %threshold value for IDCs
        compare_self,...    %ITW looks at own wealth?
        aversion,...%whether positive/negative information is taken into account
        regime);    %regime, how the environment changes. 1->low variance, 2->high variance
    
        % remember LAST STATE OF ENVIRONMENT
        pA0=pA(1,tmax);
        pB0=pB(1,tmax);
        
        choicesvec=n(:,ichoice,tmax);   %remember as 1st choice for stats
        
        % monitor performance
        [meancorr,unik,unikchar,majority_was_right,corresp]=...
            coevo_meancorrect8(pA,pB,n,tmax,nindi,ichoice,...
            istrat,ibest,nstrat,genetics);
        perfmat(w1,w2,x)=meancorr;
        bestmat(w1,w2,x)=mean(n(:,ibest,tmax));
        fitmat(w1,w2,x)=mean(n(:,ifit,tmax))-w0;
        
        end
    end
    
    waitbar(x/iterations,wb)
end

toc
close(wb) %close waitbar


