clear all
% 'achtung: no clear all'

% in this instance, all strategies are idc, but with different thresholds

'no pbsl'

% P A R A M E T E R S

iterations=1;    %how often the simulation is repeated
gen=1;              %number of simulated generations
w0=10;              %base fitness
tmax=1000            %time per generation
pA0=1/2;            %initial success rate of option A
pB0=1/2;            %initial success rate of option B
pskill=1;           %max var of (pos. or neg.) influence of skill
nindi=10000;        %total population
regime=1;           %how the environment changes. 0->no regression to the mean, 1->medium, 2->high
incr=2/100;         %increment at which the environment becomes better or worse
pincr=1;            %probability that environmental quality changes at all after each period
param=22;           %number of parameters per individual
q=.9;               %discount rate for older memory
b=1;                %benefit from choosing correctly
mutationrate=0;     %the rate at which traits mutate
mutationincr=0;     %the factor by which the probability is modified
nstrat=10;          %number of strategies
lambda=1;           %sensitiviy factor for reinforcement learners
kdoubt=2;           %threshold value for use of maj. tal. by IDCs
compare_self=0;     % Usually, ITW only screens others and then adopts the
                    % choice of the best of them. However, they could
                    % include themselves in the sample, so that they don't
                    % switch if they are best themselves. To introduce this
                    % effect, set compare_self=1; else =0.
dpA=0;
dpB=0;          % changes in pA and pB;
pA0=pA0+dpA;pB0=pB0+dpB;

'random seed'
rand('seed',17413)

% STRATEGIES:
% each individual agent is characterized by the following parameters:
%   1) history of past choices (A->1, B->0)
%   2) history of past successes (1->success, 0->failure)
%   3) the fitness
%   4) NON-FUNCTIONAL the sample size
%   5) NON-FUNCTIONAL memory (for probabilistic reinforcement learning)
%   6) bias towards choosing A (for non-probabilistic reinforcement learning)
%   7) whom ITW tracks
%   8) recentness of information of ITW (age in periods)
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
% EXPLANATION: For the population tensor, the dimensions are:
%   1) individuals
%   2) parameters (choice, success...)
%   3) time
% INDICES
% for convenience, these indices are not passed to the coevo function.
% Changes here have thus to be updated in coevo as well
ichoice=1;isucc=2;ifit=3;itally=4;imemo=5;ibias=6;
itrack=7;irecent=8;iskill=9;istrat=10;ibest=11;isource=12;
% STRATS:
ipind=13;   % 1  individual learner (threshold reinforcement learning)
ipcon=14;   % 2  conformist, sample size 3
ipoil=15;   % 3  opportunistic individual learners, sample size 3
ipoc =16;   % 4  opportunistic conformists, sample size 3
ipidc=17;   % 5  in doubt, conform, sample size 3
ipitw=18;   % 6  imitate the wealthiest, sample size 7
ip4m1=19;   % 7  scoring-type PBSL weights [4/-1], sample size 7
ip10 =20;   % 8  scoring-type PBSL weights [1/0], sample size 3
ipMcE=21;   % 9  PBSLs McElreath, sample size 3
ippct=22;   % 10 PBSLs payoff-conformist trade-off, sample size 6


% I N I T I A L I Z A T I O N
ninitial=zeros(nindi,param);        %the initial state of the population
ntsize=nindi*param*gen;
if ntsize<=50000000
    nt=zeros(nindi,param*gen);          %tensor for just the final state of the population
    nt=reshape(nt,nindi,param,gen);
    if iterations<2
        perfmat=zeros(nstrat,gen);
    elseif gen<2
        perfmat=zeros(nstrat,iterations);
    end
%     fig5=zeros(1,gen);                  %proportion mean correct
else % if the tensor becomes too big, matlab cannot handle it -> track only every 10th generation
    nt=reshape(zeros(nindi,param*gen/10),nindi,param,gen/10);
    perfmat=zeros(nstrat,gen/10);
%     fig5=zeros(1,gen/10);
end
ntn=zeros(nindi,param);             %nt now

pfix=zeros(2,tmax);               %fixed environment, for testing purposes only
% pfix(1,:)=0.99*ones(1,tmax);
% pfix(1,1:tmax)=45/100;
% pfix(1,1:20)=55/100;pfix(1,41:60)=55/100;pfix(1,81:100)=55/100;
% pfix(1,1:10)=55/100;pfix(1,21:30)=55/100;pfix(1,41:50)=55/100;
% pfix(1,61:70)=55/100;pfix(1,81:90)=55/100;pfix(1,101:110)=55/100;
% pfix(1,1:30)=55/100;pfix(1,61:90)=55/100;
% pfix(1,1:100)=10/100;
% pfix(2,1:100)=90/100;
% pfix(1,101:end)=90/100;
% pfix(2,101:end)=10/100;
% pfix(1,[51:100,151:200,251:300])=0.01;
% pfix(2,:)=1-pfix(1,:);
% pfix=[pA;pB];


% INITIAL POPULATION
ninitial(:,ichoice)=rand(nindi,1)>1/2;                %initialize random choice in 1st round
ninitial(:,iskill)=(-.5+rand(nindi,1))*pskill;        %skill is uniformly distributed with mean 0

% ~~~~~~~~~~~~ change initial population here~~~~~~~~~~~~~~~
% use this procedure to define the initial population for
% pure strategies
ind=[...
    0;          % 0
    10000;          % 1
    0;          % 2
    0;          % 3
    0;          % 4
    0;          % 5
    0;          % 6
    0;          % 7
    0;          % 8
    0];         % 9

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
% ninitial(501:1000,itally)=6;

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
if gen>10||iterations>10
    wb=waitbar(0,'progress');
end

% S T A R T   O F   T H E   S I M U L A T I O N
countA=zeros(3,1);

for x=1:iterations

    for g=1:gen
        
        % initialize simulation
        [n,pA,pB,relyOnConf] = coevo_IDC_50(...
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
        lambda,...  %sensitivity factor for reinforcement learners, higher->steeper
        genetics,...%whether pure (0) or mixed strategies (1)
        kdoubt,...  %threshold value for IDCs
        compare_self,...    %ITW looks at own wealth?
        regime);    %regime, how the environment changes. 1->low variance, 2->high variance
    
        choicesvec=n(:,ichoice,tmax);   %remember as 1st choice for stats

        % I N H E R I T A N C E
        ninitial2=ninitial;                 %copy initial state
        
        if iterations<2
            if ntsize<=50000000
                nt(:,:,g)=n(:,:,tmax);              %the final state

                %determine performance
                pident=find(pA==pB);
                nperf=squeeze(n(:,ibest,:));
                nperf(:,pident)=[];
                nperf=(nperf+1)/2;
                [unik,unikchar]=coevo_find_unique2(n,istrat,nstrat);
                for ii=1:length(unik)
                    perfmat(unik(ii),g)=mean(mean(nperf(n(:,istrat,tmax)==unik(ii),:)));
                end

            else
                if mod(g,10)==0
                    nt(:,:,g/10)=n(:,:,tmax);
                    pident=find(pA==pB);
                    nperf=squeeze(n(:,ibest,:));
                    nperf(:,pident)=[];
                    nperf=(nperf+1)/2;
                    [unik,unikchar]=coevo_find_unique2(n,istrat,nstrat);
                    for ii=1:length(unik)
                        perfmat(unik(ii),g/10)=mean(mean(nperf(n(:,istrat,tmax)==unik(ii),:)));
                    end
                end
            end
        end
        ntn=n(:,:,tmax);
        
        % COPY INHERITABLE CHARACTERS + MUTATE
        ninitial=coevo_next_gen4(nindi,g,param,ntn,mutationrate,...
            nstrat,ichoice,ifit,itally,imemo,iskill,ipind,mutationincr);
        
        % include a check for EXTINCTION event
        if mod(g,50)==0
            [unik,unikchar]=coevo_find_unique2(n,istrat,nstrat);
            if length(unik)==1;
                'EXTINCTION'
                gen=g
                break
            end
        end

        % remember LAST STATE OF ENVIRONMENT
        pA0=pA(1,tmax);
        pB0=pB(1,tmax);
        
        % distribute skill anew
        ninitial(:,iskill)=(-.5+rand(nindi,1))*pskill;
        
        %progression
        if gen>10||iterations>10
            if iterations==1
                if mod(g,10)==0
                    waitbar(g/gen,wb);
                end
            elseif iterations>1
                if mod(x,10)==0
                    waitbar(x/iterations,wb);
                end
            end
        end

    end
    
    % initiation of measure variables
    if ((g==1)&(x==1))
        [unik,unikchar]=coevo_find_unique2(n,istrat,nstrat);
        % the behavior of the strategies
        superfig=zeros(length(unik),tmax*iterations);
        % the environment
        superp=zeros(1,tmax*iterations);
%         end
    end
    
    % determine performance
    if gen<2
        nt(:,:,g)=n(:,:,tmax);              %the final state
        %determine performance
        pident=find(pA==pB);
        nperf=squeeze(n(:,ibest,:));
        nperf(:,pident)=[];
        nperf=(nperf+1)/2;
        [unik,unikchar]=coevo_find_unique2(n,istrat,nstrat);
        for ii=1:length(unik)
            perfmat(unik(ii),x)=mean(mean(nperf(n(:,istrat,tmax)==unik(ii),:)));
        end
    end
    
    % how often has this strategy been the role model for ITW?
    for s=1:length(unik)
        sourcestat(s,x)=mean(n(find(n(:,istrat,1)==unik(s)),isource,tmax));
    end
    fig1=zeros(length(unik),tmax);
    for t=1:tmax
        for j=1:length(unik)
            ind0=find((n(:,ichoice,t)==1)&(n(:,istrat,t)==unik(j)));
            ind1=find(n(:,istrat,t)==unik(j));
            fig1(j,t)=length(ind0)/length(ind1);
        end
    end
    superfig(:,((x-1)*tmax+1):x*tmax)=fig1;
    superp(((x-1)*tmax+1):x*tmax)=pA-pB;
    
    % for more than 1 iteration: use original initial state
    % as input for the next generation, except that choices
    % are inherited
    if iterations>1
        ninitial=ninitial2;
        ninitial(:,ichoice)=choicesvec;
        % distribute skill anew
        ninitial(:,iskill)=(-.5+rand(nindi,1))*pskill;
    end
    
end

% E N D   O F   S I M U L A T I O N

% in case memory would be exceeded
if ntsize>50000000
    gen=gen/10;
end

% smoothing performance
perfsmooth=perfmat;
windowSize=ceil(gen/10);
for j=1:10
    perfsmooth(j,:)=filter(ones(1,windowSize)/windowSize,1,perfmat(j,:));
end


% R E S U L T S

if iterations==1
    
    if gen==1
        fig2=coevo_plot_1gen4(n,nt,nindi,tmax,nstrat,genetics,pA,pB,...
            istrat,ichoice,itrack,irecent,ifit);
            '    strat     performance'
            [unik mean(perfmat(unik,:),2)]
    end
    
    if gen>1
        fig3=coevo_plot_gens_4(n,nt,nindi,tmax,nstrat,genetics,gen,...
            istrat,ichoice,ipind,ntsize);
    end

elseif iterations>1
    % statistics
    % frequency of the different strategies
    stratfreq=zeros(1,length(unik));
    for s=1:length(unik); stratfreq(s)=sum(n(:,istrat)==unik(s));end
    stratfreq=stratfreq/sum(stratfreq);
    % normalize sourcestat
    if length(find(unik==4))==1
        sourcestat=sourcestat/stratfreq(find(unik==4));
    end
    % create cell that contains summary statistics
    cellu=cell(2+length(unik),3);
    cellu(1,1)=cellstr('Strategy');
    cellu(1,2)=cellstr('performance');
    cellu(1,3)=cellstr('SEM');
    cellu(2:1+length(unik),1)=unikchar;
    cellu(2+length(unik),1)=cellstr('group mean');
    cellu(2:1+length(unik),2)=num2cell(mean(perfmat(unik,:),2));
    cellu(2+length(unik),2)=num2cell(sum(stratfreq.*mean(perfmat(unik,:),2)'));
    cellu(2:1+length(unik),3)=num2cell(1.96*std(perfmat(unik,:)')'/sqrt(iterations));
    cellu(2+length(unik),3)=num2cell(sum(stratfreq.*(1.96*std(perfmat(unik,:)')'/sqrt(iterations))'));
    stats=cellu
end

toc
if gen>10||iterations>10
    close(wb) %close waitbar
end

