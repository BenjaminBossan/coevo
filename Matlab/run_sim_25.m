clear all



% P A R A M E T E R S

iterations=1000;    %how often the simulation is repeated
gen=1;              %number of simulated generations
w0=10;              %base fitness
tmax=50;             %time per generation
pA0=1/2;            %initial success rate of option A
pB0=1/2;            %initial success rate of option B
pskill=0;           %max var of (pos. or neg.) influence of skill
nindi=1000;        %total population
regime=1;           %how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;         %increment at which the environment becomes better or worse
pincr=1;            %probability that environmental quality changes at all after each period
param=20;           %number of parameters per individual
q=.9;               %discount rate for older memory
b=1;                %benefit from choosing correctly
c=0;cm=0;           %cost per tally #, memory slot respectively
mutationrate=0;     %the rate at which traits mutate
mutationincr=0;     %the factor by which the probability is modified
monitorCrosscorr=1; %if this option is 1, the crosscorrelation betw. envir. and behavior is monitored
nstrat=8;           %number of strategies
lambda=1;           %sensitiviy factor for reinforcement learners
kdoubt=2;           %threshold value for use of maj. tal. by ODCs
aversion=0;         %positive (+) or negative (-) information is taken into account more
                    %strongly, with factor abs(aversion); 0,-1,+1->unbiased
hyperb=-.02          %hyperbolic discounting: the "bonus" (>0) for most recent experience
                    %if==0 -> no hyperbolic discounting
                    %0; if<0 most recent memory counts less
tallyn=7;           %tally number, the number of "screened" individuals

% 'random seed'
% rand('seed',17414)

% STRATEGIES:
% each individual agent is characterized by the following parameters:
%   1) history of past choices
%   2) history of past successes
%   3) the fitness
%   4) the tally number / memory (the inherited strategy)
%   5) memory
%   6) bias towards choosing A
%   7) whom wealth talliers track         
%   8) recentness of information
%   9) individual differences in skill
%  10) which strategy is chosen this round
%  11) whether in this period, the individual has chosen the better of the
%      two options; yes->1, no->0, draw->0.5
%  12) whenever an individual is imitated by ITW, isource gets +1
%  13) probability to use individual learning
%  14) probability to use majority tallying
%  15) probability to use success tallying
%  16) probability to use wealth imitation
%  17) probability to use hierarchical strategy 1
%      this strat uses ind. learning if successful, else maj. tallying
%  18) probability to use hierarchical strategy 2
%      this strat uses maj. tallying if successful, else ind. learning
%  19) probability to use reinforcement learning
%  20) probability to use "When in doubt, copy"
% EXPLANATION: For the population tensor, the dimensions are:
%   1) individuals
%   2) parameters (choice, success...)
%   3) time
%   strategies are encoded: 1=indi, 2=majo, 3=success, 4=wealth
% INDICES
ichoice=1;isucc=2;ifit=3;itally=4;imemo=5;ibias=6;
itrack=7;irecent=8;iskill=9;istrat=10;
ibest=11;isource=12;
ipind=13;ipmaj=14;ipsuc=15;ipwea=16;iphi1=17;iphi2=18;iprei=19;
ipdou=20;


% I N I T I A L I Z A T I O N
ninitial=zeros(nindi,param);        %the initial state of the population
ntsize=nindi*param*gen;
if ntsize<=50000000
    nt=zeros(nindi,param*gen);          %tensor for just the final state of the population
    nt=reshape(nt,nindi,param,gen);
    fig5=zeros(1,gen);                  %proportion mean correct
else
    nt=reshape(zeros(nindi,param*gen/10),nindi,param,gen/10);
    fig5=zeros(1,gen/10);
end
ntn=zeros(nindi,param);             %nt now

pfix=.5*ones(2,tmax);               %fixed environment, for testing purposes only
pfix(1,1:tmax)=0/10;
pfix(1,51:tmax)=10/10;
pfix(2,1:tmax)=8/10;
% pfix(1,151:200)=4.5/10;
% pfix(1,201:300)=3/10;
% pfix(1,251:300)=3.5/10;
% pfix(1,:)=5/10;
% pfix=[pA;pB];


% INITIAL POPULATION
ninitial(:,ichoice)=rand(nindi,1)>1/2;                %initialize random choice in 1st round
ninitial(:,iskill)=(-.5+rand(nindi,1))*pskill;        %skill is uniformly distributed with mean 0
ninitial(:,itally)=tallyn;                            %tally#

% ~~~~~~~~~~~~ change initial population here~~~~~~~~~~~~~~~
% use this procedure to define the initial population for
% pure strategies
ind=[...
    1000;          % INDividual learners
    0;          % MAJority talliers
    0;          % SUCcess talliers
    0;          % WEAlth imitators
    0;          % OMTs
    0;          % OILs
    0;          % REInforcement learners
    0];         % when in DOUbt, conform

cc=1;cind=cumsum(ind);pop=zeros(nindi,nstrat);
if cind(length(cind))~=nindi;'error in population length'
    break;end
for i=1:nstrat
    pop([cc:cind(i)],i)=1;
    cc=cc+ind(i);
end
ninitial(:,ipind:ipind+nstrat-1)=pop;
clear cc;clear pop;clear ind;clear cind;

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

% this is how the population tensor n will look like:
%
%             time
%              7\
%             /
%          4 / 1 / 0/11/ 7/ 7/-4/  3/  2/0.2/ 4/.4/.3/.2/.1/
%         3 / 1 / 0/11/ 7/ 7/-3/  3/  1/0.2/ 4/.4/.3/.2/.1/
%        2 / 0 / 1/11/ 7/ 7/-2/111/111/0.2/ 1/.4/.3/.2/.1/
%       1 / 1 / 0/10/ 7/ 7/-1/111/111/0.2/ 2/.4/.3/.2/.1/
%        /___/__/__/__/__/__/___/___/___/__/__/__/__/__/___\  parameters
%      1| c  s  f  t  m   b  t   r   s   u   p  p  p  p    /
%      2| h  u  i  a  e   i  r   e   k   s   r  r  r  r    
%      3| o  c  t  l  m   a  a   c   i   e   o  o  o  o      
%      4| i  c  t  l  o   s  c   e   l   d   b  b  b  b       
%      5| c  e  n  y  r      k   n   l   s   i  m  s  w      
%      6| e  s  e  #  y      i   t       t   n  a  u  e       
%      7|    s  s            n           r   d  j  c  a           
%       |       s            f           a   i  o  c  l                
%       V                                t                
%   individuals
%
%   t=2 ... (idem) ...
%
% dimension 1: individuals, 2: parameters, 3: strategy space, 4: time
% coding:   choice ->     1=A, 0=B
%           success ->    1=success, 0=failure
%           tally# ->     tally number
%           bias ->       bias towards choosing A
%           memory ->     same, but for reinforcement learners (=
%                         difference in propensities)
%           track info -> 1=indi, 2=majo, 3=succ, 4=wealth, ...
%           skill ->      adds or subtracts from prob. to succeed
%           probs ->      prob. of using 1] indi 2] majo 3]succ, etc.
%           strategy ->   which strategy is used eventually


% ERROR MESSAGES
% costs not implemented yet
if c>0||cm>0
    'error: costs not implemented in fitness calculation'
    bla^2
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

tic

%progression
if gen>10||iterations>10
    wb=waitbar(0,'progress');
end

% S T A R T   O F   T H E   S I M U L A T I O N
tic
countA=zeros(3,1);

for x=1:iterations

    for g=1:gen
        
        [n,pA,pB] = coevo25(...
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
        pfix,...    %a fixed environment for testing purposes
        param,...   %# of parameter values
        q,...       %oblivousness, discount for older memory
        lambda,...  %sensitivity factor for reinforcement learners, higher->steeper
        genetics,...%whether pure (0) or mixed strategies (1)
        kdoubt,...  %threshold value for ODCs
        aversion,...%whether positive/negative information is taken into account
        hyperb,...  %hyperbolic discounting: the "bonus" (>1) for most recent experience
        regime);    %regime, how the environment changes. 1->low variance, 2->high variance

        % I N H E R I T A N C E
        ninitial2=ninitial;                 %copy initial state
        if ntsize<=50000000
            nt(:,:,g)=n(:,:,tmax);              %the final state
            % monitor mean proportion of correct choices
            fig5(g)=coevo_meancorrect3(pA,pB,n,nindi,ichoice,tmax);
        else
            if mod(g,10)==0
                nt(:,:,g/10)=n(:,:,tmax);
                fig5(g/10)=coevo_meancorrect3(pA,pB,n,nindi,ichoice,tmax);
            end
        end
        ntn=n(:,:,tmax);
        
        % COPY INHERITABLE CHARACTERS + MUTATE
        ninitial=coevo_next_gen4(nindi,g,param,ntn,mutationrate,...
            nstrat,ichoice,ifit,itally,imemo,iskill,ipind,mutationincr);
        
        % include a check for EXTINCTION event
        if mod(g,50)==0
            [unik,unikchar]=coevo_find_unique(n,istrat,nstrat);
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
        [unik,unikchar]=coevo_find_unique(n,istrat,nstrat);
        meancorrstat=zeros(length(unik),iterations);
        majority_correct=zeros(1,iterations);
        wealth_info=zeros(length(unik),tmax*iterations);
        wealth_info=reshape(wealth_info,length(unik),tmax,iterations);
        wealth_age=zeros(tmax,iterations);
        correlstat=zeros(length(unik),iterations);
        correspstat=zeros(length(unik),iterations);
        sourcestat=zeros(length(unik),iterations);
        if monitorCrosscorr==1
            % the behavior of the strategies
            superfig=zeros(length(unik),tmax*iterations);
            % the environment
            superp=zeros(1,tmax*iterations);
%             % use this to calculate the delay:
%             [XCF,Lags,Bounds]=crosscorr(superp,superfig,100);
%             plot(Lags,XCF);
        end
    end
    
    % for more than 1 iteration:
%     % statistics on mean correct choices
    [meancorr,unik,unikchar,majority_was_right,corresp]=...
        coevo_meancorrect8(pA,pB,n,tmax,nindi,ichoice,...
        istrat,ibest,nstrat,genetics);
    meancorrstat(:,x)=meancorr;
    majority_correct(x)=majority_was_right;
    correspstat(:,x)=corresp;
    % how often has this strategy been the role model for wea?
    for s=1:length(unik)
        sourcestat(s,x)=mean(n(n(:,istrat,1)==unik(s),isource,tmax));
    end
    if monitorCrosscorr==1
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
    end
    
    % statistics on source and age of information from
    % wealth talliers
    if find(unik==4)~=0
        [fig_wealth_info,fig_wealth_age,unikchar2]=...
            coevo_wealth_info_age(unik,tmax,n,itrack,istrat,unikchar,irecent);
        wealth_info(:,:,x)=fig_wealth_info;
        wealth_age(:,x)=fig_wealth_age(1,:);
    end
    
    % for more than 1 iteration: use original initial state
    % as input for the next generation, except that choices
    % are inherited
    if iterations>1
        choicesvec=ninitial(:,ichoice);
        ninitial=ninitial2;
        ninitial(:,ichoice)=choicesvec;
        % distribute skill anew
        ninitial(:,iskill)=(-.5+rand(nindi,1))*pskill;
    end
    
%     % how fitness develops over time
%     % relative difference between wea and ind
%     if x==1
%         fitdiff=zeros(1,tmax);
%     end
%     meanfitind=zeros(1,tmax);
%     meanfitwea=zeros(1,tmax);
%     for tt=1:tmax
%         meanfitind(tt)=mean(n(n(:,istrat,tt)==1,ifit,tt));
%         meanfitwea(tt)=mean(n(n(:,istrat,tt)==4,ifit,tt));
%     end
%     fitdiff=fitdiff+(meanfitwea-meanfitind)./meanfitind;
%     if x==iterations
%         fitdiff=fitdiff/iterations;
%     end

%    ninitial(:,ichoice)=rand(nindi,1)>1/2;
%     if mean(n(101:1000,ichoice,tmax))>4/5
%         countA(1)=countA(1)+1;
%     elseif mean(n(101:1000,ichoice,tmax))<4/5
%         countA(2)=countA(2)+1;
%     else
%         countA(3)=countA(3)+1;
%     end
%     
%     if (mod(x,25)==1)&&(x<102)
%         mean(meancorrstat(2,1:x))
%         countA
%     end
end

% E N D   O F   S I M U L A T I O N

% in case memory would be exceeded
if ntsize>50000000
    gen=gen/10;
end

% R E S U L T S

% all unique strategies

if iterations==1
    
    if gen==1
        fig1=coevo_plot_1gen2(n,nt,nindi,tmax,nstrat,genetics,pA,pB,...
            istrat,ichoice,itrack,irecent,ifit);
        if genetics==0
            '    strat     mean correct'
            [unik meancorrstat]
        elseif genetics==1
            '    strat     mean correct'
            [unik repmat(meancorrstat,length(unik),1)]
        end
        [switchfreq,smat]=coevo_switchfreq(...
            pA,pB,n,tmax,nindi,ichoice,istrat,nstrat,genetics);
        switchfrequency=switchfreq
    end
    
    if gen>1
        fig1=coevo_plot_gens2(n,nt,nindi,tmax,nstrat,genetics,gen,...
            istrat,ichoice,fig5,ipind,ntsize);
    end

elseif iterations>1
    % statistics
    if genetics==0
        % frequency of the different strategies
        stratfreq=zeros(1,length(unik));
        for s=1:length(unik); stratfreq(s)=sum(n(:,istrat)==unik(s));end
        stratfreq=stratfreq/sum(stratfreq);
        % normalize sourcestat
        if length(find(unik==4))==1
            sourcestat=sourcestat/stratfreq(find(unik==4));
        end
        % create cell that contains summary statistics
        cellu=cell(2+length(unik),5);
        cellu(1,1)=cellstr('Strategy');
        cellu(1,2)=cellstr('mean corr');
        cellu(1,3)=cellstr('% maj. right');
%         cellu(1,4)=cellstr('correla°');
        cellu(1,4)=cellstr('corresp.');
        cellu(1,5)=cellstr('source');
        cellu(2:1+length(unik),1)=unikchar;
        cellu(2+length(unik),1)=cellstr('group mean');
        cellu(2:1+length(unik),2)=num2cell(mean(meancorrstat,2));
        cellu(2+length(unik),2)=num2cell(sum(stratfreq.*mean(meancorrstat,2)'));
        cellu(2+length(unik),3)=num2cell(mean(majority_correct)); 
%         cellu(2:1+length(unik),4)=num2cell(mean(correlstat,2));
%         cellu(2+length(unik),4)=num2cell(sum(stratfreq.*mean(correlstat,2)'));
        cellu(2:1+length(unik),4)=num2cell(mean(correspstat,2));
        cellu(2+length(unik),4)=num2cell(sum(stratfreq.*mean(correspstat,2)'));
        cellu(2:1+length(unik),5)=num2cell(mean(sourcestat,2));
        cellu(2+length(unik),5)=num2cell(sum(stratfreq.*mean(sourcestat,2)'));
        stats=cellu
        
    elseif genetics==1
        used_strategies=unikchar;
        probabilities=ninitial(1,[ipind:ipind+nstrat-1]);
        probabilities=probabilities(find(probabilities~=0));
        used_strategies(:,2)=cellstr(num2str(probabilities'))
        mean_correct=mean(meancorrstat(1,:),2)
        geomean_correct=geomean(meancorrstat(1,:)')'
        majority_chose_correctly=mean(majority_correct)
    end
        
end

toc
if gen>10||iterations>10
    close(wb) %close waitbar
end

