function [n,pA,pB] = coevo26(...
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
        lambda,...  %sensitivity factor for reinforcement learners; higher->steeper
        genetics,...%whether pure (0) or mixed strategies (1)
        kdoubt,...  %threshold value for ODCs
        aversion,...%whether positive/negative info is taken into account
        regime);    %regime, how the environment changes. 1->low variance, 2->high variance


% INDICES
ichoice=1;isucc=2;ifit=3;itally=4;imemo=5;ibias=6;
itrack=7;irecent=8;iskill=9;istrat=10;
ibest=11;isource=12;
ipind=13;ipmaj=14;ipsuc=15;ipwea=16;ipomt=17;ipoil=18;iprei=19;
ipidc=20;ipomn=21;

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

%TRACK whom wealth imitators copy and recentness of this information
n(:,itrack,:) =8*ones(nindi,tmax);  %8 codes for the non-wealth tallier
n(:,itrack,1) =9*ones(nindi,1);     %9 codes for initial guesses
n(:,irecent,1)=zeros(nindi,1);      %0 codes for non-wealth tallying


%create a RANDOM ENVIRONMENT according to specification
[pA,pB] = randomenvironment2(tmax,regime,incr,pincr,pA0,pB0);
% %use this fixed environment only for testing purposes:
% pA=pfix(1,:);
% pB=pfix(2,:);
% 'achtung pfix'

% pA=pA-.1;pB=pB-.1;

% vector containing the SKILL
skillvec=ninitial(:,iskill);

%random matrix used for determining success
randsucmatA=rand(nindi,tmax);
randsucmatB=rand(nindi,tmax);

% determine which STRATEGY is used in EACH PERIOD,
% according to the PROBABILITY distribution
n=coevo_determine_strat2(n,nindi,nstrat,istrat,ipind,tmax);

for t=1:tmax
    
    % S U C C E S S
    % update success
    n(:,isucc,t)=n(:,ichoice,t).*(randsucmatA(:,t)-skillvec<pA(t))+...
        mod(n(:,ichoice,t)+1,2).*(randsucmatB(:,t)-skillvec<pB(t));
    % update BIAS towards A or B for strategies relying on individual learning
    if aversion==0      %wins and losses count equally
        n(:,ibias,t+1)=q*n(:,ibias,t)+...   %discount older bias
            ((n(:,ichoice,t)==n(:,isucc,t))*2-1);
    else
        n(:,ibias,t+1)=q*n(:,ibias,t)+...   %discount older bias
            ((n(:,ichoice,t)==1)&(n(:,isucc,t)==1))-... %choice A + suc.
            ((n(:,ichoice,t)==0)&(n(:,isucc,t)==1))-... %choice B + suc.
            aversion*((n(:,ichoice,t)==1)&(n(:,isucc,t)==0))+... %choice A + fail
            aversion*((n(:,ichoice,t)==0)&(n(:,isucc,t)==0)); %choice B + fail
    end
    
    % bias+1 if A succeeded or B failed, -1 vice versa
    
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
    
    % INDIDIVUAL LEARNING
    scvec=scvec+(n(:,istrat,t)==1).*...     %if indiv. learning this round
        n(:,ibias,t+1);                     %score corresponds to bias
    
    % MAJORITY TALLYING
    
    % matrix containing random individuals who are sampled
    randtmat=ceil(nindi*rand(nindi,max(n(:,itally,1))));
    
    for i=1:nindi
        if n(i,istrat,t)==2                                 %if individual is majority tallier
            for ii=1:n(i,itally,t)                          %screen a # of indis according to tally #
                scvec(i)=scvec(i)+...                       
                    (2*(n(randtmat(i,ii),ichoice,t)==1)-1); %favor corresponding choice
            end
        end
    end
        
    % SUCCESS TALLYING
    if aversion==0      %if wins and losses count equally
        for ii=1:n(i,itally,t)
            scvec=scvec+(n(:,istrat,t)==3).*...
                (2*(n(randtmat(:,ii),ichoice,t)==n(randtmat(:,ii),isucc,t))-1);
        end
    else
        for ii=1:n(i,itally,t)
            scvec=scvec+(n(:,istrat,t)==3).*(...
                ((n(randtmat(:,ii),ichoice,t)==1)&(n(randtmat(:,ii),isucc,t)==1))-... %choice A + suc
                ((n(randtmat(:,ii),ichoice,t)==0)&(n(randtmat(:,ii),isucc,t)==1))-... %choice B + suc
                aversion*((n(randtmat(:,ii),ichoice,t)==1)&(n(randtmat(:,ii),isucc,t)==0))+... %choice A + fail
                aversion*((n(randtmat(:,ii),ichoice,t)==0)&(n(randtmat(:,ii),isucc,t)==0))); %choice B + fail
        end
    end
    
%     _____________________________________________________
    
    % IMITATE THE WEALTHIEST
    talmax=max(n(:,itally,1));
    %snapshot of present population
    nnow=n(:,:,t);
    %add a random decimal to fitness so that no 2 fitnesses are the same
    %if you don't do this, draws in fitness will always be sorted the
    %same way, which might lead to biases in whom ITW imitates
    nnow(:,ifit)=nnow(:,ifit)+b*rand(nindi,1);
    % before sorting, remember the initial position
    nnow=[nnow [1:nindi]'];
    %sort individuals according to fitness (ascending)
    nnow=sortrows(nnow,ifit);
    
    %create the randomly queued population  
    randtmat=1+nindi*rand(nindi,talmax);
    randtmat2=randtmat';
    
    % find coordinate of the highest ranking individual that is tallied
    c1=find((randtmat==repmat(max(randtmat')',1,talmax))');
    % since the tallied individuals are sorted according to their
    % fitness, choosing the first individual is equal to choosing
    % the fittest individual.
    % choose the same as this individual
    scvec=scvec+(n(:,istrat)==4).*(2*nnow(floor(randtmat2(c1)),ichoice)-1);
    
    % for a tally number of 1, the routine does not work, since max(Y')=/=Y'
    % for size(Y)=(X,1), thus we take this:
    if talmax==1
        scvec(n(:,istrat)==4)=(2*nnow(floor(randtmat2(n(:,istrat)==4)),ichoice)-1);
    end
    
    % track source and age of INFORMATION of wealth imitators
    n(:,itrack,t+1)=(n(:,istrat,t)==4).*...
        nnow(floor(randtmat2(c1)),istrat);
    
    % AGE
    n(:,irecent,t+1)=n(:,irecent,t)+(n(:,istrat,t)==4);
    n(find(n(:,itrack,t+1)~=4),irecent,t+1)=1;
    
    % SOURCE
    % if information comes from wealth imitator, check whom SHE has
    % imitated
    n(:,itrack,t+1)=(n(:,itrack,t+1)==4).*...
        nnow(floor(randtmat2(c1)),itrack)+...
        (n(:,itrack,t+1)~=4).*n(:,itrack,t+1);
    % keep track of every time an individual is imitated by a wealth
    % imitator
    [s1 s2]=size(nnow);
    for jj=1:nindi
        n(nnow(floor(randtmat2(c1(jj))),s2),isource,t)=...
            n(nnow(floor(randtmat2(c1(jj))),s2),isource,t)+...
            (n(jj,istrat,1)==4);
    end
    %_______________________________________________________
    
    % HIERARCHICAL 1: opportunistic majority tallier (OMT)
    % if successful, use ind. learning, else majority tallying
    scvec=scvec+(n(:,istrat,t)==5).*...     %if OMT this round
        (n(:,isucc,t)==1).*...              %if successful this round
        n(:,ibias,t+1);                     %score corresponds to bias
    randtmat=ceil(nindi*rand(nindi,max(n(:,itally,1))));
    for i=1:nindi
        if n(i,istrat,t)==5
            for ii=1:n(i,itally,t)                          %screen a # of indis according to tally #
                scvec(i)=scvec(i)+...                       
                    (n(i,isucc,t)==0).*...                  %if not successful this round
                    (2*(n(randtmat(i,ii),ichoice,t)==1)-1); %favor corresponding choice
            end
        end
    end
    
    % HIERARCHICAL 2: opportunistic individual learner (OIL)
    % if successful, use majority tallying, else individual learning
    scvec=scvec+(n(:,istrat,t)==6).*...     %if OIL this round
        (n(:,isucc,t)==0).*...              %if not successful this round
        n(:,ibias,t+1);                     %score corresponds to bias
    randtmat=ceil(nindi*rand(nindi,max(n(:,itally,1))));
    for i=1:nindi
        if n(i,istrat,t)==6
            for ii=1:n(i,itally,t)                          %screen a # of indis according to tally #
                scvec(i)=scvec(i)+...                       
                    (n(i,isucc,t)==1).*...                  %if successful this round
                    (2*(n(randtmat(i,ii),ichoice,t)==1)-1); %favor corresponding choice
            end
        end
    end
    
    % HIERARCHICAL 3: when in DOUBT, conform (use majority tallying)
    % a strategy is in doubt if its bias does not favor any of the
    % options in particual
    % threshold: kdoubt
    scvec=scvec+(n(:,istrat,t)==8).*...     %if IDC this round
        (abs(n(:,ibias,t))>kdoubt).*...     %if one option is strongly favored
        n(:,ibias,t+1);                     %score corresponds to bias
    randtmat=ceil(nindi*rand(nindi,max(n(:,itally,1))));
    for i=1:nindi
        if n(i,istrat,t)==8
            for ii=1:n(i,itally,t)                          %screen a # of indis according to tally #
                scvec(i)=scvec(i)+...                       
                    (abs(n(i,ibias,t))<=kdoubt).*...             %if undecided this round
                    (2*(n(randtmat(i,ii),ichoice,t)==1)-1); %favor corresponding choice
            end
        end
    end
        
        
    % PROBABILISTIC REINFORCEMENT LEARNING (REI)
    randthresh=rand(nindi,1);   %threshold for choice of A or B
    probrei=zeros(nindi,1);     %vector of the probabilities to choose A
    % the difference between the propensity to choose A and the propensity
    % to choose b:
    
    if aversion==0      %wins and losses count equally
        n(:,imemo,t+1)=q*n(:,imemo,t)+...   %discount older bias
            ((n(:,ichoice,t)==n(:,isucc,t))*2-1);
    elseif aversion<0   %only negative info taken into account
        n(:,imemo,t+1)=q*n(:,imemo,t)+...
            abs(aversion)*((n(:,ichoice,t)==0)&(n(:,isucc,t)==0))-...
            abs(aversion)*((n(:,ichoice,t)==1)&(n(:,isucc,t)==0))-...
            ((n(:,ichoice,t)==0)&(n(:,isucc,t)==1))+...
            ((n(:,ichoice,t)==1)&(n(:,isucc,t)==1));
    elseif aversion>0   %only positive info taken into account
        n(:,imemo,t+1)=q*n(:,imemo,t)-...
            aversion*((n(:,ichoice,t)==0)&(n(:,isucc,t)==1))+...
            aversion*((n(:,ichoice,t)==1)&(n(:,isucc,t)==1))+...
            ((n(:,ichoice,t)==0)&(n(:,isucc,t)==0))-...
            ((n(:,ichoice,t)==1)&(n(:,isucc,t)==0));
    end
    
    probrei=1./(1+exp(lambda*(-n(:,imemo,t+1))));
    scvec=scvec+(n(:,istrat,t)==7).*(2*(probrei>randthresh)-1);
    
    
    % OMNISCIENT STRATEGY (OMN)
    if t<tmax %(necessary; if not, pA(tmax+1) is called
        if pA(t+1)==pB(t+1) %if A and B equally good, choose randomly
            scvec=scvec+(n(:,istrat,t)==9).*(2*(rand(nindi,1)>.5*ones(nindi,1))-1);
        else %if A better than B, choose A and v.v.
            scvec=scvec+(n(:,istrat,t)==9).*repmat(2*(pA(t+1)>pB(t+1))-1,nindi,1);
        end
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