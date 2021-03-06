function [n,pA,pB] = coevo51(...
        tmax,...    %tmax time of simulation
        nindi,...   %nindi # of individuals
        nstrat,...  %# of strategy types
        tallyn,...  %sample size itw
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
if length(ninitial(:, istrat)) == 1
    n = reshape(repmat(ninitial, 1, tmax), nindi, param, tmax);
else
    n(:,:,1)=ninitial;              %the initial state of the population

    n(:,iskill,:)=repmat(ninitial(:,iskill),1,tmax);    % skill stays same
    for s=1:nstrat
        n(:,ipind+s-1,:)=repmat(ninitial(:,ipind+s-1),1,tmax);
    end

    % strategy stays the same
    n(:,istrat,:)=repmat(ninitial(:,istrat),1,tmax);    % skill stays same

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

    % % determine which STRATEGY is used in EACH PERIOD,
    % % according to the PROBABILITY distribution
    % if sum(mod(unique(n(:,ipind:ipind+nstrat-1,1)),1))==0
    %     %if there are only pure strategies
    %     for jj=1:nindi
    %         n(jj,istrat,1:tmax)=find(n(jj,ipind:ipind+nstrat-1,1)==1);
    %     end
    % else
    %     %if there are mixed strategies
    %     n=coevo_determine_strat2(n,nindi,nstrat,istrat,ipind,tmax);
    % end

    for t=1:tmax

        % random vector with size of population
        randIndi=rand(nindi,1);

        % S U C C E S S
        % update success
        n(:,isucc,t)=n(:,ichoice,t).*(randsucmatA(:,t)-skillvec<pA(t))+...
            mod(n(:,ichoice,t)+1,2).*(randsucmatB(:,t)-skillvec<pB(t));

        % CHOICE FREQUENCIES
        % frequency of A choices in last period:
        x=sum(n(:,ichoice,t))/nindi;
        % freq of successful A choices
        xAs=sum(n(:,ichoice,t)&n(:,isucc,t))/nindi;
        % freq of failed A choices
        xAf=sum(n(:,ichoice,t)&(n(:,isucc,t)==0))/nindi;
        % freq of successful B choices
        xBs=sum((n(:,ichoice,t)==0)&n(:,isucc,t))/nindi;
        % freq of failed B choices
        xBf=sum((n(:,ichoice,t)==0)&(n(:,isucc,t)==0))/nindi;


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

        % 1 INDIDIVUAL LEARNING
        scvec=scvec+(n(:,istrat,t)==1).*...     %if indiv. learning this round
            n(:,ibias,t+1);                     %score corresponds to bias


        % 2 CONFORMISM

        scvec=scvec+(n(:,istrat,t)==2).*(2*(randIndi<((3-2*x)*x^2))-1);


        % 3 OPPORTUNISTIC INDIVIDUAL LEARNERS

        % if successful last period, use conformism
        index1=((n(:,istrat,t)==3)&(n(:,isucc,t)==1));
        scvec=scvec+index1.*(2*(randIndi<((3-2*x)*x^2))-1);

        % if unsuccessful last period, use individual learning
        scvec=scvec+((n(:,istrat,t)==3)&(n(:,isucc,t)==0)).*n(:,ibias,t+1);


        % 4 OPPORTUNISTIC CONFORMISTS

        % if unsuccessful last period, use conformism
        index1=((n(:,istrat,t)==4)&(n(:,isucc,t)==0));
        scvec=scvec+index1.*(2*(randIndi<((3-2*x)*x^2))-1);

        % if successful last period, use individual learning
        scvec=scvec+((n(:,istrat,t)==4)&(n(:,isucc,t)==1)).*n(:,ibias,t+1);


        % 5 IN DOUBT, CONFORM

        % if there is doubt, use conformism
        index1=((n(:,istrat,t)==5)&(abs(n(:,ibias,t))<kdoubt));
        scvec=scvec+index1.*(2*(randIndi<((3-2*x)*x^2))-1);

        % if certain, use individual learning
        scvec=scvec+(n(:,istrat,t)==5).*(abs(n(:,ibias,t))>=kdoubt).*n(:,ibias,t+1);


        % 6 IMITATE THE WEALTHIEST

        if length(find(n(:,istrat,t)==6))>0 %if there are ITW users at all
            talmax=7;
            %snapshot of present population
            nnow=squeeze(n(:,:,t));
            %add a random decimal to fitness so that no 2 fitnesses are the same
            %if you don't do this, draws in fitness will always be sorted the
            %same way, which might lead to biases in whom ITW imitates
            nnow(:,ifit)=nnow(:,ifit)+b*rand(nindi,1);
            % before sorting, remember the initial position
            nnow=[nnow [1:nindi]'];
            %sort individuals according to fitness (ascending)
            nnow=sortrows(nnow,ifit);

            %create the randomly queued population
            if compare_self==0
                randtmat=1+nindi*rand(nindi,talmax);
            else
                randtmat=1+nindi*rand(nindi,talmax-1);
                % if ITW also compares with its own payoff, it is as if each ITW
                % user had herself within her sample.
                randtmat=[[1:nindi]' randtmat];
            end
            randtmat2=randtmat';

            % find coordinate of the highest ranking individual that is tallied
            c1=find((randtmat==repmat(max(randtmat')',1,talmax))');
            % since the tallied individuals are sorted according to their
            % fitness, choosing the first individual is equal to choosing
            % the fittest individual.
            % choose the same as this individual
            scvec=scvec+(n(:,istrat)==6).*(2*nnow(floor(randtmat2(c1)),ichoice)-1);

            % for a tally number of 1, the routine does not work, since max(Y')=/=Y'
            % for size(Y)=(X,1), thus we take this:
            if talmax==1
                scvec(n(:,istrat)==6)=(2*nnow(floor(randtmat2(n(:,istrat)==6)),ichoice)-1);
            end

            % track source and age of INFORMATION of wealth imitators
            n(:,itrack,t+1)=(n(:,istrat,t)==6).*...
                nnow(floor(randtmat2(c1)),istrat);

            % AGE
            n(:,irecent,t+1)=n(:,irecent,t)+(n(:,istrat,t)==6);
            n(find(n(:,itrack,t+1)~=6),irecent,t+1)=1;

            % SOURCE
            % if information comes from wealth imitator, check whom SHE has
            % imitated
            n(:,itrack,t+1)=(n(:,itrack,t+1)==6).*...
                nnow(floor(randtmat2(c1)),itrack)+...
                (n(:,itrack,t+1)~=6).*n(:,itrack,t+1);
            % keep track of every time an individual is imitated by a wealth
            % imitator
            [s1 s2]=size(nnow);
            for jj=1:nindi
                n(nnow(floor(randtmat2(c1(jj))),s2),isource,t)=...
                    n(nnow(floor(randtmat2(c1(jj))),s2),isource,t)+...
                    (n(jj,istrat,1)==6);
            end
        end
        %_______________________________________________________


        % 7 PBSL [4/-1]

        scvec=scvec+(n(:,istrat,t)==7).*...
            (-1+2*(randIndi<xAs^7+21*xAf^5*xAs*(xAs+xBf)+7*xAs^6*(xBf+xBs)+21*xAs^5*(xBf+xBs)^2+35*xAs^4*(xBf+xBs)^3+xBf^6*(xBf+7*xBs)+21*xAs^2*xBf^3*(xBf^2+5*xBf*xBs+10*xBs^2)+(7*xAs*xBf^4*(2*xBf^2+12*xBf*xBs+15*xBs^2))/2+35*xAs^3*xBf*(xBf^3+4*xBf^2*xBs+6*xBf*xBs^2+4*xBs^3)+(35*xAf^4*xAs*(2*xAs^2+6*xBf^2+3*xAs*(2*xBf+xBs)))/2+35*xAf^3*(xAs^4+4*xAs*xBf^3+xBf^4+4*xAs^3*(xBf+xBs)+6*xAs^2*xBf*(xBf+2*xBs))+21*xAf^2*(xAs^5+xBf^5+5*xAs^4*(xBf+xBs)+10*xAs^3*(xBf+xBs)^2+10*xAs^2*xBf^2*(xBf+3*xBs)+5*xAs*xBf^3*(xBf+4*xBs))+7*xAf*(xAs^6+6*xAs^5*(xBf+xBs)+15*xAs^4*(xBf+xBs)^2+xBf^5*(xBf+3*xBs)+6*xAs*xBf^4*(xBf+5*xBs)+20*xAs^3*xBf*(xBf^2+3*xBf*xBs+3*xBs^2)+15*xAs^2*xBf^2*(xBf^2+4*xBf*xBs+6*xBs^2))));


        % 8 PBSL [1/0]
        scvec=scvec+(n(:,istrat,t)==8).*...
            (-1+2*(randIndi<(xAf^3+2*xAs^3+xBf^3+3*xAf^2*(2*xAs+xBf)+6*xAs^2*(xBf+xBs)+6*xAs*xBf*(xBf+xBs)+3*xAf*(2*xAs^2+xBf^2+2*xAs*(2*xBf+xBs)))/2));


        % 9 PBSL McElreath
        scvec=scvec+(n(:,istrat,t)==9).*...
            (-1+2*(randIndi<xAf^3+3*xAf^2*(xAs+xBf)+3*xAf*xAs*(xAs+2*xBf)+xAs*(xAs^2+3*xAs*(xBf+xBs)+3*xBf*(xBf+2*xBs))));


        % 10 PBSL PAYOFF-CONFORMISM TRADE-OFF
        scvec=scvec+(n(:,istrat,t)==10).*...
            (-1+2*(randIndi<xAf^6+6*xAf^5*(xAs+xBf)+15*xAf^4*(xAs+xBf)^2+10*xAf^3*(2*xAs^3+6*xAs^2*xBf+xBf^3+6*xAs*xBf*(xBf+xBs))+15*xAf^2*xAs*(xAs^3+4*xAs^2*(xBf+xBs)+6*xAs*xBf*(xBf+2*xBs)+2*xBf^2*(2*xBf+3*xBs))+6*xAf*xAs*(xAs^4+5*xAs^3*(xBf+xBs)+5*xBf^3*(xBf+2*xBs)+5*xAs^2*(2*xBf^2+4*xBf*xBs+xBs^2)+5*xAs*xBf*(2*xBf^2+6*xBf*xBs+3*xBs^2))+xAs*(xAs^5+6*xAs^4*(xBf+xBs)+15*xAs^3*(xBf+xBs)^2+6*xBf^3*(xBf^2+5*xBf*xBs+10*xBs^2)+10*xAs^2*(2*xBf^3+6*xBf^2*xBs+6*xBf*xBs^2+xBs^3)+15*xAs*xBf*(xBf^3+4*xBf^2*xBs+6*xBf*xBs^2+2*xBs^3))));



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
end