% determines the next generation according to the fitness of the
% individuals. first fitness is normalized, then individuals are
% randomly drawn from the population and a random function deter-
% mines whether this individual reproduces. This function works
% so that the individual with the highest fitness has a chance of
% 'probrepbest' to reproduce, the worst individual of zero.

% this specific instance inherits and mutates the variable weights
% of gains and losses of PBSLs. Mutations do not affect strategy type

function ninitial=coevo_next_gen_pbsl_01(nindi,g,param,nt,mutationrate,...
    nstrat,ichoice,ifit,itally,imemo,iskill,ipind,mutationincr,weightscope,...
    igain,iloss);

% note: the variable "nt" is the variable "ntn" in the runsim code

ninitial=zeros(nindi,param);

% extract relevant information
ntf=nt(:,ifit);    %fitness
ntc=nt(:,ichoice); %choice
ntt=nt(:,itally);  %tally #
ntk=nt(:,iskill);  %skill
nts=nt(:,ipind:ipind+nstrat-1);  %strategy
ntgain=nt(:,igain);    %weight of gains
ntloss=nt(:,iloss);    %weight of losses


% REPLICATION

% vector which will contain indices of offspring
offspringind=zeros(sum(nt(:,ifit)),1);
noffspring=length(offspringind);

counter=0;
counter2=0;
for i=1:nindi
    counter=counter+1;
    for j=1:ntf(i)
        counter2=counter2+1;
        offspringind(counter2)=counter;
    end
end

% queue offspring in random fashion
randvec2=randperm(noffspring);

% take the first x individuals with x=pop size
randvec=randvec2(1:nindi);

for i=1:nindi
    ninitial(i,ichoice)=ntc(offspringind(randvec(i)));
    %next gen copies last choice as initial choice

    ninitial(i,itally)=ntt(offspringind(randvec(i)));
    %next gen copies tally #

    ninitial(i,ipind:ipind+nstrat-1)=...
        nts(offspringind(randvec(i)),:);
    %next gen copies strategy vector
    
    ninitial(i,igain)=ntgain(offspringind(randvec(i)));
    ninitial(i,iloss)=ntloss(offspringind(randvec(i)));
    %next gen copies weights of gains and losses
end



% MUTATION

if mutationrate>0
    %whether mutation takes place
    randmat=rand(nindi,1)<mutationrate;
    randmat2=rand(nindi,1)<mutationrate;
    %scope of the mutation
    mutgain=floor((2*weightscope+1)*rand(nindi,1))-weightscope;
    mutloss=floor((2*weightscope+1)*rand(nindi,1))-weightscope;
    %apply mutations
    ninitial(randmat,igain)=mutgain(randmat);
    ninitial(randmat2,iloss)=mutloss(randmat2);
end

%divide weights of gains and losses by greatest common divisor
gvec=ninitial(:,igain);
lvec=ninitial(:,iloss);
gcdvec=gcd(gvec,lvec);
gcdvec(gcdvec==0)=1;    % gcd of (0,0) is 0, so cannot divide by it
ninitial(:,[igain iloss])=[gvec./gcdvec lvec./gcdvec];



    