% determines the next generation according to the fitness of the
% individuals. first fitness is normalized, then individuals are
% randomly drawn from the population and a random function deter-
% mines whether this individual reproduces. This function works
% so that the individual with the highest fitness has a chance of
% 'probrepbest' to reproduce, the worst individual of zero.

function ninitial=coevo_next_gen4(nindi,g,param,nt,mutationrate,...
    nstrat,ichoice,ifit,itally,istrat,iskill,ipind,mutationincr);

ninitial=zeros(nindi,param);

% extract relevant information
ntf=nt(:,ifit);    %fitness
ntc=nt(:,ichoice); %choice
ntt=nt(:,itally);  %tally #
ntk=nt(:,iskill);  %skill
nts=nt(:,istrat);  %strategy


% REPLICATION

% vector which will contain indices of offspring
offspringind=zeros(sum(nt(:,ifit)),1);
noffspring=length(offspringind);

% # of offspring proportional to fitness
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

% take as many individuals as population should contain
randvec=randvec2(1:nindi);

for i=1:nindi
    ninitial(i,ichoice)=ntc(offspringind(randvec(i)));
    %next gen copies last choice as initial choice

    ninitial(i,itally)=ntt(offspringind(randvec(i)));
    %next gen copies tally #

    ninitial(i,istrat)=nts(offspringind(randvec(i)));
end



    
