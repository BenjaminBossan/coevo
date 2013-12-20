% determines the next generation according to the fitness of the
% individuals. first fitness is normalized, then individuals are
% randomly drawn from the population and a random function deter-
% mines whether this individual reproduces. This function works
% so that the individual with the highest fitness has a chance of
% 'probrepbest' to reproduce, the worst individual of zero.

function ninitial=coevo_next_gen4(nindi,g,param,nt,mutationrate,...
    nstrat,ichoice,ifit,itally,imemo,iskill,ipind,mutationincr);

ninitial=zeros(nindi,param);

% extract relevant information
ntf=nt(:,ifit);    %fitness
ntc=nt(:,ichoice); %choice
ntt=nt(:,itally);  %tally #
ntk=nt(:,iskill);  %skill
nts=nt(:,ipind:ipind+nstrat-1);  %strategy


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
    
%     really????
%     ninitial(i,iskill)=ntk(offspringind(randvec(i)));
%     %next gen copies skill level

    ninitial(i,ipind:ipind+nstrat-1)=...
        nts(offspringind(randvec(i)),:);
    %next gen copies strategy vector
end



% MUTATION

if mutationrate>0
    %whether mutation takes place
    randmat=rand(nindi,nstrat);
    %whether mutation leads to increase or decrease
    riseorfall=sign(rand(nindi,nstrat)-.5);
    for s=1:nstrat
        ninitial(:,ipind+s-1)=ninitial(:,ipind+s-1)+...              %old probs
            ninitial(:,ipind+s-1).*(randmat(:,s)<mutationrate).*...  %whether mutation occurs
            riseorfall(:,s)*mutationincr;                            %increase or decrease
    end
    % normalize:
    norma=sum(ninitial(:,ipind:ipind+nstrat-1),2);
    norma=repmat(norma,1,nstrat);
    ninitial(:,ipind:ipind+nstrat-1)=...
        ninitial(:,ipind:ipind+nstrat-1)./norma;
end



    
