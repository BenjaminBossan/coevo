% this function serves to determine which strategy is chosen each
% period according to a probability function that is given. The
% choice of strategy is then plugged into the corresponding slot

function n=coevo_determine_strat2(n,nindi,nstrat,istrat,ipind,tmax);

randstratvec=rand(nindi,tmax);

for t=1:tmax
    for s=1:nstrat
        n(:,istrat,t)=n(:,istrat,t)+s*ones(nindi,1).*(n(:,istrat,t)==0).*...
            (sum(n(:,ipind:ipind+s-1,t),2)>=randstratvec(:,t));
    end
end