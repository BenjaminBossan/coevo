% the probability of an environment becoming more likely to lead to success
% in the next round is dependent on the probability that it lead to success
% in the current round: pX+=f(pX). Specifically, we define:
% pX+=min(0,max(1,0.5-r*(pX-0.5)))
% the value of r is in [0,infinity[
% regime 1 corresponds to r = 1
% regime 2 corresponds to r = 0
% regime 3 corresponds to r = 2
% in contrast, if regime<0, we implement a regime that only
% knows pX={0,1}~=pY, with switch probability=1/abs(regime)
% this is more close to the older models

function [pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);

incr2=incr/2;

% generate random vector for environmental change
% determines whether there is a switch at all
penvvecA=rand(1,tmax);
penvvecB=rand(1,tmax);
penvvecA=penvvecA<=pincr;
penvvecB=penvvecB<=pincr;
% determines whether A or B increases
envvecA=rand(1,tmax);
envvecB=rand(1,tmax);

% success rate of hunting grounds
% initialization
pA=zeros(1,tmax);
pB=pA;
pA(1)=pA0;
pB(1)=pB0;

% reversion parameter
r=regime;

% min and max p that can be attained
minpA=pA0-((pA0-mod(pA0,incr))/incr)*incr;
maxpA=pA0+((pA0-mod(pA0,incr))/incr)*incr;
minpB=pB0-((pB0-mod(pB0,incr))/incr)*incr;
maxpB=pB0+((pB0-mod(pB0,incr))/incr)*incr;


for t=2:tmax
    if penvvecA(t)==1
        if envvecA(t-1)<min(1,(max(0,.5-r*(pA(t-1)-.5))))
            pA(t)=min(pA(t-1)+incr,maxpA);
        else
            pA(t)=max(pA(t-1)-incr,minpA);
        end
    else
        pA(t)=pA(t-1);
    end
end
for t=2:tmax
    if penvvecB(t)==1
        if envvecB(t-1)<min(1,(max(0,.5-r*(pB(t-1)-.5))))
            pB(t)=min(pB(t-1)+incr,maxpB);
        else
            pB(t)=max(pB(t-1)-incr,minpB);
        end
    else
        pB(t)=pB(t-1);
    end
end
% the following steps are necessary because matlab seems,
% sometimes, to add miniscule amounts to pA or pB (1e-17)
% so that two seemingly equal values are not equal anymore
% then it may happen that pA=0.6000>0.6000=pB. I don't know
% whence this comes from, honestly...
% pA=round(pA/incr2)*incr2;
% pB=round(pB/incr2)*incr2;

pA=10000*pA;pA=round(pA);pA=pA/10000;
pB=10000*pB;pB=round(pB);pB=pB/10000;

    