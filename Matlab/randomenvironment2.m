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

function [pA,pB] = randomenvironment2(tmax,regime,incr,pincr,pA0,pB0);

incr2=incr/2;

% generate random vector for environmental change
% clear envvec
% clear envvecA
% clear envvecB
envvec=rand(1,tmax);
envvec=envvec<=pincr;

% success rate of hunting grounds
pA=zeros(1,tmax);
pB=pA;
pA(1)=pA0;
pB(1)=pB0;
envvecA=rand(1,tmax);
envvecB=rand(1,tmax);

if regime==1    % patch more unlikely to increase in quality of alreay good
    for t=1:tmax-1
        if envvec(t)==1
            if envvecA(t)>pA(t)
                pA(t+1)=min(pA(t)+incr,1-incr);
            else
                pA(t+1)=max(pA(t)-incr,incr);
            end
        else
            pA(t+1)=pA(t);
        end
    end
    for t=1:tmax-1
        if envvec(t)==1
            if envvecB(t)>pB(t)
                pB(t+1)=min(pB(t)+incr,1-incr);
            else
                pB(t+1)=max(pB(t)-incr,incr);
            end
        else
            pB(t+1)=pB(t);
        end
    end
    pA=round(pA/incr2)*incr2;
    pB=round(pB/incr2)*incr2;
elseif regime==0    % constant probability of increase in quality
    for t=1:tmax-1
        if envvec(t)==1
            if envvecA(t)>.5
                pA(t+1)=min(pA(t)+incr,1-incr);
            else
                pA(t+1)=max(pA(t)-incr,incr);
            end
        else
            pA(t+1)=pA(t);
        end
    end
    for t=1:tmax-1
        if envvec(t)==1
            if envvecB(t)>.5
                pB(t+1)=min(pB(t)+incr,1-incr);
            else
                pB(t+1)=max(pB(t)-incr,incr);
            end
        else
            pB(t+1)=pB(t);
        end
    end
    pA=round(pA/incr2)*incr2;
    pB=round(pB/incr2)*incr2;
elseif regime==2    % as regime 1, but even more severe regression to the mean
    for t=1:tmax-1
        if envvec(t)==1
            if envvecA(t)<min(1,max(0,0.5-2*(pA(t)-0.5)))
                pA(t+1)=min(pA(t)+incr,1-incr);
            else
                pA(t+1)=max(pA(t)-incr,incr);
            end
        else
            pA(t+1)=pA(t);
        end
    end
    for t=1:tmax-1
        if envvec(t)==1
            if envvecB(t)<min(1,max(0,0.5-2*(pB(t)-0.5)))
                pB(t+1)=min(pB(t)+incr,1-incr);
            else
                pB(t+1)=max(pB(t)-incr,incr);
            end
        else
            pB(t+1)=pB(t);
        end
    end
    pA=round(pA/incr2)*incr2;
    pB=round(pB/incr2)*incr2;
elseif regime<0
    switchvec=2*(rand(1,tmax)>(1/abs(regime)))-1;
    if (pA0==1/2)&(pB0==1/2)
        pA(1)=rand>1/2;
    else
        pA(1)=2*pA0-1;
    end
    for t=2:tmax
        pA(t)=pA(t-1)*switchvec(t);
    end
    pA=max(0,pA);
    pB=1-pA;
end