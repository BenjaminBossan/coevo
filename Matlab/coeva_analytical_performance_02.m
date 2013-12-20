% this programm just computes the behavior of the strategies
% without any kind of evolution.
% The goal is to find out the PERFORMANCE of the different strategies

clear all

% PARAMETERS
tmax=1000000;                    % periods per generation
nStrats=11;                 % total possible number of strategies
regime=1;                   % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;                 % increment at which the environment becomes better or worse
pincr=1;                    % probability that environmental quality changes at all after each period
dpA=0.25;dpB=0.25;                % shift of pA and pB
pA0=.5;pB0=.5;
indLearn=0;                 % if == 0, ind. learning is win-stay lose-shift, else its constant

% STRATEGIES
%  1 individual learners
%  2 conformists tn 3 [1/1]
%  3 contrarian tn 3 [-1/-1]
%  4 PBSL [1/-1]    equal weights
%  5 PBSL [3/-1]    gains-biased
%  6 PBSL [1/-3]    loss-averse
%  7 PBSL [1/0]     only gains
%  8 PBSL [0/-1]    only losses
%  9 PBSL [3/1]     pbsl-conformist
% 10 PBSL [-1/-3]   pbsl-contrarian, hates wins, hates losses even more
% 11 PBSL a la McElreath

% INITIALIZATION
xInit=zeros(nStrats,1);

xInit(1)=.2;
xInit(2)=.8;
% xInit(3)=.9;
% xInit(4)=1;
% xInit(5)=1;
% xInit(6)=.6;
% xInit(7)=1/7;
% xInit(8)=1;
% xInit(9)=.1;
% xInit(10)=1;
% xInit(11)=1;


% xInit=ones(nStrats,1)/nStrats;

c=zeros(nStrats,tmax);      % percent correct choices matrix
x=xInit;                    % frequency of the strategies

    
% ROUTINE
    
% ENVIRONMENT
% routine to determine pA and pB
[pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);
% possible shift
pA=min(1,max(0,pA+dpA));pB=min(1,max(pB+dpB,0));
% we define C is the better and F as the worse choice
pC=pA;
pF=pB;
% remember last pA and pB
pA0=pA(tmax);pB0=pB(tmax);

% in very first generation
c(:,1)=1/2;

tic

% ROUTINE
for t=2:tmax

    cmean=sum(x(:,1).*c(:,t-1));

    % individual learners
    if indLearn==0
        cind=(1-pF(t))/(2-pC(t)-pF(t));
    else
        cind=((1-2*(1-indLearn))*sign(pC(t)-pF(t))+1)/2;
    end
    c(1,t)=cind;

    % 02 conformists
    cconf3=(3-2*cmean)*cmean^2;
    c(2,t)=cconf3;
    
    % 03 contrarian
    c(3,t)=(-1+cmean).^2.*(1+2.*cmean);
    
    %  4 PBSL [1/-1]    equal weights
    c(4,t)=-((1-pF(t)+cmean.*(-1+pC(t)+pF(t))).^2.*(-1-2.*pF(t)+2.*cmean.*(-1+pC(t)+pF(t))));   
    
    %  5 PBSL [3/-1]    gains-biased
    c(5,t)=1+cmean.^2.*(-1+pC(t)).^2.*(-3+cmean.*(2+pC(t)))-3.*(-1+cmean).*(-1+cmean.*(cmean+2.*(1+cmean.*(-2+pC(t))).*pC(t))).*pF(t)-3.*(-1+cmean).^2.*(-1+2.*cmean.*pC(t)).*pF(t).^2+(-1+cmean).^3.*pF(t).^3;
    
    %  6 PBSL [1/-3]    loss-averse
    c(6,t)=1+cmean.^2.*(-1+pC(t)).*(3-3.*pC(t)+cmean.*(-2+pC(t).*(4+pC(t))))-6.*(-1+cmean).*cmean.*(1+cmean.*(-1+pC(t))).*(-1+pC(t)).*pF(t)-3.*(-1+cmean).^2.*cmean.*(-1+2.*pC(t)).*pF(t).^2+(-1+cmean).^3.*pF(t).^3;
    
    %  7 PBSL [1/0]     only gains
    c(7,t)=(1+cmean.*pC(t).*(3+cmean.*pC(t).*(-3+cmean.*pC(t)))-3.*(-1+cmean).*(-1+cmean.^2.*pC(t).^2).*pF(t)-3.*(-1+cmean).^2.*(-1+cmean.*pC(t)).*pF(t).^2+(-1+cmean).^3.*pF(t).^3)/2;
    
    %  8 PBSL [0/-1]    only losses
    c(8,t)=(2+6.*cmean.^2.*(-1+pC(t))+cmean.^3.*(4-6.*pC(t)+pC(t).^3)-3.*(-1+cmean).*cmean.*(2.*(-1+pC(t))+cmean.*(2+(-4+pC(t)).*pC(t))).*pF(t)-3.*(-1+cmean).^2.*cmean.*pC(t).*pF(t).^2+(-1+cmean).^3.*pF(t).^3)/2;
    
    %  9 PBSL [3/1]     pbsl-conformist
    c(9,t)=cmean.*(3.*(-1+cmean).*cmean.*pC(t).^2.*pF(t)+cmean.*(3-2.*cmean+3.*(-1+cmean).*pF(t))+3.*(-1+cmean).*pC(t).*(-1+cmean+(2+cmean.*(-4+pF(t))-pF(t)).*pF(t)));
    
    % 10 PBSL [-1/-3]   pbsl-contrarian, hates wins, hates losses even more
    c(10,t)=(-1+cmean).*(-1+cmean.*(-1+cmean.*(2-3.*pC(t).^2)+3.*cmean.*pC(t).^2.*pF(t)+3.*(-1+cmean).*(-1+pC(t)).*pF(t).^2));
    
    % 11 PBSL a la McElreath
    c(11,t)=cmean.*(-3.*(-1+cmean).*cmean.*pC(t).^2.*pF(t)+cmean.*(3-2.*cmean+3.*(-1+cmean).*pF(t))-3.*(-1+cmean).^2.*pC(t).*(-1+pF(t).^2));
    

end

% calculation of PERFORMANCE
c2=c-.5;
pdif=sign(repmat(pC-pF,nStrats,1));
ctimesp=c2.*pdif;
ctimesp2=ctimesp;
ctimesp2(find(pdif==0))=[];
ctimesp2=reshape(ctimesp2,nStrats,length(ctimesp2)/nStrats);
% performance
perf=mean(ctimesp,2)+.5;
perf2=mean(ctimesp2,2)+.5;
% standard deviaton of performance
perfstd=std((ctimesp)')';
% standard error
perferr=perfstd*1.96/sqrt(tmax);

toc

[[1:nStrats]' perf2 perferr]
