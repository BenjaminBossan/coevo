% this programm just computes the behavior of the strategies
% without any kind of evolution.
% The goal is to find out the PERFORMANCE of the different strategies
% This particular instance tests out all pbsl score type strategies with
% weights ranging from - to + WEIGHTSCOPE (including redundant ones)
% it is assumed that we have homogeneous populations
% The sample size of the social learners can either be 3 or 7.
% Note that with 7, computations will be slowed down by several
% orders of magnitude.

clear all

% PARAMETERS
tmax=20                    % periods per generation
regime=1;                   % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;                 % increment at which the environment becomes better or worse
pincr=1;                    % probability that environmental quality changes at all after each period
dpA=0;
dpB=0;                     % shift of pA and pB
pA0=.5;pB0=.5;
weightScope=7              % the min and max weight that is tested
choiceLimit=0.001;          % puts a limit on the min/max proportion, so that they are not 0 or 1
tallyn=7                    % tally number (=sample size)

x=zeros(2*weightScope+1,(2*weightScope+1)*tmax);      % percent of A choices
x=reshape(x,2*weightScope+1,2*weightScope+1,tmax);      %tensor of A choices over time
perfmat=x;  % performance matrix
perfInd=zeros(1,tmax);
perfMcE=zeros(1,tmax);

% individual learners & PBSL McElreath
xind=zeros(1,tmax);
xMcE=zeros(1,tmax);
    
% ROUTINE
    
% ENVIRONMENT
% routine to determine pA and pB
[pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);
% possible shift dp
pA=min(1,max(0,pA+dpA));pB=min(1,max(pB+dpB,0));

% 'CARE HERE'
% [pA pB]=randomenvironment_AR_01(tmax,pA0,pB0,0);

tic

%     % individual learners
%     if indLearn==0
%         cind=(1-pF(t))/(2-pC(t)-pF(t));
%     else
%         cind=((1-2*(1-indLearn))*sign(pC(t)-pF(t))+1)/2;
%     end
%     c(1,t)=cind;
    
% first choice random
x(:,:,1)=0.5;

wb=waitbar(0,'progress');
counter=0;
counterMax=(2*weightScope+1)^2;

% ROUTINE
for i=1:2*weightScope+1
    for j=1:2*weightScope+1
        wG=i-weightScope-1; %weight gains
        wL=j-weightScope-1; %weight losses
        
        if tallyn==7
            for t=1:tmax-1
                x(i,j,t+1)=(pB(t).^7*(-1+x(i,j,t)).^7*(-1+sign(wG))-42*(-1+pA(t))*(-1+pB(t))*pB(t).^5*(-1+x(i,j,t)).^6*x(i,j,t)*(-1+sign(wG))-7*pA(t)*pB(t).^6*(-1+x(i,j,t)).^6*x(i,j,t)*(-1+sign(wG))+210*(-1+pA(t)).^2*(-1+pB(t)).^2*pB(t).^3*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wG))+210*(-1+pA(t))*pA(t)*(-1+pB(t))*pB(t).^4*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wG))+21*pA(t).^2*pB(t).^5*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wG))-140*(-1+pA(t)).^3*(-1+pB(t)).^3*pB(t)*(-1+x(i,j,t)).^4*x(i,j,t).^3*(-1+sign(wG))-630*(-1+pA(t)).^2*pA(t)*(-1+pB(t)).^2*pB(t).^2*(-1+x(i,j,t)).^4*x(i,j,t).^3*(-1+sign(wG))-420*(-1+pA(t))*pA(t).^2*(-1+pB(t))*pB(t).^3*(-1+x(i,j,t)).^4*x(i,j,t).^3*(-1+sign(wG))-35*pA(t).^3*pB(t).^4*(-1+x(i,j,t)).^4*x(i,j,t).^3*(-1+sign(wG))-140*(-1+pA(t)).^3*pA(t)*(-1+pB(t)).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(wG))-630*(-1+pA(t)).^2*pA(t).^2*(-1+pB(t)).^2*pB(t)*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(wG))-420*(-1+pA(t))*pA(t).^3*(-1+pB(t))*pB(t).^2*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(wG))-35*pA(t).^4*pB(t).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(wG))+210*(-1+pA(t)).^2*pA(t).^3*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wG))+210*(-1+pA(t))*pA(t).^4*(-1+pB(t))*pB(t)*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wG))+21*pA(t).^5*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wG))-42*(-1+pA(t))*pA(t).^5*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(wG))-7*pA(t).^6*pB(t)*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(wG))+pA(t).^7*x(i,j,t).^7*(1+sign(wG))+7*pA(t)*(-1+pB(t)).^6*(-1+x(i,j,t)).^6*x(i,j,t)*(1+sign(wG-6*wL))+21*(-1+pB(t)).^5*pB(t).^2*(-1+x(i,j,t)).^7*(1+sign(-2*wG-5*wL))+21*pA(t).^2*(-1+pB(t)).^5*(-1+x(i,j,t)).^5*x(i,j,t).^2*(1+sign(2*wG-5*wL))-35*(-1+pB(t)).^4*pB(t).^3*(-1+x(i,j,t)).^7*(1+sign(-3*wG-4*wL))-42*(-1+pA(t))*pA(t)*(-1+pB(t)).^5*(-1+x(i,j,t)).^5*x(i,j,t).^2*(1+sign(wG-4*wL))-105*pA(t).^2*(-1+pB(t)).^4*pB(t)*(-1+x(i,j,t)).^5*x(i,j,t).^2*(1+sign(wG-4*wL))+35*pA(t).^3*(-1+pB(t)).^4*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(3*wG-4*wL))+35*(-1+pB(t)).^3*pB(t).^4*(-1+x(i,j,t)).^7*(1+sign(-4*wG-3*wL))-105*(-1+pA(t))*(-1+pB(t)).^4*pB(t).^2*(-1+x(i,j,t)).^6*x(i,j,t)*(1+sign(-2*wG-3*wL))-140*pA(t)*(-1+pB(t)).^3*pB(t).^3*(-1+x(i,j,t)).^6*x(i,j,t)*(1+sign(-2*wG-3*wL))-105*(-1+pA(t))*pA(t).^2*(-1+pB(t)).^4*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(2*wG-3*wL))-140*pA(t).^3*(-1+pB(t)).^3*pB(t)*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(2*wG-3*wL))+35*pA(t).^4*(-1+pB(t)).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(4*wG-3*wL))-21*(-1+pB(t)).^2*pB(t).^5*(-1+x(i,j,t)).^7*(1+sign(-5*wG-2*wL))+140*(-1+pA(t))*(-1+pB(t)).^3*pB(t).^3*(-1+x(i,j,t)).^6*x(i,j,t)*(1+sign(-3*wG-2*wL))+105*pA(t)*(-1+pB(t)).^2*pB(t).^4*(-1+x(i,j,t)).^6*x(i,j,t)*(1+sign(-3*wG-2*wL))+105*(-1+pA(t)).^2*pA(t)*(-1+pB(t)).^4*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(wG-2*wL))+420*(-1+pA(t))*pA(t).^2*(-1+pB(t)).^3*pB(t)*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(wG-2*wL))+210*pA(t).^3*(-1+pB(t)).^2*pB(t).^2*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(wG-2*wL))-140*(-1+pA(t))*pA(t).^3*(-1+pB(t)).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(3*wG-2*wL))-105*pA(t).^4*(-1+pB(t)).^2*pB(t)*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(3*wG-2*wL))+21*pA(t).^5*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(5*wG-2*wL))+210*(-1+pA(t)).^2*pA(t).^2*(-1+pB(t)).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(2*wG-wL))+420*(-1+pA(t))*pA(t).^3*(-1+pB(t)).^2*pB(t)*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(2*wG-wL))+105*pA(t).^4*(-1+pB(t))*pB(t).^2*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(2*wG-wL))-105*(-1+pA(t))*pA(t).^4*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(4*wG-wL))-42*pA(t).^5*(-1+pB(t))*pB(t)*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(4*wG-wL))+7*pA(t).^6*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(6*wG-wL))-(-1+pB(t)).^7*(-1+x(i,j,t)).^7*(-1+sign(wL))+7*(-1+pA(t))*(-1+pB(t)).^6*(-1+x(i,j,t)).^6*x(i,j,t)*(-1+sign(wL))+42*pA(t)*(-1+pB(t)).^5*pB(t)*(-1+x(i,j,t)).^6*x(i,j,t)*(-1+sign(wL))-21*(-1+pA(t)).^2*(-1+pB(t)).^5*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wL))-210*(-1+pA(t))*pA(t)*(-1+pB(t)).^4*pB(t)*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wL))-210*pA(t).^2*(-1+pB(t)).^3*pB(t).^2*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wL))+35*(-1+pA(t)).^3*(-1+pB(t)).^4*(-1+x(i,j,t)).^4*x(i,j,t).^3*(-1+sign(wL))+420*(-1+pA(t)).^2*pA(t)*(-1+pB(t)).^3*pB(t)*(-1+x(i,j,t)).^4*x(i,j,t).^3*(-1+sign(wL))+630*(-1+pA(t))*pA(t).^2*(-1+pB(t)).^2*pB(t).^2*(-1+x(i,j,t)).^4*x(i,j,t).^3*(-1+sign(wL))+140*pA(t).^3*(-1+pB(t))*pB(t).^3*(-1+x(i,j,t)).^4*x(i,j,t).^3*(-1+sign(wL))+35*(-1+pA(t)).^4*(-1+pB(t)).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(wL))+420*(-1+pA(t)).^3*pA(t)*(-1+pB(t)).^2*pB(t)*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(wL))+630*(-1+pA(t)).^2*pA(t).^2*(-1+pB(t))*pB(t).^2*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(wL))+140*(-1+pA(t))*pA(t).^3*pB(t).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(wL))-21*(-1+pA(t)).^5*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wL))-210*(-1+pA(t)).^4*pA(t)*(-1+pB(t))*pB(t)*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wL))-210*(-1+pA(t)).^3*pA(t).^2*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wL))+7*(-1+pA(t)).^6*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(wL))+42*(-1+pA(t)).^5*pA(t)*pB(t)*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(wL))-(-1+pA(t)).^7*x(i,j,t).^7*(1+sign(wL))-7*(-1+pA(t))*pB(t).^6*(-1+x(i,j,t)).^6*x(i,j,t)*(1+sign(-6*wG+wL))+105*(-1+pA(t)).^2*(-1+pB(t))*pB(t).^4*(-1+x(i,j,t)).^5*x(i,j,t).^2*(1+sign(-4*wG+wL))+42*(-1+pA(t))*pA(t)*pB(t).^5*(-1+x(i,j,t)).^5*x(i,j,t).^2*(1+sign(-4*wG+wL))-210*(-1+pA(t)).^3*(-1+pB(t)).^2*pB(t).^2*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(-2*wG+wL))-420*(-1+pA(t)).^2*pA(t)*(-1+pB(t))*pB(t).^3*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(-2*wG+wL))-105*(-1+pA(t))*pA(t).^2*pB(t).^4*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(-2*wG+wL))-210*(-1+pA(t)).^2*(-1+pB(t)).^3*pB(t).^2*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(2*wG+wL))-420*(-1+pA(t))*pA(t)*(-1+pB(t)).^2*pB(t).^3*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(2*wG+wL))-105*pA(t).^2*(-1+pB(t))*pB(t).^4*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(2*wG+wL))-210*(-1+pA(t)).^3*pA(t).^2*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(2*wG+wL))-420*(-1+pA(t)).^2*pA(t).^3*(-1+pB(t))*pB(t)*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(2*wG+wL))-105*(-1+pA(t))*pA(t).^4*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(2*wG+wL))+105*(-1+pA(t))*(-1+pB(t)).^2*pB(t).^4*(-1+x(i,j,t)).^6*x(i,j,t)*(-1+sign(4*wG+wL))+42*pA(t)*(-1+pB(t))*pB(t).^5*(-1+x(i,j,t)).^6*x(i,j,t)*(-1+sign(4*wG+wL))+105*(-1+pA(t)).^2*pA(t).^4*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(4*wG+wL))+42*(-1+pA(t))*pA(t).^5*pB(t)*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(4*wG+wL))-7*(-1+pB(t))*pB(t).^6*(-1+x(i,j,t)).^7*(-1+sign(6*wG+wL))-7*(-1+pA(t))*pA(t).^6*x(i,j,t).^7*(1+sign(6*wG+wL))-21*(-1+pA(t)).^2*pB(t).^5*(-1+x(i,j,t)).^5*x(i,j,t).^2*(1+sign(-5*wG+2*wL))+140*(-1+pA(t)).^3*(-1+pB(t))*pB(t).^3*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(-3*wG+2*wL))+105*(-1+pA(t)).^2*pA(t)*pB(t).^4*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(-3*wG+2*wL))-105*(-1+pA(t)).^4*(-1+pB(t)).^2*pB(t)*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(-wG+2*wL))-420*(-1+pA(t)).^3*pA(t)*(-1+pB(t))*pB(t).^2*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(-wG+2*wL))-210*(-1+pA(t)).^2*pA(t).^2*pB(t).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(-wG+2*wL))+105*(-1+pA(t)).^2*(-1+pB(t)).^4*pB(t)*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wG+2*wL))+420*(-1+pA(t))*pA(t)*(-1+pB(t)).^3*pB(t).^2*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wG+2*wL))+210*pA(t).^2*(-1+pB(t)).^2*pB(t).^3*(-1+x(i,j,t)).^5*x(i,j,t).^2*(-1+sign(wG+2*wL))+105*(-1+pA(t)).^4*pA(t)*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wG+2*wL))+420*(-1+pA(t)).^3*pA(t).^2*(-1+pB(t))*pB(t)*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wG+2*wL))+210*(-1+pA(t)).^2*pA(t).^3*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(wG+2*wL))-140*(-1+pA(t)).^3*pA(t).^3*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(3*wG+2*wL))-105*(-1+pA(t)).^2*pA(t).^4*pB(t)*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(3*wG+2*wL))+21*(-1+pA(t)).^2*pA(t).^5*x(i,j,t).^7*(1+sign(5*wG+2*wL))-35*(-1+pA(t)).^3*pB(t).^4*(-1+x(i,j,t)).^4*x(i,j,t).^3*(1+sign(-4*wG+3*wL))+105*(-1+pA(t)).^4*(-1+pB(t))*pB(t).^2*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(-2*wG+3*wL))+140*(-1+pA(t)).^3*pA(t)*pB(t).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(-2*wG+3*wL))+105*(-1+pA(t)).^4*pA(t).^2*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(2*wG+3*wL))+140*(-1+pA(t)).^3*pA(t).^3*pB(t)*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(2*wG+3*wL))-35*(-1+pA(t)).^3*pA(t).^4*x(i,j,t).^7*(1+sign(4*wG+3*wL))-35*(-1+pA(t)).^4*pB(t).^3*(-1+x(i,j,t)).^3*x(i,j,t).^4*(1+sign(-3*wG+4*wL))+42*(-1+pA(t)).^5*(-1+pB(t))*pB(t)*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(-wG+4*wL))+105*(-1+pA(t)).^4*pA(t)*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(-wG+4*wL))-42*(-1+pA(t))*(-1+pB(t)).^5*pB(t)*(-1+x(i,j,t)).^6*x(i,j,t)*(-1+sign(wG+4*wL))-105*pA(t)*(-1+pB(t)).^4*pB(t).^2*(-1+x(i,j,t)).^6*x(i,j,t)*(-1+sign(wG+4*wL))-42*(-1+pA(t)).^5*pA(t)*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(wG+4*wL))-105*(-1+pA(t)).^4*pA(t).^2*pB(t)*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(wG+4*wL))+35*(-1+pA(t)).^4*pA(t).^3*x(i,j,t).^7*(1+sign(3*wG+4*wL))-21*(-1+pA(t)).^5*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t).^5*(1+sign(-2*wG+5*wL))-21*(-1+pA(t)).^5*pA(t).^2*x(i,j,t).^7*(1+sign(2*wG+5*wL))-7*(-1+pA(t)).^6*pB(t)*(-1+x(i,j,t))*x(i,j,t).^6*(1+sign(-wG+6*wL))+7*(-1+pB(t)).^6*pB(t)*(-1+x(i,j,t)).^7*(-1+sign(wG+6*wL))+7*(-1+pA(t)).^6*pA(t)*x(i,j,t).^7*(1+sign(wG+6*wL)))/2;
                x(i,j,t+1)=min(x(i,j,t+1),1-choiceLimit);
                x(i,j,t+1)=max(x(i,j,t+1),choiceLimit);
            end
            
        elseif tallyn==3
        
            if wG>0&wL>0
                for t=1:tmax-1
                x(i,j,t+1)=(x(i,j,t)*(6*(-1+pA(t)).^2*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t)-12*(-1+pA(t))*pA(t)*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t)+12*(-1+pA(t))*pA(t)*pB(t)*(-1+x(i,j,t))*x(i,j,t)-6*pA(t).^2*pB(t)*(-1+x(i,j,t))*x(i,j,t)-2*(-1+pA(t)).^3*x(i,j,t).^2+6*(-1+pA(t)).^2*pA(t)*x(i,j,t).^2-6*(-1+pA(t))*pA(t).^2*x(i,j,t).^2+2*pA(t).^3*x(i,j,t).^2+3*pA(t)*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*(1+sign(wG-2*wL))+3*pA(t).^2*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t)*(1+sign(2*wG-wL))-3*(-1+pA(t))*pB(t).^2*(-1+x(i,j,t)).^2*(1+sign(-2*wG+wL))-3*(-1+pA(t)).^2*pB(t)*(-1+x(i,j,t))*x(i,j,t)*(1+sign(-wG+2*wL))))/2;
                x(i,j,t+1)=min(x(i,j,t+1),1-choiceLimit);
                x(i,j,t+1)=max(x(i,j,t+1),choiceLimit);
                end
            elseif wG>0&wL<0
                for t=1:tmax-1
                x(i,j,t+1)=(2*(-1+pB(t)).^3*(-1+x(i,j,t)).^3-6*(-1+pA(t))*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*x(i,j,t)+6*pA(t)*(-1+pB(t)).^2*(-1+x(i,j,t)).^2*x(i,j,t)-12*pA(t)*(-1+pB(t))*pB(t)*(-1+x(i,j,t)).^2*x(i,j,t)-12*(-1+pA(t))*pA(t)*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^2+6*pA(t).^2*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^2-6*pA(t).^2*pB(t)*(-1+x(i,j,t))*x(i,j,t).^2+2*pA(t).^3*x(i,j,t).^3-3*(-1+pB(t))*pB(t).^2*(-1+x(i,j,t)).^3*(-1+sign(2*wG+wL))-3*(-1+pA(t))*pA(t).^2*x(i,j,t).^3*(1+sign(2*wG+wL))+3*(-1+pB(t)).^2*pB(t)*(-1+x(i,j,t)).^3*(-1+sign(wG+2*wL))+3*(-1+pA(t)).^2*pA(t)*x(i,j,t).^3*(1+sign(wG+2*wL)))/2;
                x(i,j,t+1)=min(x(i,j,t+1),1-choiceLimit);
                x(i,j,t+1)=max(x(i,j,t+1),choiceLimit);
                end
            elseif wG<0&wL>0
                for t=1:tmax-1
                x(i,j,t+1)=(-2*pB(t).^3*(-1+x(i,j,t)).^3+12*(-1+pA(t))*(-1+pB(t))*pB(t)*(-1+x(i,j,t)).^2*x(i,j,t)-6*(-1+pA(t))*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t)+6*pA(t)*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t)+6*(-1+pA(t)).^2*(-1+pB(t))*(-1+x(i,j,t))*x(i,j,t).^2-6*(-1+pA(t)).^2*pB(t)*(-1+x(i,j,t))*x(i,j,t).^2+12*(-1+pA(t))*pA(t)*pB(t)*(-1+x(i,j,t))*x(i,j,t).^2-2*(-1+pA(t)).^3*x(i,j,t).^3-3*(-1+pB(t))*pB(t).^2*(-1+x(i,j,t)).^3*(-1+sign(2*wG+wL))-3*(-1+pA(t))*pA(t).^2*x(i,j,t).^3*(1+sign(2*wG+wL))+3*(-1+pB(t)).^2*pB(t)*(-1+x(i,j,t)).^3*(-1+sign(wG+2*wL))+3*(-1+pA(t)).^2*pA(t)*x(i,j,t).^3*(1+sign(wG+2*wL)))/2;
                x(i,j,t+1)=min(x(i,j,t+1),1-choiceLimit);
                x(i,j,t+1)=max(x(i,j,t+1),choiceLimit);
                end
            elseif wG<0&wL<0
                for t=1:tmax-1
                x(i,j,t+1)=((-1+x(i,j,t))*(-2-2*x(i,j,t)+3*pA(t)*x(i,j,t)-6*pA(t)*pB(t)*x(i,j,t)+3*pB(t).^2*x(i,j,t)+4*x(i,j,t).^2-3*pA(t)*x(i,j,t).^2-3*pA(t).^2*x(i,j,t).^2-3*pB(t)*x(i,j,t).^2+12*pA(t)*pB(t)*x(i,j,t).^2-3*pB(t).^2*x(i,j,t).^2+3*x(i,j,t)*(pB(t)*x(i,j,t)+pA(t).^2*pB(t)*x(i,j,t)+pA(t)*(-1+pB(t)*(2-4*x(i,j,t))+pB(t).^2*(-1+x(i,j,t))+x(i,j,t)))*sign(wG-2*wL)+3*pA(t).^2*(-1+pB(t))*x(i,j,t).^2*sign(2*wG-wL)-3*pB(t).^2*x(i,j,t)*sign(-2*wG+wL)+3*pA(t)*pB(t).^2*x(i,j,t)*sign(-2*wG+wL)+3*pB(t).^2*x(i,j,t).^2*sign(-2*wG+wL)-3*pA(t)*pB(t).^2*x(i,j,t).^2*sign(-2*wG+wL)))/2;
                x(i,j,t+1)=min(x(i,j,t+1),1-choiceLimit);
                x(i,j,t+1)=max(x(i,j,t+1),choiceLimit);
                end
            elseif wG==0
                for t=1:tmax-1
                x(i,j,t+1)=(1+(-1-pB(t).^3*(-1+x(i,j,t)).^3+3*pA(t)*pB(t).^2*(-1+x(i,j,t)).^2*x(i,j,t)-6*(-1+pA(t))*x(i,j,t).^2-(4-6*pA(t)+pA(t).^3)*x(i,j,t).^3+3*pB(t)*(-1+x(i,j,t))*x(i,j,t)*(-2+2*pA(t)+2*x(i,j,t)-4*pA(t)*x(i,j,t)+pA(t).^2*x(i,j,t)))*sign(wL))/2;
                x(i,j,t+1)=min(x(i,j,t+1),1-choiceLimit);
                x(i,j,t+1)=max(x(i,j,t+1),choiceLimit);
                end
            elseif wL==0
                for t=1:tmax-1
                x(i,j,t+1)=(1+(pB(t).^3*(-1+x(i,j,t)).^3-3*pB(t).^2*(-1+x(i,j,t)).^2*(-1+pA(t)*x(i,j,t))-3*pB(t)*(-1+x(i,j,t))*(-1+pA(t).^2*x(i,j,t).^2)+pA(t)*x(i,j,t)*(3-3*pA(t)*x(i,j,t)+pA(t).^2*x(i,j,t).^2))*sign(wG))/2;
                x(i,j,t+1)=min(x(i,j,t+1),1-choiceLimit);
                x(i,j,t+1)=max(x(i,j,t+1),choiceLimit);
                end
            end
        end
        
        counter=counter+1;
        waitbar(counter/counterMax);
    end
%     waitbar(i/(2*weightScope+1));
end

% for good measure, calculate individual learners and PBSL McElreath also
xind(1)=1/2;xMcE(1)=1/2;
for t=2:tmax
    xind(t)=(1-pB(t))/(2-pB(t)-pA(t));
    xMcE(t)=xMcE(t-1)*(-3*(-1+xMcE(t-1))*xMcE(t-1)*pA(t).^2*pB(t)+xMcE(t-1)*(3-2*xMcE(t-1)+3*(-1+xMcE(t-1))*pB(t))-3*(-1+xMcE(t-1)).^2*pA(t)*(-1+pB(t).^2));
    xMcE(t)=min(xMcE(t),1-choiceLimit);
    xMcE(t)=max(xMcE(t),choiceLimit);
end
close(wb)

% calculation of PERFORMANCE

perfmat2=perfmat;
perfInd2=perfInd;
perfMcE2=perfMcE;
x2=x;
for t=1:tmax
    perfmat2(:,:,t)=pA(t).*x(:,:,t)+pB(t).*(1-x(:,:,t));
end

perfInd2=pA.*xind+pB.*(1-xind);
perfMcE2=pA.*xMcE+pB.*(1-xMcE);

ind=find(pA==pB); %indices of when pA and pB are equal (->draw);
delp=pA-pB; %delta p
delp(ind)=[];   %remove draws
delp=sign(delp);    % A or B better
x(:,:,ind)=[];  %remove draws from choices
xind(ind)=[];
xMcE(ind)=[];
perfmat(:,:,ind)=[];
perfInd(ind)=[];
perfMcE(ind)=[];
for t=1:tmax-length(ind)
    perfmat(:,:,t)=(delp(t)>0).*x(:,:,t)+(delp(t)<0).*(1-x(:,:,t));
end

perfInd=(delp>0).*xind+(delp<0).*(1-xind);
perfMcE=(delp>0).*xMcE+(delp<0).*(1-xMcE);

[{'ind. learners'} {'pbsl McElreath'}]
[mean(perfInd) mean(perfMcE)]

'best performance:'
max(max(mean(perfmat,3)))
'weight gains, losses'
[s1 s2]=find(mean(perfmat,3)==max(max(mean(perfmat,3))));
[s1-weightScope-1 s2-weightScope-1]

x1=reshape(repmat([1:2*weightScope+1],2,1),2*(2*weightScope+1),1);
x1=(2*weightScope+1)*x1+.5;
y1=reshape(repmat([0 1 1 0]',weightScope+1,1),2*(2*weightScope+2),1);
y1(end)=[];
y1(end)=[];
figure
plot(x1,y1,'k')
hold on
pm1=reshape(mean(perfmat,3),(2*weightScope+1)^2,1);
plot(pm1,'k.')

fheatmap_01(mean(perfmat,3));
title('performance')

% % matrix that contains the degree of conformist bias
% % it is assumed that default parameters are used
% confTable=zeros(2*weightScope+1,2*weightScope+1);
% for gg=1:2*weightScope+1
%     wG=gg-weightScope-1;
%     for ll=1:2*weightScope+1
%         wL=ll-weightScope-1;
%         if dpA==0
%             confTable(gg,ll)=-0.2500000000000003+0.06291418750000007*sign(wG)-0.012029390625000025*sign(wG-2.*wL)+0.012029390624999999*sign(2*wG-wL)+0.06291418750000001*sign(wL)+0.08130642187499999*sign(2*wG+wL)+0.08130642187499998*sign(wG+2*wL);
%         elseif dpA==-.25
%             confTable(gg,ll)=-0.25+0.03145709374999996*sign(wG)-0.013546007812500009*sign(wG-2.*wL)+0.0046533984375000115*sign(2*wG-wL)+0.14124628125000005*sign(wL)+0.0311241328125*sign(2*wG+wL)+0.09047308593749998*sign(wG+2*wL);
%         elseif dpA==.25
%             confTable(gg,ll)=-0.2500000000000002+0.14124628125000005*sign(wG)-0.004653398437500018*sign(wG-2.*wL)+0.01354600781250001*sign(2*wG-wL)+0.031457093750000026*sign(wL)+0.09047308593749998*sign(2*wG+wL)+0.031124132812499976*sign(wG+2*wL);
%         end
%     end
% end

toc