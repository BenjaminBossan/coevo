% this programm just computes the behavior of the strategies
% without any kind of evolution.
% The goal is to find out the PERFORMANCE of the different strategies
% This particular instance tests scoring PBSL with weights
% [1/0] (-> x) and weights [4/-1] (-> y)
clear all

'random seed'
rand('seed',17413)

% PARAMETERS
tmax=250                    % periods per generation
regime=1                   % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;                 % increment at which the environment becomes better or worse
pincr=1;                    % probability that environmental quality changes at all after each period
dpA=0;
dpB=0;                    % shift of pA and pB
pA0=1/2;pB0=1/2;
choiceLimit=1/10000        % puts a limit on the min/max proportion, so that they are not 0 or 1

% INITIALIZE
x=zeros(6,tmax);        % percent of A choices
y=x;
perfmatx=x;              % performance matrix
perfmaty=x;
    
% ROUTINE
    
% ENVIRONMENT
% routine to determine pA and pB
[pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);
% possible shift dp
pA=min(1,max(0,pA+dpA));pB=min(1,max(pB+dpB,0));


tic

% first choice random
x(:,:,1)=1/2;
y(:,:,1)=1/2;

wb=waitbar(0,'progress');

% ROUTINE
for t=1:tmax-1
    
    x(1,t+1)=(1+2*pB(t)*(-1+x(1,t))+pB(t)^2*(-1+x(1,t))^2+2*pA(t)*x(1,t)-pA(t)^2*x(1,t)^2)/2;
    x(2,t+1)=(1+pB(t)^3*(-1+x(2,t))^3+3*pA(t)*x(2,t)-3*pA(t)^2*x(2,t)^2+pA(t)^3*x(2,t)^3-3*pB(t)^2*(-1+x(2,t))^2*(-1+pA(t)*x(2,t))-3*pB(t)*(-1+x(2,t))*(-1+pA(t)^2*x(2,t)^2))/2;
    x(3,t+1)=(1+pB(t)^4*(-1+x(3,t))^4+4*pA(t)*x(3,t)-6*pA(t)^2*x(3,t)^2+4*pA(t)^3*x(3,t)^3-pA(t)^4*x(3,t)^4-6*pB(t)^2*(-1+x(3,t))^2*(-1+2*pA(t)*x(3,t))-4*pB(t)^3*(-1+x(3,t))^3*(-1+2*pA(t)*x(3,t))+4*pB(t)*(-1+x(3,t))*(-1+pA(t)*x(3,t))^2*(1+2*pA(t)*x(3,t)))/2;
    x(4,t+1)=(1+pB(t)^5*(-1+x(4,t))^5+5*pA(t)*x(4,t)-10*pA(t)^2*x(4,t)^2+10*pA(t)^3*x(4,t)^3-5*pA(t)^4*x(4,t)^4+pA(t)^5*x(4,t)^5-5*pB(t)^4*(-1+x(4,t))^4*(-1+3*pA(t)*x(4,t))-5*pB(t)*(-1+x(4,t))*(-1+pA(t)*x(4,t))^3*(1+3*pA(t)*x(4,t))+10*pB(t)^3*(-1+x(4,t))^3*(1-4*pA(t)*x(4,t)+2*pA(t)^2*x(4,t)^2)+10*pB(t)^2*(-1+x(4,t))^2*(1-3*pA(t)*x(4,t)+2*pA(t)^3*x(4,t)^3))/2;
    x(5,t+1)=(1+pB(t)^6*(-1+x(5,t))^6+6*pA(t)*x(5,t)-15*pA(t)^2*x(5,t)^2+20*pA(t)^3*x(5,t)^3-15*pA(t)^4*x(5,t)^4+6*pA(t)^5*x(5,t)^5-pA(t)^6*x(5,t)^6-6*pB(t)^5*(-1+x(5,t))^5*(-1+4*pA(t)*x(5,t))+6*pB(t)*(-1+x(5,t))*(-1+pA(t)*x(5,t))^4*(1+4*pA(t)*x(5,t))+15*pB(t)^4*(-1+x(5,t))^4*(1-6*pA(t)*x(5,t)+5*pA(t)^2*x(5,t)^2)-15*pB(t)^2*(-1+x(5,t))^2*(-1+pA(t)*x(5,t))^2*(-1+2*pA(t)*x(5,t)+5*pA(t)^2*x(5,t)^2)+20*pB(t)^3*(-1+x(5,t))^3*(1-6*pA(t)*x(5,t)+6*pA(t)^2*x(5,t)^2))/2;
    x(6,t+1)=(1+pB(t)^7*(-1+x(6,t))^7+7*pA(t)*x(6,t)-21*pA(t)^2*x(6,t)^2+35*pA(t)^3*x(6,t)^3-35*pA(t)^4*x(6,t)^4+21*pA(t)^5*x(6,t)^5-7*pA(t)^6*x(6,t)^6+pA(t)^7*x(6,t)^7-7*pB(t)^6*(-1+x(6,t))^6*(-1+5*pA(t)*x(6,t))-7*pB(t)*(-1+x(6,t))*(-1+pA(t)*x(6,t))^5*(1+5*pA(t)*x(6,t))+21*pB(t)^5*(-1+x(6,t))^5*(1-8*pA(t)*x(6,t)+9*pA(t)^2*x(6,t)^2)+21*pB(t)^2*(-1+x(6,t))^2*(-1+pA(t)*x(6,t))^3*(-1+2*pA(t)*x(6,t)+9*pA(t)^2*x(6,t)^2)-35*pB(t)^4*(-1+x(6,t))^4*(-1+9*pA(t)*x(6,t)-15*pA(t)^2*x(6,t)^2+5*pA(t)^3*x(6,t)^3)-35*pB(t)^3*(-1+x(6,t))^3*(-1+8*pA(t)*x(6,t)-12*pA(t)^2*x(6,t)^2+5*pA(t)^4*x(6,t)^4))/2;
    
    y(1,t+1)=1+pB(t)^2*(-1+y(1,t))^2+(-1+pA(t))*y(1,t)-(-1+pA(t))*pA(t)*y(1,t)^2-pB(t)*(2-3*y(1,t)+y(1,t)^2);
    y(2,t+1)=1+pB(t)^3*(-1+y(2,t))^3-3*(-1+pA(t))^2*y(2,t)^2+(2-3*pA(t)+pA(t)^3)*y(2,t)^3-3*pB(t)^2*(-1+y(2,t))^2*(-1+2*pA(t)*y(2,t))-3*pB(t)*(-1+y(2,t))*(-1+2*pA(t)*(1-2*y(2,t))*y(2,t)+y(2,t)^2+2*pA(t)^2*y(2,t)^2);
    y(3,t+1)=1+pB(t)^4*(-1+y(3,t))^4-3*(-1+pA(t))^2*y(3,t)^2+2*(-1+pA(t))^2*(1+2*pA(t))*y(3,t)^3-pA(t)*(2-3*pA(t)+pA(t)^3)*y(3,t)^4-4*pB(t)^3*(-1+y(3,t))^3*(-1+3*pA(t)*y(3,t))+3*pB(t)^2*(-1+y(3,t))^2*(2-y(3,t)^2+2*pA(t)*y(3,t)*(-4+3*y(3,t)))+2*pB(t)*(-1+y(3,t))*(2-3*y(3,t)^2+y(3,t)^3+6*pA(t)^3*y(3,t)^3+6*pA(t)*y(3,t)*(-1+2*y(3,t))-3*pA(t)^2*y(3,t)^2*(1+3*y(3,t)));
    y(4,t+1)=1-(3*pB(t)^5*(-1+y(4,t))^5)/2+10*(-1+pA(t))^3*y(4,t)^3-5*(-1+pA(t))^3*(3+pA(t))*y(4,t)^4-((-1+pA(t))^3*(-12-11*pA(t)+3*pA(t)^2)*y(4,t)^5)/2-5*pB(t)^4*(-1+y(4,t))^4*(1+(-2+4*pA(t))*y(4,t))+5*pB(t)^3*(-1+y(4,t))^3*(-1+(6-12*pA(t))*y(4,t)+3*(-1+4*pA(t)^2)*y(4,t)^2)-(5*pB(t)*(-1+y(4,t))*(-1+(-4+8*pA(t))*y(4,t)+(6-12*pA(t)^2)*y(4,t)^2+(4-48*pA(t)+72*pA(t)^2-24*pA(t)^3)*y(4,t)^3+(-5+32*pA(t)-36*pA(t)^2+8*pA(t)^4)*y(4,t)^4))/2+30*pB(t)^2*(-1+y(4,t))^2*y(4,t)*(1-y(4,t)+3*pA(t)^2*(1-2*y(4,t))*y(4,t)+2*pA(t)^3*y(4,t)^2+pA(t)*(-2+3*y(4,t)^2));
    y(5,t+1)=1-5*pB(t)^6*(-1+y(5,t))^6+10*(-1+pA(t))^3*y(5,t)^3-15*(-1+pA(t))^3*(1+pA(t))*y(5,t)^4+6*(-1+pA(t))^3*(1+3*pA(t)+pA(t)^2)*y(5,t)^5+5*(-3+pA(t))*(-1+pA(t))^3*pA(t)^2*y(5,t)^6-6*pB(t)^5*(-1+y(5,t))^5*(4+5*(-1+pA(t))*y(5,t))+15*pB(t)^4*(-1+y(5,t))^4*(-3-8*(-1+pA(t))*y(5,t)+2*(-2+5*pA(t)^2)*y(5,t)^2)-10*pB(t)^3*(-1+y(5,t))^3*(4+18*(-1+pA(t))*y(5,t)+(18-36*pA(t)^2)*y(5,t)^2+(-5-12*pA(t)+30*pA(t)^2)*y(5,t)^3)+30*(-1+pA(t))*pB(t)*(-1+y(5,t))*y(5,t)*(-1+2*(1+pA(t))*y(5,t)+(-1-7*pA(t)+2*pA(t)^2)*y(5,t)^2+4*pA(t)*(1+pA(t)-pA(t)^2)*y(5,t)^3+pA(t)^2*(-3+pA(t)+pA(t)^2)*y(5,t)^4)-15*pB(t)^2*(-1+y(5,t))^2*(1+8*(-1+pA(t))*y(5,t)-6*(-2+3*pA(t)^2)*y(5,t)^2-2*(3+9*pA(t)-18*pA(t)^2+2*pA(t)^3)*y(5,t)^3+(1+6*pA(t)-20*pA(t)^3+10*pA(t)^4)*y(5,t)^4);
    y(6,t+1)=1-6*pB(t)^7*(-1+y(6,t))^7-35*(-1+pA(t))^4*y(6,t)^4+21*(-1+pA(t))^4*(4+pA(t))*y(6,t)^5+7*(-1+pA(t))^4*(-10-7*pA(t)+2*pA(t)^2)*y(6,t)^6-(-1+pA(t))^4*(-20-24*pA(t)+3*pA(t)^2+6*pA(t)^3)*y(6,t)^7+(7*pB(t)^6*(-1+y(6,t))^6*(-10+(6+9*pA(t))*y(6,t)))/2+21*pB(t)^5*(-1+y(6,t))^5*(-4+5*(1+pA(t))*y(6,t)+15*(-1+pA(t))*pA(t)*y(6,t)^2)-35*pB(t)^4*(-1+y(6,t))^4*(3-3*(2+pA(t))*y(6,t)-30*(-1+pA(t))*pA(t)*y(6,t)^2+(2-15*pA(t)+20*pA(t)^3)*y(6,t)^3)-70*pB(t)^3*(-1+y(6,t))^3*(1-3*y(6,t)-18*(-1+pA(t))*pA(t)*y(6,t)^2+(3-18*pA(t)+20*pA(t)^3)*y(6,t)^3+(-1-2*pA(t)+30*pA(t)^2-40*pA(t)^3+10*pA(t)^4)*y(6,t)^4)+(21*pB(t)^2*(-1+y(6,t))^2*(-2-5*(-2+pA(t))*y(6,t)+60*(-1+pA(t))*pA(t)*y(6,t)^2-10*(2-9*pA(t)+8*pA(t)^3)*y(6,t)^3-10*(-1-6*pA(t)+36*pA(t)^2-40*pA(t)^3+10*pA(t)^4)*y(6,t)^4+(2-75*pA(t)+240*pA(t)^2-200*pA(t)^3+30*pA(t)^5)*y(6,t)^5))/2+(7*(-1+pA(t))*pB(t)*(-1+y(6,t))*y(6,t)*(3*pA(t)^4*(20-27*y(6,t))*y(6,t)^4+9*pA(t)^5*y(6,t)^5+2*(-1+y(6,t))^3*(3+9*y(6,t)+8*y(6,t)^2)+3*pA(t)^3*y(6,t)^3*(-40+20*y(6,t)+23*y(6,t)^2)+pA(t)^2*y(6,t)^2*(-40+360*y(6,t)-420*y(6,t)^2+109*y(6,t)^3)-2*pA(t)*y(6,t)*(-15+20*y(6,t)+90*y(6,t)^2-150*y(6,t)^3+58*y(6,t)^4)))/2;
    
    x(:,t+1)=min(x(:,t+1),1-choiceLimit);
    x(:,t+1)=max(x(:,t+1),choiceLimit);
    
    y(:,t+1)=min(y(:,t+1),1-choiceLimit);
    y(:,t+1)=max(y(:,t+1),choiceLimit);
    
    if mod(t,10)==0
        waitbar(t/tmax);
    end
end

close(wb)

% calculation of PERFORMANCE
ind=find(pA==pB); %indices of when pA and pB are equal (->draw);
delp=pA-pB; %delta p
delp(ind)=[];   %remove draws
delp=sign(delp);    % A or B better
x2=x;   %save complete form
y2=y;   %save complete form
x(:,ind)=[];  %remove draws from choices
y(:,ind)=[];  %remove draws from choices
perfmatx2=perfmatx; %save complete form
perfmaty2=perfmaty; %save complete form
perfmatx(:,ind)=[];
perfmaty(:,ind)=[];
for t=1:tmax-length(ind)
    perfmatx(:,t)=repmat((delp(t)>0),6,1).*x(:,t)+repmat((delp(t)<0),6,1).*(1-x(:,t));
    perfmaty(:,t)=repmat((delp(t)>0),6,1).*y(:,t)+repmat((delp(t)<0),6,1).*(1-y(:,t));
end

'weights [1/0] -- [4/-1]'
[[2:7]' mean(perfmatx,2) [2:7]' mean(perfmaty,2)]

'best performance for weights [1/0] -- [4/-1]'
[max(mean(perfmatx,2)) 1+find(mean(perfmatx,2)==max(mean(perfmatx,2))) max(mean(perfmaty,2)) 1+find(mean(perfmaty,2)==max(mean(perfmaty,2)))]



toc