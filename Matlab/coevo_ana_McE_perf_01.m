% this programm just computes the behavior of the strategies
% without any kind of evolution.
% The goal is to find out the PERFORMANCE of the different strategies
% This particular instance tests PBSLs McElreath for
% sample sizes 2 to 7.

clear all

'random seed'
rand('seed',17413)

% PARAMETERS
tmax=100000                    % periods per generation
regime=1;                   % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;                 % increment at which the environment becomes better or worse
pincr=1;                    % probability that environmental quality changes at all after each period
dpA=0
dpB=0                    % shift of pA and pB
pA0=.5;pB0=.5;
indLearn=0;                 % if == 0, ind. learning is win-stay lose-shift, else its constant
choiceLimit=0.001          % puts a limit on the min/max proportion, so that they are not 0 or 1

x=zeros(6,tmax);        % percent of A choices
perfmat=x;              % performance matrix
    
% ROUTINE
    
% ENVIRONMENT
% routine to determine pA and pB
[pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);
% possible shift dp
pA=min(1,max(0,pA+dpA));pB=min(1,max(pB+dpB,0));


tic

% first choice random
x(:,:,1)=0.5;

wb=waitbar(0,'progress');

% ROUTINE
for t=1:tmax-1
    x(1,t+1)=x(1,t)*(1+pA(t)+pB(t)*(-1+x(1,t))-pA(t)*x(1,t));
    x(2,t+1)=x(2,t)*(-3*pA(t)*(-1+pB(t).^2)*(-1+x(2,t)).^2+(3+3*pB(t)*(-1+x(2,t))-2*x(2,t))*x(2,t)-3*pA(t).^2*pB(t)*(-1+x(2,t))*x(2,t));
    x(3,t+1)=x(3,t)*((1+pB(t)*(-1+x(3,t)))*(3+3*pB(t)*(-1+x(3,t))-2*x(3,t))*x(3,t)+3*pA(t).^2*(-1+2*pB(t))*(-1+x(3,t)).^2*x(3,t)-4*pA(t).^3*pB(t)*(-1+x(3,t))*x(3,t).^2+2*pA(t)*(-1+x(3,t)).^2*(2+2*pB(t).^3*(-1+x(3,t))+x(3,t)-3*pB(t).^2*x(3,t)));
    x(4,t+1)=x(4,t)*(10*pA(t).^3*pB(t)*(-4+5*pB(t))*(-1+x(4,t)).^2*x(4,t).^2-5*pA(t).^4*pB(t)*(-1+x(4,t))*x(4,t).^3+10*pA(t).^2*(-1+pB(t))*(-1+x(4,t)).^2*x(4,t)*(1+pB(t)+5*pB(t).^2*(-1+x(4,t))-x(4,t)-7*pB(t)*x(4,t))-5*pA(t)*(-1+x(4,t)).^3*(1-pB(t).^4+(3+pB(t).^2*(-12+pB(t)*(8+pB(t))))*x(4,t))+x(4,t).^2*(10+10*pB(t).^2*(-1+x(4,t)).^2+3*x(4,t)*(-5+2*x(4,t))-5*pB(t)*(-1+x(4,t))*(-4+3*x(4,t))));
    x(5,t+1)=x(5,t)*(-15*pA(t).^4*pB(t)*(-6+5*pB(t))*(-1+x(5,t)).^2*x(5,t).^3-6*pA(t).^5*pB(t)*(-1+x(5,t))*x(5,t).^4+6*pA(t)*(-1+pB(t))*(-1+x(5,t)).^3*(1+pB(t)+pB(t).^2+pB(t).^3+pB(t).^4+(3+pB(t)*(3+pB(t)*(-27+(13-2*pB(t))*pB(t))))*x(5,t)+(1+pB(t)*(1+pB(t)*(16+(-14+pB(t))*pB(t))))*x(5,t).^2)+(1+pB(t)*(-1+x(5,t)))*x(5,t).^2*(10+10*pB(t).^2*(-1+x(5,t)).^2+3*x(5,t)*(-5+2*x(5,t))-5*pB(t)*(-1+x(5,t))*(-4+3*x(5,t)))+10*pA(t).^3*(-1+x(5,t)).^2*x(5,t).^2*(1-x(5,t)+6*(-1+pB(t))*pB(t)*(1+3*x(5,t)))+15*pA(t).^2*(-1+pB(t))*(-1+x(5,t)).^2*x(5,t)*(1+5*pB(t).^3*(-1+x(5,t)).^2-x(5,t).^2-pB(t).^2*(-1+x(5,t))*(-11+7*x(5,t))-pB(t)*(1+x(5,t))*(-1+7*x(5,t))));
    x(6,t+1)=x(6,t)*(-21*pA(t).^5*pB(t)*(-12+11*pB(t))*(-1+x(6,t)).^2*x(6,t).^4-7*pA(t).^6*pB(t)*(-1+x(6,t))*x(6,t).^5-35*pA(t).^4*pB(t)*(-1+x(6,t)).^2*x(6,t).^3*(-9+pB(t)*(27+19*pB(t)*(-1+x(6,t))-45*x(6,t))+27*x(6,t))-7*pA(t)*(-1+x(6,t)).^4*(-1+pB(t).^6-2*(2+pB(t).^3*(-30+(-15+pB(t))*(-3+pB(t))*pB(t)))*x(6,t)+(-10+pB(t).^2*(90+pB(t)*(-180+pB(t)*(135+(-36+pB(t))*pB(t)))))*x(6,t).^2)+x(6,t).^3*(35+35*pB(t).^3*(-1+x(6,t)).^3-21*pB(t).^2*(-1+x(6,t)).^2*(-5+4*x(6,t))-2*x(6,t)*(42+5*x(6,t)*(-7+2*x(6,t)))+7*pB(t)*(-1+x(6,t))*(15+2*x(6,t)*(-12+5*x(6,t))))-35*pA(t).^3*(-1+pB(t))*(-1+x(6,t)).^2*x(6,t).^2*((-1+x(6,t)).^2+19*pB(t).^3*(-1+x(6,t)).^2-pB(t).^2*(-1+x(6,t))*(-17+53*x(6,t))+pB(t)*(1+x(6,t)*(-26+37*x(6,t))))-21*pA(t).^2*(-1+pB(t))*(-1+x(6,t)).^3*x(6,t)*(1+11*pB(t).^4*(-1+x(6,t)).^2+(3-4*x(6,t))*x(6,t)-pB(t).^3*(-1+x(6,t))*(-19+64*x(6,t))+pB(t)*(1+(3-34*x(6,t))*x(6,t))+pB(t).^2*(1+x(6,t)*(-57+86*x(6,t)))));
    
    x(:,t+1)=min(x(:,t+1),1-choiceLimit);
    x(:,t+1)=max(x(:,t+1),choiceLimit);
    
    if mod(t,10)==0
        waitbar(t/tmax);
    end
end

close(wb)

% calculation of PERFORMANCE
perfmat2=perfmat;
% for t=1:tmax
%     perfmat2(:,t)=pA(t).*x(:,t)+pB(t).*(1-x(:,t));
% end
for i=1:6
    perfmat2(i,:)=pA.*x(i,:)+pB.*(1-x(i,:));
end

ind=find(pA==pB); %indices of when pA and pB are equal (->draw);
delp=pA-pB; %delta p
delp(ind)=[];   %remove draws
delp=sign(delp);    % A or B better
x2=x;   %save complete form
x(:,ind)=[];  %remove draws from choices
perfmat(:,ind)=[];
for t=1:tmax-length(ind)
    perfmat(:,t)=(delp(t)>0).*x(:,t)+(delp(t)<0).*(1-x(:,t));
end

[[2:7]' mean(perfmat,2) mean(perfmat2,2)]

'best performance:'
max(mean(perfmat,2))
'sample size of that PBSL McElreath'
1+find(mean(perfmat,2)==max(mean(perfmat,2)))


toc