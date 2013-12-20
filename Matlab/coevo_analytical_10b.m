% this is an attempt to make the model analytical. It is assumed that the
% performance of individual learners is constant. Path dependence is
% ignored. Social learners have 7 individuals as sample. Choice proportions
% are always in equilibrium.

clear all

% PARAMETERS
tmax=50;                    % periods per generation
gen=100000                % generations
selCoef=1               % maximum selection coefficient
nStrats=11;                 % total possible number of strategies
vertTransm=1;               % is there vertical tranmission from parents to offspring?
                            % 1->yes, 0->no, 0.5->yes but not for
                            % conformists
regime=1;                   % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;                 % increment at which the environment becomes better or worse
pincr=1;                    % probability that environmental quality changes at all after each period
dpA=0;dpB=0;                % shift of pA and pB
indLearn=0;                 % if == 0, ind. learning is win-stay lose-shift, else its constant
mutationRate=0/100         % probability of mutation occurring

tic

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

xInit=zeros(nStrats,1);

% xInit(1)=.99;
% xInit(2)=.01;
% xInit(3)=1;
% xInit(4)=1;
% xInit(5)=.5;
% xInit(6)=1;
% xInit(7)=1;
% xInit(8)=1;
% xInit(9)=1;
% xInit(10)=1;
% xInit(11)=.5;

xInit=ones(nStrats,1)/nStrats;

% TESTED PARAMETER RANGE
% incrVec=[.01:.01:.1];
% pincrVec=[.05:.05:1];
% incrVec=.02;
% pincrVec=1;
% incrVec=.02
% pincrVec=1
% incrVec=.02;
% pincrVec=1;
% regimeVec=[0 .1 .25 .5:.5:3.5];
regimeVec=[0:.1:3];
% regimeVec=1
% tmaxVec=[10 20 30 40 50 75 100 150 200:50:550];
% indVec=[.51:.02:.99];
% indVec=0
tmaxVec=50
% tmaxVec=[5 50 100 250];
% tmaxVec=50;
% dpAvec=.0
% dpAvec=[-.25:.01:.25];

% define the FOCAL PARAMETER VECTORS
% u2=dpBvec;
% u2=pincrVec;
% u1=incrVec;
% u1=dpAvec;
u1=regimeVec;
u2=tmaxVec;

xlab='t_{max}';
ylab='reversion facotr';

nRepit=length(u1)*length(u2);

m=zeros(nStrats,length(u1)*length(u2));
m=reshape(m,nStrats,length(u1),length(u2));
mv=m;   %variance

wb=waitbar(0,'progress');

for j=1:length(u1)
    
    % dpA=dpAvec(j);
    % dpB=dpA;
    % incr=u1(j);
    regime=regimeVec(j);
    
    for i=1:length(u2)
        
        tmax=tmaxVec(i);
        % pincr=u2(i);
        
        pA0=1/2;pB0=1/2;            % initial values of pA and pB, so that pA=/=pB forall t
            
        [x c w maxGen pA pB]=coevo_ana_10b(...
        xInit,...                   % initial conditions
        tmax,...                    % periods per generation
        gen,...                     % generations
        selCoef,...                 % maximum selection coefficient
        nStrats,...                 % total possible number of strategies
        vertTransm,...              % whether there is vertical transmission
        indLearn,...                % ind. learning performance
        mutationRate,...            % probability of mutation occurring
        regime,...                  % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
        pA0,...
        pB0,...                     % initial values of pA and pB
        incr,...                    % increment at which the environment becomes better or worse
        pincr,...                   % probability that environmental quality changes at all after each period
        dpA,...                     % shift of pA and pB
        dpB);

        xLast=mean(x(:,maxGen-9999:maxGen),2);
        xVar=var(x(:,maxGen-9999:maxGen)')';
        m(:,j,i)=xLast;
        mv(:,j,i)=xVar;
        
        if length(u1)>1
            waitbar(j/length(u1),wb);
        else
            waitbar(i/length(u2),wb);
        end
        
    end
end

% FIGURES


if nRepit>=6    
    
    for s=1:nStrats
%         % the focal strategy
%         mc=squeeze(m(s,:,:));
%         figure
%         % colormap from white->low freq to black->high freq
%         colormap([1-[0:.01:1]' 1-[0:.01:1]' 1-[0:.01:1]']);
%         % add 0 and 1 to data so that full scale [0;1] is used
%         mc2=[mc;zeros(1,length(u2))];
%         mc2(length(u1)+1,length(u2))=1;
%         lastEntry=2*u1(length(u1))-u1(length(u1)-1);
%         % plot "heat map"
%         surf(u2,[u1 lastEntry],zeros(length(u1)+1,length(u2)),mc2);
%         axis([min(u2) max(u2) min(u1) max(u1) 0 1])
%         set(gca,'fontsize',14)
%         xlabel(xlab,'fontsize',14)
%         ylabel(ylab,'fontsize',14)
%         % view from above
%         view([0 90])
%         axis square
%         colorbar
%         title(s)
        
        % the focal strategy
        mc=squeeze(m(s,:,:));
        figure
        % colormap from white->low freq to black->high freq
        colormap([1-[0:.01:1]' 1-[0:.01:1]' 1-[0:.01:1]']);
        mc=[mc;mc(end,:)];
        mc=[mc mc(:,end)];
        v1=[u1 u1(end)+1];v2=[u2 u2(end)+1];
        % add 0 and 1 to data so that full scale [0;1] is used
        mc2=[zeros(1,length(v2));mc];
        mc2(1,1)=1;
        % plot "heat map"
        surf(v2,[v1(1)-v1(2) v1],zeros(length(v1)+1,length(v2)),mc2);
        axis([min(u2) max(u2) min(u1) max(u1) 0 1])
        set(gca,'fontsize',14)
        xlabel(xlab,'fontsize',14)
        ylabel(ylab,'fontsize',14)
        % view from above
        view([0 90])
        axis square
        colorbar
        title(s)
    end
    
elseif nRepit==1
    figure
    hold on
    plot([0:maxGen],x')
    axis([0 maxGen 0 1])
    set(gca,'fontsize',14)
    
%     figure
%     scatter(x(4,2:gen),(w(4,2:gen)-1)/selCoef,'k.')
    
end

toc

close(wb)



