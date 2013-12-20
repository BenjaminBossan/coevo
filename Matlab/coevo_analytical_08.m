% this is an attempt to make the model analytical. It is assumed that the
% performance of individual learners is constant. Path dependence is
% ignored. Social learners have 3 individuals as sample. Choice proportions
% are always in equilibrium.

clear all

% PARAMETERS
tmax=50;                    % periods per generation
gen=100000                % generations
selCoef=0.5               % maximum selection coefficient
nStrats=19;                 % total possible number of strategies
vertTransm=1               % is there vertical tranmission from parents to offspring?
                            % 1->yes, 0->no, 0.5->yes but not for
                            % conformists
regime=1;                   % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;                 % increment at which the environment becomes better or worse
pincr=1;                    % probability that environmental quality changes at all after each period
dpA=0;dpB=0;                % shift of pA and pB
indLearn=0;                 % if == 0, ind. learning is win-stay lose-shift, else its constant
mutationRate=1/100;         % probability of mutation occurring

tic

% STRATEGIES
%  1 individual learners
%  2 conformists tn 3
%  3 conformists tn 7
%  4 OILs tn 3
%  5 OILs tn 7
%  6 OCs tn 3
%  7 OCs tn 7
%  8 pbsl equal weight tn 3
%  9 pbsl equal weight tn 7
% 10 pbsl gains/losses 3/1 tn 3
% 11 pbsl gains/losses 3/1 tn 7
% 12 pbsl gains/losses 1/3 tn 3
% 13 pbsl gains/losses 1/3 tn 7
% 14 pbsl conf gains/losses 3/1 tn 3
% 15 pbsl conf gains/losses 3/1 tn 7
% 16 pbsl only gains tn 3
% 17 pbsl only gains tn 7
% 18 pbsl McElreath tn 3
% 19 pbsl hybrid tn 3


xInit=zeros(nStrats,1);

% 'xInit'
% xInit(1)=9/10;
% xInit(2)=1/20;
% xInit(3)=1/20;

% 'no pbsls'
xInit(1)=6/7;
xInit(2)=1/7;
% xInit(3)=1/7;
% xInit(4)=1/7;
% xInit(5)=1/7;
% xInit(6)=1/7;
% xInit(7)=1/7;
% xInit(8)=0;
% xInit(9)=0;
% xInit(10)=0;
% xInit(11)=0;
% xInit(12)=0;
% xInit(13)=0;
% xInit(14)=0;
% xInit(15)=0;
% xInit(16)=0;
% xInit(17)=0;
% xInit(18)=0;
% xInit(19)=0;

% 'xInit'
% xInit(1)=.9;
% xInit(2:nStrats)=.1/(nStrats-1)*ones(1,nStrats-1);

% xInit=ones(nStrats,1)/nStrats;

% TESTED PARAMETER RANGE
incrVec=[.01:.01:.1];
incrVec=.02
pincrVec=[.1:.1:1];
pincrVec=1
% incrVec=.02
% pincrVec=1
% incrVec=.02;
% pincrVec=1;
% regimeVec=[0 .1 .25 .5 1 1.5 2 2.5 3];
% regimeVec=1
% regimeVec=1
% tmaxVec=[10 20 30 40 50 75 100 150 250 500];
% indVec=[.51:.02:.99];
% indVec=0
% tmaxVec=50
% tmaxVec=[5 50 100 250];
% tmaxVec=50;
% dpAvec=[-.25:.01:.25];
% dpBvec=[-.25:.05:.25];

% define the FOCAL PARAMETER VECTORS
% u2=dpBvec;
u2=pincrVec;
u1=incrVec;
% u1=dpAvec;
% u2=tmaxVec;

xlab='p_{incr}';
ylab='k_{incr}';

nRepit=length(u1)*length(u2);

m=zeros(nStrats,length(u1)*length(u2));
m=reshape(m,nStrats,length(u1),length(u2));
mv=m;   %variance

wb=waitbar(0,'progress');

for j=1:length(u1)
    
    incr=u1(j);
    
    for i=1:length(u2)
        
        pincr=u2(i);
        
        pA0=1/2;pB0=1/2;            % initial values of pA and pB, so that pA=/=pB forall t
            
        [x c w maxGen pA pB]=coevo_ana_08(...
        xInit,...                   % initial conditions
        tmax,...                    % periods per generation
        gen,...                     % generations
        selCoef,...                 % maximum selection coefficient
        nStrats,...                 % total possible number of strategies
        vertTransm,...              % whether there is vertical transmission
        indLearn,...                % ind. learning performance
        mutationRate,...            % prob. of mutation occurring
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
        % the focal strategy
        mc=squeeze(m(s,:,:));
        figure
        % colormap from white->low freq to black->high freq
        colormap([1-[0:.01:1]' 1-[0:.01:1]' 1-[0:.01:1]']);
        % add 0 and 1 to data so that full scale [0;1] is used
        mc2=[mc;zeros(1,length(u2))];
        mc2(length(u1)+1,length(u2))=1;
        lastEntry=2*u1(length(u1))-u1(length(u1)-1);
        % plot "heat map"
        surf(u2,[u1 lastEntry],zeros(length(u1)+1,length(u2)),mc2);
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