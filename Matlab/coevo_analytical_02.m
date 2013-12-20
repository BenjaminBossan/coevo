% this is an attempt to make the model analytical. It is assumed that the
% performance of individual learners is constant. Path dependence is
% ignored. Social learners have 3 individuals as sample. Choice proportions
% are always in equilibrium.

clear all

% PARAMETERS
cind=0.6;                   % performance individual learners
tmax=50;                    % periods per generation
gen=100000;                 % generations
selCoef=1/10;               % maximum selection coefficient
wind=1;                     % fitness individual learners
probEnvChange=1/10;         % probability of environmental change
pC=.54;pF=.46;                % probabilities that correct/false options yield success
nStrats=7;                  % total possible number of strategies
windowSize=2500;            % window size for moving average filter
stabCrit=selCoef/100;       % stability criterion to estimate equilibria
genBonus=1000;              % # of generations that simulation continues to
                            % run after stability criterion is reached

tic

'cind adjusted to win-stay lose-shift'

% STRATEGIES
% 1 ind: individual
% 2 con: conformists
% 3 pbsl1: payoff-biased social learners equal weight
% 4 pbsl3: payoff-biased social learners, more weight on gains
% 5 mce: McElreath payoff-biased social learners
% 6 OIL
% 7 OC

xInit(1,1)=.9;
xInit(2,1)=1/60;
xInit(3,1)=1/60;
xInit(4,1)=1/60;
xInit(5,1)=1/60;
xInit(6,1)=1/60;
xInit(7,1)=1/60;

% xInit=1/nStrats;

% TESTED PARAMETER RANGE

% cindVec=[0.51:.01:.99];
% cindVec=.65;

% probEnvChangeVec=[.01:.02:.99];
% probEnvChangeVec=.1;

% pCvec=[0:.2:1];
% pFvec=[0:.2:1];
pCvec=.54;pFvec=0.46;

% tmaxVec=[10:20:100];
% average # of rounds before env. change occurs
nRoundVec=[2:1:50];

dpVec=[.1:.02:1];

% define the FOCAL PARAMETER VECTORS
% u1=cindVec;
% u2=probEnvChangeVec;
u2=dpVec;
u1=nRoundVec;
% u2=pCvec;
% u1=pFvec;
xlab='p_A - p_B';
ylab='p_{env}';

nRepit=length(u1)*length(u2);

m=zeros(nStrats,length(u1)*length(u2));
m=reshape(m,nStrats,length(u1),length(u2));
mv=m;   %variance

wb=waitbar(0,'progress');

for j=1:length(u1)
    
    probEnvChange=1/nRoundVec(j);
    % pC=pCvec(j);

    for i=1:length(u2)

        pC=0.5+dpVec(i)/2;
        pF=1-pC;
        % pF=pFvec(i);
        
        cind=(1-pF)/(1-pC+1-pF);
        
        if pC>pF

            [x c w maxGen]=coevo_ana_03(...
            xInit,...                   % initial conditions
            cind,...                    % performance individual learners
            tmax,...                    % periods per generation
            gen,...                     % generations
            selCoef,...                 % maximum selection coefficient
            wind,...                    % fitness individual learners
            probEnvChange,...           % probability of environmental change
            pC,pF,...                   % probabilities that correct/false options yield success
            nStrats,...                 % total possible number of strategies
            windowSize,...              % window size for moving average filter
            stabCrit,...                % stability criterion to estimate equilibria
            genBonus);                  % # of generations that simulation continues to
                                        % run after stability criterion is reached

            xLast=mean(x(:,maxGen-999:maxGen),2);
            xVar=var(x(:,maxGen-999:maxGen)')';
            m(:,j,i)=xLast;
            mv(:,j,i)=xVar;
            
        end
    end
    
    waitbar(j/length(u1),wb);
    
end

% FIGURES

if nRepit>=6

%     figure
% 
%     plot3(repmat(cindVec,length(probEnvChangeVec),1)',repmat(probEnvChangeVec,length(cindVec),1),mc,'k')
%     hold on
%     plot3(repmat(cindVec,length(probEnvChangeVec),1),repmat(probEnvChangeVec,length(cindVec),1)',mc','k')
%     set(gca,'fontsize',14)
%     xlabel('performance individual learners','fontsize',14)
%     ylabel('probability of environmental change','fontsize',14)
%     zlabel('frequency of conformists','fontsize',14)
% 
%     figure
% 
%     plot3(repmat(cindVec,length(probEnvChangeVec),1),repmat(probEnvChangeVec,length(cindVec),1)',mc')
%     colormap gray
%     set(gca,'fontsize',14)
%     xlabel('performance individual learners','fontsize',14)
%     ylabel('probability of environmental change','fontsize',14)
%     zlabel('frequency of conformists','fontsize',14)
    
    for s=1:nStrats
        mc=squeeze(m(s,:,:));
        figure
        colormap([1-[0:.01:1]' 1-[0:.01:1]' 1-[0:.01:1]']);
        mc2=[mc ;zeros(1,length(u2))];
        mc2(length(u1)+1,length(u2))=1;
        lastEntry=2*u1(length(u1))-u1(length(u1)-1);
        surf(u2,[u1 lastEntry],zeros(length(u1)+1,length(u2)),mc2)
        axis([min(u2) max(u2) min(u1) max(u1) 0 1])
        set(gca,'fontsize',14)
        ylabel(ylab,'fontsize',14)
        xlabel(xlab,'fontsize',14)
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
end

close(wb)

toc