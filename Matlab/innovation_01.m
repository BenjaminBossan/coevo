% how the rate of innovation and social learning contribute to the spread
% of innovations
% open questions: can you only advance one innovation at a time?
%                 can you only socially learn from those exactly one step
%                   ahead?

% PARAMETER

tmax=1000;       % max generation time
nindi=1000;      % population size
lr=0.0;         % loss rate of learned innovations
nparam=1;       % # of parameters of population tensor
innoCap=7;     % maximal level of innovation
nIter=1;       % # of iterations
learnMode=2;    % how to learn socially: 1->direct bias, 2->conformism

% indices
iinno=1;

% INITIALIZE

% irVec=[1/10000:1/10000:2/1000];
% slVec=[1/10:2/100:2];
% slVec=[1:2:21];

irVec=1/100000;
slVec=3;

tMat=zeros(length(irVec),length(slVec));
tMatVar=tMat;

tic

nRepit=nIter*length(irVec)*length(slVec);

if nRepit>100
    wb=waitbar(0,'progress');
end


for ii=1:length(irVec)
    for ss=1:length(slVec)
        [t1 t2 innoVec]=coevo_innovation_01(...
            tmax,nindi,irVec(ii),slVec(ss),...
            lr,nparam,innoCap,nIter,learnMode,nRepit);
        tMat(ii,ss)=t1;
        tMatVar(ii,ss)=t2;
    end
    if nRepit>100
        waitbar(ii/length(irVec),wb);
    end
end

if nRepit>100
    close(wb)
end


toc

if nRepit>8

    figure

    plot3(repmat(irVec,length(slVec),1)',repmat(slVec,length(irVec),1),log(tMat),'k')
    hold on
    plot3(repmat(irVec,length(slVec),1),repmat(slVec,length(irVec),1)',log(tMat)','k')

    xlabel('innovation rate')
    ylabel('social learning parameter')
    zlabel('log time needed')
    set(gca,'fontsize',14)

    figure

    surf(repmat(irVec,length(slVec),1)',repmat(slVec,length(irVec),1),log(tMat))
    colormap gray

    xlabel('innovation rate')
    ylabel('social learning parameter')
    zlabel('log time needed')
    set(gca,'fontsize',14)
    
elseif nRepit==1
    
    figure
    plot(innoVec')
    
end


