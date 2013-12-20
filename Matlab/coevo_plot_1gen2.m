% plot function in case of just one generation
% for 1) pure strategies, shows how each strategy has chosen on average in
% each period
% for 2) mixed strategies, shows how all indidivuals have chosen

function fig1=coevo_plot_1gen2(n,nt,nindi,tmax,nstrat,genetics,pA,pB,...
    istrat,ichoice,itrack,irecent,ifit)

% in case of PURE STRATEGIES, find amount of pure strategies
if genetics==0
    
    counter=0;    
    %find unique strategies
    [unik,unikchar]=coevo_find_unique(n,istrat,nstrat);
    
    
    %# of A choices of individuals per strategy over periods
    fig1=zeros(length(unik),tmax);
    for t=1:tmax
        for j=1:length(unik)
            ind0=find((n(:,ichoice,t)==1)&(n(:,istrat,t)==unik(j)));
            ind1=find(n(:,istrat,t)==unik(j));
%             fig1(j,t)=sum(n(ind0,ichoice,t))/max(1,sum(n(:,istrat,t)))*unik(j);
            fig1(j,t)=length(ind0)/length(ind1);
        end
    end


    % plot the actual figure
    % in case of pure strategies: each strategy seperately
    figure
    if sum(find(fig1(1,:)~=0))>0
        plot(fig1(1,:)','r.')
    end
    hold on
    if length(unik)>1
        if sum(find(fig1(2,:)~=0))>0
            plot(fig1(2,:)','b.')
        end
        if length(unik)>2
            if sum(find(fig1(3,:)~=0))>0
                plot(fig1(3,:)','m.')
            end
            if length(unik)>3
                if sum(find(fig1(4,:)~=0))>0
                    plot(fig1(4,:)','k.')
                end
                if length(unik)>4
                    if sum(find(fig1(5,:)~=0))>0
                        plot(fig1(5,:)','g.')
                    end
                    if length(unik)>5
                        if sum(find(fig1(6,:)~=0))>0
                            plot(fig1(6,:)','.c')
                        end
                        if length(unik)>6
                            if sum(find(fig1(7,:)~=0))>0
                                plot(fig1(7,:)','rx')
                            end
                            if length(unik)>7
                                if sum(find(fig1(8,:)~=0))>0
                                    plot(fig1(8,:)','bx')
                                end
                                if length(unik)>8
                                    if sum(find(fig1(9,:)~=0))>0
                                        plot(fig1(9,:)','mx')
                                    end
                                    if length(unik)>9
                                        if sum(find(fig1(10,:)~=0))>0
                                            plot(fig1(10,:)','kx')
                                        end
                                        if length(unik)>10
                                            if sum(find(fig1(11,:)~=0))>0
                                                plot(fig1(11,:)','gx')
                                            end
                                            if length(unik)>11
                                                if sum(find(fig1(12,:)~=0))>0
                                                    plot(fig1(12,:)','cx')
                                                end
                                                if length(unik)>12
                                                    if sum(find(fig1(13,:)~=0))>0
                                                        plot(fig1(13,:)','ro')
                                                    end
                                                    if length(unik)>13
                                                        if sum(find(fig1(14,:)~=0))>0
                                                            plot(fig1(14,:)','bo')
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
%     grid on

%     % plot the choice of the 10% fittest individuals
%     figfit=zeros(1,tmax);
%     tmin=1;
%     for t=1:tmax
%         ind0=find((n(:,ichoice,t)==1)&(n(:,ifit,t)>quantile(n(:,ifit,t),.9)));
%         ind1=find(n(:,ifit,t)>quantile(n(:,ifit,t),.9));
%         if length(ind1)>0
%             figfit(t)=length(ind0)/length(ind1);
%         else
%             tmin=tmin+1;
%             figfit(t)=1/2;
%         end
%     end
%     plot([tmin:tmax],figfit(tmin:tmax),'ro')
%     unikchar(length(unik)+1)={'fittest 10%'};
%     unikchar=unikchar';
    
    legend(unikchar)
    
    plot(.5*ones(1,tmax),'k-')      % 50% line
    % plot of the best hunting ground choice (A->1, B->0)
    pmarkov=(1-pB)./(2-pA-pB);      %steady state vector according to markov chain analysis
%     plot(pmarkov,'k')
    pbetterA=pA.*(1-pB)./(pA+pB-2*pA.*pB);  % percent time A is better
    plot(.5+pA-pB,'b')
    axis([1 tmax 0 1])
    title('percent of A choices according to strategy')   
    
    hold off
%     % to find the # of times that the best choice switched:
%     pvec2=(pbetterA-.5)>=0;
%     pvec3=pvec2(1:tmax-1)-pvec2(2:tmax);
%     switchtimes=sum(abs(pvec3));

    if find(unik==4)~=0     %if there are wealth talliers
        
        % create a figure that tracks the SOURCE of information
        fig_wealth_info=zeros(length(unik)+1,tmax);
        % count for each period how often a certain strategy was imitated
        for t=1:tmax
            for s=1:length(unik)
                fig_wealth_info(s,t)=sum(find(n(:,itrack,t)==unik(s))~=0);
            end
            %initial guesses
            fig_wealth_info(length(unik)+1,t)=sum(intersect...
                (find(n(:,itrack,t)==9),find(n(:,istrat,t)==4))~=0);
        end
        % remove wealth imitators as possible sources
        fig_wealth_info(find(unik==4),:)=[];        
        % adjust legend
        unikchar2=unikchar;
        unikchar2(find(unik==4))=[];
        unikchar2(find(unik==length(unik)))=cellstr('initial guess');
        figure
        plot((fig_wealth_info/max(max(fig_wealth_info)))')
        title('which strategies were imitated by wealth imitators')
        legend(unikchar2)
        
        % create a figure that tracks the AGE of information
        fig_wealth_age=zeros(2,tmax);
        fig_wealth_age(1,:)=squeeze(mean(n(:,irecent,:),1));
        fig_wealth_age(2,:)=squeeze(median(n(:,irecent,:),1));
        figure
        plot(fig_wealth_age(1,:)')
        title('mean age of information')
%         legend({['mean','median']})
        
    end
    
end

% in case of MIXED STRATEGIES
if genetics==1
    
    fig1=zeros(1,tmax);
    for t=1:tmax
        ind0=find((n(:,ichoice,t)==1));
        fig1(t)=length(ind0)/nindi;
    end
    figure
    plot(fig1(1,:)','r.')
    hold on
    plot(.5*ones(1,tmax),'k-')      % 50% line
    % plot of the best hunting ground choice (A->1, B->0)
    pmarkov=(1-pB)./(2-pA-pB);      %steady state vector according to markov chain analysis
    plot(.5+pA-pB,'k')
    pbetterA=pA.*(1-pB)./(pA+pB-2*pA.*pB);  % percent time A is better
%     plot(pbetterA,'b')
    axis([1 tmax 0 1])
    title('percent of A choices according to strategy')
    
%     % plot the choice of the 10% fittest individuals
%     figfit=zeros(1,tmax);
%     tmin=1;
%     for t=1:tmax
%         ind0=find((n(:,ichoice,t)==1)&(n(:,ifit,t)>quantile(n(:,ifit,t),.9)));
%         ind1=find(n(:,ifit,t)>quantile(n(:,ifit,t),.9));
%         if length(ind1)>0
%             figfit(t)=length(ind0)/length(ind1);
%         else
%             tmin=tmin+1;
%             figfit(t)=1/2;
%         end
%     end
%     plot([tmin:tmax],figfit(tmin:tmax),'ro')
    
%     % plot the choice of the individuals that did not play ITW
%     fignotitw=zeros(1,tmax);
%     tmin=1;
%     for t=1:tmax
%         ind0=find((n(:,ichoice,t)==1)&(n(:,istrat,t)~=4));
%         ind1=find(n(:,istrat,t)~=4);
%         fignotitw(t)=length(ind0)/length(ind1);
%     end
%     plot([tmin:tmax],fignotitw(tmin:tmax),'rx')
    
    hold off
end

% % how often one option became better than the other
% pvec2=(pbetterA-.5)>=0;
% pvec3=pvec2(1:tmax-1)-pvec2(2:tmax);
% switch_probability_per_round=sum(abs(pvec3))/tmax;
% 
% % in how many periods one option was more than twice as good as the other
% pAtwiceB=sum(find(pA./pB>2)~=0);
% pBtwiceA=sum(find(pA./pB<1/2)~=0);
% percent_periods_X_better_than_2Y=(pAtwiceB+pBtwiceA)/tmax;



