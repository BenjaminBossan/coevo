function fig1=coevo_plot_gens_4(n,nt,nindi,tmax,nstrat,genetics,gen,...
    istrat,ichoice,ipind,ntsize)

if genetics==0
    
    counter=0;
    
    [unik,unikchar]=coevo_find_unique2(nt,istrat,nstrat);
    
    %just strategies
    nts=squeeze(nt(:,istrat,:));
    
    fig1=zeros(length(unik),gen);
    
    counter=0;
    for s=1:length(unik)
        if sum(find(nts(:,s,1)==unik(s)))>0
            counter=counter+1;
            for g=1:gen                
                fig1(counter,g)=sum(nts(:,g)==unik(s));
            end
        end
    end
    
%     % running average over 20 generations
%     windowSize=10;
%     meancorrflat=filter(ones(1,windowSize)/windowSize,1,meancorr);
    
    figure
    hold on
    plot([0:gen-1],fig1'/nindi)
    
    xlabel('generation')
    ylabel('frequency of the strategies')
    
%     if gen>50
%         plot([windowSize+1:gen],meancorrflat(windowSize+1:gen),'m--')
%     else
%         plot(meancorr,'m--')
%     end
    
    legend(unikchar)
    axis([0 gen 0 1])
    hold off
end
% 
% if genetics==1
%     
%     % the probabilities of the strategies
%     nts=squeeze(nt(:,ipind:ipind+nstrat-1,:));
%     
%     % find unique genes for strategies
%     counter=0;
%     for s=1:nstrat
%         if sum(find(nts(:,s,:)~=0)~=0)>1
%             counter=counter+1;
%             unik(counter)=s;
%             if s==1
%                 unikchar(counter)=cellstr('individual');
%             elseif s==2
%                 unikchar(counter)=cellstr('majority');
%             elseif s==3
%                 unikchar(counter)=cellstr('success');
%             elseif s==4
%                 unikchar(counter)=cellstr('wealth');
%             elseif s==5
%                 unikchar(counter)=cellstr('ind or maj');
%             elseif s==6
%                 unikchar(counter)=cellstr('maj or ind');
%             end
%         end
%     end
%     
%     fig1=zeros(length(unik),gen);
%     
%     counter=0;
%     for s=1:length(unik)
%         for g=1:gen                
%             fig1(s,g)=mean(nts(:,unik(s),g));
%         end
%     end
%     
%     meancorr=fig5;
%     % running average over 20 generations
%     windowSize=20;
%     meancorrflat=filter(ones(1,windowSize)/windowSize,1,meancorr);
%     
%     
% %     figure
%     hold on
%     plot(fig1')
%     
%     title('talliers over time,w0=10,tmax=50,tally#=7,incr=.02,pincr=1,c=0,b=1')
%     xlabel('generations')
%     ylabel('solid:individuals,dashed:mean(running average) fitness')
%     
%     if gen>50
%         plot([10:g-10],meancorrflat(10:g-10),'m--')
%     else
%         plot(meancorr,'m--')
%     end
%     
%     legend(unikchar)
%     
%     axis([0 gen 0 1])
%     hold off
%     
% end








