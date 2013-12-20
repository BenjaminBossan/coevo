function [fig_wealth_info,fig_wealth_age,unikchar2] = ...
    coevo_wealth_info_age(unik,tmax,n,itrack,istrat,unikchar,irecent)

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

% fig_wealth_age=zeros(2,tmax);
% fig_wealth_age(1,:)=squeeze(mean(n(:,irecent,:),1));
% fig_wealth_age(2,:)=squeeze(median(n(:,irecent,:),1));

fig_wealth_age=zeros(2,tmax);
% find wealth imitators and measure their mean information age
for t=1:tmax
    fig_wealth_age(1,t)=squeeze(mean(n(n(:,istrat,t)==4,irecent,t),1));
end




% tip: to find the number of periods in which the actual age is significantly
% greater than what would be expected from random sampling, use:

% xx=3/fraction_wealth;tt=[1:tmax];randage=(1-xx.^(tt-1))/(1-xx);
% pvec=ones(1,tmax);
% for i=3:tmax
% pvec(i)=mybinomtest(sum(find(wealth_age(i,:)>randage(i))~=0),iterations,1/2);
% end
% hold on
% plot(mean(wealth_age,2))
% plot(randage,'r--')
% in_how_many_perids_significant_difference=sum(find(pvec<.05/tmax)~=0)
% higher=0;lower=0;
% ind0=find(pvec<.05/tmax);
% for j=1:length(ind0)
% higher=higher+...
% ((sum(wealth_age(ind0(j),:)>randage(ind0(j))))>...
% (sum(wealth_age(ind0(j),:)<randage(ind0(j)))));
% lower=lower+...
% ((sum(wealth_age(ind0(j),:)>randage(ind0(j))))<...
% (sum(wealth_age(ind0(j),:)<randage(ind0(j)))));
% end
% higher
% lower

% randage would be the age expected from random sampling
% we need to find p-values less than 0.05/tmax because we apply
% bonferroni correction (p_use = p_desired / repetitions)
% higher and lower state in how many cases the actual age was
% higher or lower than the age from random sampling