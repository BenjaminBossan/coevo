% how many correct choices were made this round
% this function only monitors how the group as a whole
% performed and does not distinguish between strategies.
% the function is only used to show the average performance
% over generations.

function fig4=coevo_meancorrect3(pA,pB,n,nindi,ichoice,tmax)

pAbest=repmat(sign(pA-pB),nindi,1);         %=1 if A better, else -1
ntc=squeeze(n(:,ichoice,:));               %choices
ntc=2*ntc-1;
equalAB=sum(find(pAbest(1,:)==0)~=0);       %rounds in which neither was better
fig2=min(max(0,pAbest.*ntc),1);
fig3=sum(fig2,2)/(tmax-equalAB);            %percent rounds with correct choice
fig4=mean(fig3);
