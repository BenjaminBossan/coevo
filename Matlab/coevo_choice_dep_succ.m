% plots choice of a strategy making visible how it depends
% on the history of success and losses
% illustrates how the individual learning strategies work

i=2;

ind1=find(n(i,ichoice,:)==1);
ind0=find(n(i,ichoice,:)==0);
vec1=squeeze(cumsum(n(i,isucc,ind1)-.5,3));
vec0=squeeze(cumsum(n(i,isucc,ind0)-.5,3));
vecsum=zeros(1,tmax);
c1=0;c0=0;
for t=1:tmax
    if n(i,ichoice,t)==1
        c1=c1+1;
        vecsum(t)=vec1(c1);
    else
        c0=c0+1;
        vecsum(t)=vec0(c0);
    end
end
figure
hold on
plot(ind1,vec1,'.')
plot(ind0,vec0,'r.')
legend('choice for A','choice for B')
plot(vecsum,'k')
plot(ind1,vec1,'.')
plot(ind0,vec0,'r.')
xlabel('period')
ylabel('success with A or B')