% the information content of the tallied individuals
% higher tn causes higher information content

tic

nindi=100;       % # of individuals
tn=3;             % tally number
tmax=50;          % max. age of information
tt=[1:tmax];      % time vector
c=0.7;            % correspondence

rand('seed',17413);
ttmax=25000;regime=1;pincr=1;
[pA,pB] = randomenvironment2(ttmax,regime,0.02,pincr,0.5,0.5); 
pab=pA>=pB; 
mat=.5*ones(tmax,ttmax); 
for tau=1:tmax 
    for tt=1:ttmax-tau 
        mat(tau,tt)=(pab(tt)==pab(tt+tau)); 
    end 
end
reliability=mean(mat,2);

val=reliability'; % value of information
val=val*c+(1-val)*(1-c);    % probability that information is accurate
                            % either because correct + no switch
                            % or false + n switches (n odd);

M=0:nindi;        % # of majority talliers
N=nindi-M;        % # of individual learners


xmax=tmax+1;
x=zeros(xmax,nindi+1);
for iter=1:xmax
    if iter==1
        x(iter,:)=N;
    else
        x(iter,:)=x(iter-1,:)+N/nindi.*(nindi-x(iter-1,:));
    end
end
x=diff(x);
x(tmax,:)=M-sum(x);
% the last x is just the residual
val(length(val)+1)=.5;
% the last value is just the guess

xx=[N;x];
xx=xx/nindi;
v1=[0:tn];
v2=1-v1;
prob=ones(nindi+1,1);
% for t=1:tmax+1
%     for i=1:nindi+1
%         % j is the number of maj. talliers among the tallied group
%         % for j=0:
%         prob(i)=prob(i)-xx(t,i)*binopdf2(tn,tn,(nindi-i+1)/nindi)*sum(binopdf2([0:floor(tn/2)],tn,c));
%         % for j=tn:
%         prob(i)=prob(i)-xx(t,i)*binopdf2(0,tn,(nindi-i+1)/nindi)*sum(binopdf2([0:floor(tn/2)],tn,val(t)));
%         for j=1:tn-1;
%             for k=1:min(tn+1-v1(j+1),tn+1-v2(j+1));
%                 prob(i)=prob(i)-...
%                     xx(t,i)*...                             %probability that information is this old
%                     binopdf2(tn-k,tn,(nindi-i+1)/nindi)*...    %probability of that many ind
%                     binopdf2(k,tn,(nindi-i+1)/nindi)*...       %probability of that many maj
%                     sum(binopdf2([0:ceil(tn/2)-k],tn,c))*... %through individual learners
%                     sum(binopdf2([0:k],tn,val(t)));        %through majority talliers
%             end
%         end
%     end
% end

for t=1:tmax+1
    for i=1:nindi+1
        % j is the number of maj. talliers among the tallied group
        % for j=0:
        prob(i)=prob(i)-xx(t,i)*binopdf2(tn,tn,(nindi-i+1)/nindi)*sum(binopdf2([0:floor(tn/2)],tn,c));
        % for j=tn:
        prob(i)=prob(i)-xx(t,i)*binopdf2(0,tn,(nindi-i+1)/nindi)*sum(binopdf2([0:floor(tn/2)],tn,val(t)));
        for j=1:tn-1;
            for k=1:min(tn+1-v1(j+1),tn+1-v2(j+1));
                prob(i)=prob(i)-...
                    xx(t,i)*...                             %probability that information is this old
                    binopdf2(tn-k,tn,(nindi-i+1)/nindi)*...    %probability of that many ind      %probability of that many maj
                    sum(binopdf2([0:ceil(tn/2)-k],tn-j,c))*... %through individual learners
                    sum(binopdf2([0:k],j,val(t)));        %through majority talliers
            end
        end
    end
end

figure
plot(M/nindi,prob,'k')

toc

