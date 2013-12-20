% the information content of the tallied individuals
% higher tn causes higher information content

clear all

tic

wb=waitbar(0,'progress');

rand('seed',17413);
ttmax=50000;
regime=1;
pincr=1;
[pA,pB] = randomenvironment2(ttmax,regime,0.02,pincr,0.5,0.5); 
pab=pA>=pB; 
mat=.5*ones(50,ttmax); 
for tau=1:50 
    for tt=1:ttmax-tau 
        mat(tau,tt)=(pab(tt)==pab(tt+tau)); 
    end
end 
reliability=mean(mat,2);

nindi=100;        % # of individuals
yymax=7;          % how often the whole procedure is iterated till convergence
tn=7;             % tally number
tmax=50;          % max. age of information
tt=[1:tmax];      % time vector
c=.5738            % correspondence
val=reliability'; % value of information
val=val*c+(1-val)*(1-c);    % probability that information is accurate
                            % either because correct + no switch
                            % or false + n switches (n odd);
if mod(tn,2)==0
    'error: only works for odd tally numbers'
    break
end
                            
                            
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
% the last x becomes the residual
x(tmax,:)=M-sum(x);
% the last value of information is just the guess, 50%
val(length(val)+1)=.5;
val=repmat(val,nindi+1,1);

% normalize x
xx=[N;x];
xx=xx/nindi;



%probability of success
prob=zeros(nindi+1,tmax+1);
probyy=zeros(nindi+1,yymax);

yy=1;



for t=1:tmax+1

    for i=1:nindi+1

        % j is the number of maj. talliers among the tallied group
        % for j=0:
        prob(i,t)=prob(i,t)+xx(t,i)*...
            binopdf2(tn,tn,(nindi-i+1)/nindi)*...
            sum(binopdf2([ceil(tn/2):tn],tn,c));
        % for j=tn:
        prob(i,t)=prob(i,t)+xx(t,i)*...
            binopdf2(0,tn,(nindi-i+1)/nindi)*...
            sum(binopdf2([ceil(tn/2):tn],tn,val(i,t)));
        for j=1:tn-1        % # of majority talliers
            k=tn-j;         % # of individual learners
            for jj=0:j
                for kk=0:k
                    if jj+kk>=ceil(tn/2)
                        prob(i,t)=prob(i,t)+xx(t,i)*...             %prob. that info. is t periods old
                            binopdf2(k,tn,(nindi-i+1)/nindi)*...    %prob. of having k ind. learners
                            binopdf2(kk,k,c)*...                    %prob. that of these k, kk choose correctly
                            binopdf2(jj,j,val(i,t));                  %prob. that of j maj. talliers, jj choose correctly
                    end
                end
            end        
        end
    end
end

'achtung'
% prob=ones(nindi+1,tmax+1)/(tmax+2);
prob=prob/10;

probyy(:,yy)=sum(prob');
waitbar(yy/yymax,wb);

% iteration of the run

for yy=2:yymax

    pp=sum(prob');

    reliability(tmax+1)=.5;
    val=repmat(pp,tmax+1,1)';
    % probability that information is accurate
                                % either because correct + no switch
                                % or false + n switches (n odd);
    % the last value of information is just the guess, 50%
    for t=1:tmax+1
        val(:,t)=val(:,t)*reliability(t)+(1-val(:,t))*(1-reliability(t));
    end

    if mod(tn,2)==0
        'error: only works for odd tally numbers'
        break
    end

    %probability of success
    prob=zeros(nindi+1,tmax+1);
    for t=1:tmax+1

        for i=1:nindi+1

            % j is the number of maj. talliers among the tallied group
            % for j=0:
            prob(i,t)=prob(i,t)+xx(t,i)*...
                binopdf2(tn,tn,(nindi-i+1)/nindi)*...
                sum(binopdf2([ceil(tn/2):tn],tn,c));
            % for j=tn:
            prob(i,t)=prob(i,t)+xx(t,i)*...
                binopdf2(0,tn,(nindi-i+1)/nindi)*...
                sum(binopdf2([ceil(tn/2):tn],tn,val(i,t)));
            for j=1:tn-1        % # of majority talliers
                k=tn-j;         % # of individual learners
                for jj=0:j
                    for kk=0:k
                        if jj+kk>=ceil(tn/2)
                            prob(i,t)=prob(i,t)+xx(t,i)*...             %prob. that info. is t periods old
                                binopdf2(k,tn,(nindi-i+1)/nindi)*...    %prob. of having k ind. learners
                                binopdf2(kk,k,c)*...                    %prob. that of these k, kk choose correctly
                                binopdf2(jj,j,val(i,t));                  %prob. that of j maj. talliers, jj choose correctly
                        end
                    end
                end        
            end
        end
    end
    probyy(:,yy)=sum(prob');
    waitbar(yy/yymax,wb);
end

close(wb);

figure
plot(M/nindi,probyy)
hold on
plot([0 1],[c c],'k--')
axis([0 1 0 1]);
title(tn)
xlabel('frequency of majority tallier')
ylabel('mean correct choices (%)')

toc