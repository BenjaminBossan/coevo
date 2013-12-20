function [tMean tVar innoVec] = coevo_innovation_01(...
    tmax,...
    nindi,...
    ir,...
    sl,...
    lr,...
    nparam,...
    innoCap,...
    nIter,...
    learnMode,...
    nRepit);


% indices
iinno=1;

tVec=zeros(1,nIter);

for rep=1:nIter

    % INITIALIZATION

    n=zeros(nindi,tmax*nparam);
    n=reshape(n,nindi,nparam,tmax);

%     innoVec=zeros(innoCap,t);

    % ALGORITHM

    for t=2:tmax

        % innovation
        n(:,iinno,t)=n(:,iinno,t-1)+(rand(nindi,1)<ir);

        % loss of innovation
        n(:,iinno,t)=n(:,iinno,t)-(rand(nindi,1)<lr);
        
        if learnMode==1
            % social learning, direct bias
            unik=[min(n(:,iinno,t)):max(n(:,iinno,t))];
            for j=min(unik):max(unik)-1
                indLearn=find(n(:,iinno,t)==j);
                indTeach=find(n(:,iinno,t)==j+1);
                n(indLearn,iinno,t)=...
                    n(indLearn,iinno,t)+...
                    (rand(length(indLearn),1)<(length(indTeach)/nindi).^(1/sl));
            end
            
        elseif learnMode==2
        
            % social learning, conformism
            unik=[min(n(:,iinno,t)):max(n(:,iinno,t))];
            for j=min(unik):max(unik)-1
                indLearn=find(n(:,iinno,t)==j);
                indTeach=find(n(:,iinno,t)==j+1);
                xVar=length(indTeach)/nindi;
                learnProb=sum(binopdf(ceil(sl/2),sl,xVar));
                n(indLearn,iinno,t)=...
                    n(indLearn,iinno,t)+...
                    (rand(length(indLearn),1)<learnProb);
            end
        end

        % caps on innovation
        n(:,iinno,t)=min(innoCap,max(0,n(:,iinno,t)));

        % record innovation status
        for j=1:innoCap
            innoVec(j,t)=sum(n(:,iinno,t)>=j);
        end

        if sum(n(:,iinno,t)==innoCap)>=nindi*.5
            break
        end
    end
    
    tVec(rep)=t;
    
end

tMean=mean(tVec);
tVar=var(tVec);

if nRepit==1
    innoVec=zeros(innoCap,t);
    for ii=1:innoCap
        for tt=1:t
            innoVec(ii,tt)=sum(n(:,iinno,tt)>=ii);
        end
    end
else
    innoVec=0;
end



