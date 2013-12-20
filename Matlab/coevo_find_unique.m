% find unique strategies
% input: n, istrat, nstrat
% output unik (list of unique strategy types),
%        unikchar (name of these strategies)

function [unik,unikchar]=coevo_find_unique(nt,istrat,nstrat)

counter=0;
%find unique strategies
%extract strategies
nts=squeeze(nt(:,istrat,:));
for s=1:nstrat
    % if strategy s is not never chosen
    if sum(find(nts(:,1)==s))~=0
        counter=counter+1;
        unik(counter)=s;
    end
end
unik=unik';

unikchar=cell(length(unik),1);
for l=1:length(unik)
    if unik(l)==1
        unikchar(l)=cellstr('individual');
    elseif unik(l)==2
        unikchar(l)=cellstr('conformist');
    elseif unik(l)==3
        unikchar(l)=cellstr('PBSL 1/1');
    elseif unik(l)==4
        unikchar(l)=cellstr('ITW');
    elseif unik(l)==5
        unikchar(l)=cellstr('OMT');
    elseif unik(l)==6
        unikchar(l)=cellstr('OIL');
    elseif unik(l)==7
        unikchar(l)=cellstr('reinforcement');
    elseif unik(l)==8
        unikchar(l)=cellstr('IDC');
    elseif unik(l)==9
        unikchar(l)=cellstr('omniscient');
    elseif unik(l)==10
        unikchar(l)=cellstr('PBSL 3/1');
    elseif unik(l)==11
        unikchar(l)=cellstr('PBSL 1/3');
    elseif unik(l)==12
        unikchar(l)=cellstr('PBSL only gains');
    elseif unik(l)==13
        unikchar(l)=cellstr('PBSL-conformist');
    elseif unik(l)==14
        unikchar(l)=cellstr('PBSL McElreath');
    end
end