% check whether there are only pure strategies
% if yes -> genetics = 0
% if not -> genetics = 1
% as soon as mutations are possible, genetics are not pure

function genetics=coevo_check_genetics2(mutationrate,nstrat,ninitial,ipind,nindi,gen)

% if there are mutations, there are no pure strategies
if gen<2
    genetics=0;
else
    if mutationrate==0
        genetics=0;
    else
        genetics=1;
    end
end

counter=0;

while genetics==0
    
    % stop if each individual was checked
    counter=counter+1;
    if counter>nindi
        break
    end
    
    % if it is not true for each strategy that it is chosen
    % with probability 0 or 1, then not all strategies are pure
    
    for s=1:nstrat
        if mod(ninitial(counter,ipind-1+s),1)~=0
            genetics=1;
        end
    end
    
end