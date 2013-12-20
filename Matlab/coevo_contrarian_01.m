% how many of the contrarian views in the population are expressed
% by social learners and how many by individual learners
% input: superfig
% fsl: frequency of social learners

function coevo_contrarian_01(superfig,fsl);

% freq. of individual learners
fil=1-fsl;

% mean behavior of the population
sfmean=fsl*superfig(2,:)+fil*superfig(1,:);

% from this, derive the contrarian option
sfcontra=1-sfmean>1/2;

% percent of contrarian views expressed by social learners:
slcontra=superfig(2,:).*sfcontra+(1-superfig(2,:)).*(1-sfcontra);

% and by individual learners
ilcontra=superfig(1,:).*sfcontra+(1-superfig(1,:)).*(1-sfcontra);

absolute_prop_of_contrarian_views_sl=mean(slcontra)

absolute_prop_of_contrarian_views_il=mean(ilcontra)

perc_contrarian_views_by_il_of_all_contrarian_views=fil*mean(ilcontra)/...
    (fil*mean(ilcontra)+fsl*mean(slcontra))*100

leverage_of_il_compared_to_their_frequency=perc_contrarian_views_by_il_of_all_contrarian_views/fil/100