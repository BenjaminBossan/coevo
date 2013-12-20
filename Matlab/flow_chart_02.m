% selection gradients diplayed in a flow chart, for 3 types

% freq of wealth, majority
% x=1/9*[1,2,3,1,4,1,2,4,5,7,3,6,0,0];
% y=1/9*[1,2,3,4,1,7,5,4,2,1,0,0,3,6];

% freq maj vs wea vs OMT a.o.
x=1/9*[1,2,3,1,4,1,2,4,5,7,3,6,0,0];
y=1/9*[1,2,3,4,1,7,5,4,2,1,0,0,3,6];

% freq of majority, success
% x=1/9*[1,2,3,1,4,1,2,4,5,7,0,0,3,6,3,6];
% y=1/9*[1,2,3,4,1,7,5,4,2,1,3,6,0,0,6,3];

% freq of wealth, majo, for pskill=1 or pskill=2
% x=1/9*[1 2 3 1 4 1 2 4 5 7 3 6 0 0];
% y=1/9*[1 2 3 4 1 7 5 4 2 1 0 0 3 6];


% coordinates, e.g. x=wealth, y=majo, 1-x-y=indi
z=0; %irrelevant

f=[x;y;1-x-y]';

% DATASETS

% test matrix
% M=...
%     [1 1 1;
%     1 1 1;
%     2 1 1;
%     1 1 1;
%     1 1 1;
%     1 1 1;
%     1 1 1;
%     1 1 1;
%     1 1 1;
%     1 1 1;
%     1 0 1;
%     1 0 1;
%     0 1 1;
%     0 1 1];

% fitness, or % correct choices: wealth-majority-individual
% M=...
%     [6675 6984 5965;...
%     6849 7313 5959;...
%     6884 7152 5948;...
% %     6000 0000 0000;... % for testing purposes
%     7013 7335 5960;...
%     6740 7215 5943;...
%     5373 4950 5929;...
%     6232 5908 5922;...
%     6006 5666 5946;...
%     6995 7275 5974;...
%     6827 7059 5954;...
%     
%     6605 0    5942;...
%     6678 0    5957;...
%     0    7228 5943;...
%     0    5895 5922];

% fitness, but not % correct choices: wealth-majority-individual
% M=...
%     [5200 5230 5129;
%     5181 5222 5095;
%     5232 5259 5130;
%     5236 5268 5126;
%     5176 5222 5084;
%     5080 5035 5125;
%     5168 5124 5115;
%     5136 5107 5117;
%     5201 5225 5107;
%     5177 5195 5081;
%     5175    0 5100;
%     5204    0 5128;
%     0    5204 5074;
%     0    5158 5142];

% fitness, or % correct choices: majo-success-indi
% M=...
%     [6367 5808 5627;...
%     6635 5823 5638;...
%     6923 5799 5614;...
%     6475 5797 5613;...
%     7019 5804 5623;...
%     6654 5834 5639;...
%     6798 5836 5632;...
%     7249 5818 5629;...
%     6915 5777 5621;...
%     3280 5798 5646;...
%     0    5803 5619;...
%     0    5809 5625;...
%     6688    0 5612;...
%     4762    0 5628;...
%     7135 5843    0;...
%     5159 5763    0];

% percent correct choices: wealth-majo-indi, pskill=1

% M=...
%     [5948 6265 5581;...
%     6104 6483 5580;...
%     6086 6219 5568;...
%     6148 6330 5588;...
%     5995 6326 5584;...
%     4249 2940 5570;...
%     5189 4301 5552;...
%     4798 3789 5606;...
%     6067 6174 5551;...
%     5805 5745 5574;...
%     5275 0    5312;...
%     5305 0    5327;...
%     0    6002 5295;...
%     0    3873 5307];

% percent correct choices: wealth-majo-indi, pskill=2

% M=...
%     [5291 5644 5287;...
%     5413 5804 5294;...
%     5546 5748 5335;...
%     5579 5806 5304;...
%     5323 5632 5282;...
%     4463 3233 5294;...
%     4724 3774 5327;...
%     4136 3050 5264;...
%     5438 5477 5290;...
%     5256 5020 5292;...
%     5271 0    5300;...
%     5317 0    5340;...
%     0    6000 5294;...
%     0    4020 5305];

% % maj vs wea vs OILs
% M=...
%     [6787 6635 6235;
%     6201 6230 5966;
%     5283 5530 5583;
%     6411 6419 6075;
%     5282 5521 5576;
%     5575 5802 5714;
%     5255 5511 5560;
%     5178 5326 5554;
%     4930 5129 5417;
%     5182 5269 5552;
%     6355    0 6037;
%     5003    0 5459;
%     0    6617 6242;
%     0    6338 6026];

% maj vs wea vs OMTs
% M=...
%     [6375 6482 6300;
%     5443 5732 5672;
%     5062 5402 5433;
%     6084 6248 6095;
%     5005 5336 5386;
%     5431 5682 5680;
%     5761 5999 5917;
%     4958 5119 5349;
%     5090 5285 5439;
%     4934 5035 5346;
%     5429    0 5674;
%     4922    0 5337;
%     0    6536 6364;
%     0    6369 6230];

% wea vs maj vs ind pskill=1
% M=...
%     [6054 6459 5676;
%     6284 6752 5721;
%     5972 6066 5653;
%     6218 6475 5683;
%     6149 6584 5701;
%     5212 5083 5656;
%     5244 4967 5609;
%     5297 5112 5651;
%     5960 5950 5691;
%     5421 5177 5623;
%     5927    0 5650;
%     6015    0 5682;
%     0    6786 5697;
%     0    5448 5698];


% % wea vs maj vs ind pskill=2 redone, corrected
% M=...
%     [5512 5920 5408;
%     5507 5893 5370;
%     5501 5470 5463;
%     5407 5416 5402;
%     5485 5818 5385;
%     4927 4902 5295;
%     4933 4860 5205;
%     5132 5099 5358;
%     4878 4720 5277;
%     4991 4870 5300;
%     5323    0 5341;
%     5355    0 5348;
%     0    6253 5448;
%     0    5040 5348];

% wea vs maj vs ind, tmax=250
M=...
    [6511 7089 6054;
    6721 7342 6064;
    6823 7169 6062;
    6926 7343 6056;
    6611 7256 6059;
    5344 4926 6056;
    6186 5981 6078;
    5756 5451 6038;
    6691 7081 6052;
    6549 6829 6050;
    6430    0 6069;
    6412    0 6062;
    0    7381 6054;
    0    6150 6063];
    
M=single(M);
% for some odd reason, matlab works with a different precision
% for the first row and all other rows. For instance, if M=ones(3,length),
% only the first point, derived from the first row, will point in any
% direction whereas the other points are neutral, as they should be.
% by converting M, this problem is solved.


%~~~~~~~~~~~~~~~~END DATA~~~~~~~~~~~~~~~~~~~~
    

meanM=zeros(length(M),1);
for j=1:length(M)
    meanM(j)=sum(M(j,:).*f(j,:));
end

% difference (equation) matrix
dM=zeros(length(M),3);
for j=1:length(M)
    dM(j,1)=         x(j)*(M(j,1)-meanM(j))/meanM(j);
    dM(j,2)=         y(j)*(M(j,2)-meanM(j))/meanM(j);
    dM(j,3)=(1-x(j)-y(j))*(M(j,3)-meanM(j))/meanM(j);
end


% REQUIRED??
% dM=dM-repmat(min(dM')',1,3);

% vector projected onto
x1=ones(length(M),2);
y1=ones(length(M),2);
z1=ones(length(M),2);

% the vectors point straight at x (resp. y) at each point
for i=1:length(M)
    x1(i,1)=1-x(i);
    x1(i,2)=-y(i);
    y1(i,1)=-x(i);
    y1(i,2)=1-y(i);
    z1(i,1)=-x(i);
    z1(i,2)=-y(i);
    
    x1(i,:)=x1(i,:)./norm(x1(i,:));
    y1(i,:)=y1(i,:)./norm(y1(i,:));
    z1(i,:)=z1(i,:)./norm(z1(i,:));
end

% project the data onto vectors pointing at x, y, or z=1-x-y
u=zeros(length(M),2);
% u(:,1)=x1(:,1).*dM(:,1)+y1(:,1).*dM(:,2)+z1(:,1).*dM(:,3);
% u(:,2)=x1(:,2).*dM(:,1)+y1(:,2).*dM(:,2)+z1(:,2).*dM(:,3);
u(:,1)=dM(:,1);
u(:,2)=dM(:,2);

% plot all
figure
% contour(x,y,z);
hold on
quiver(x,y,u(:,1),u(:,2),2/3)
plot([0 1],[1 0],'k')
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
axis square
hold off

% to plot evolutionary path (running average):
% plot(smooth(fig1(1,:)/nindi,20),smooth(fig1(2,:)/nindi,20))
