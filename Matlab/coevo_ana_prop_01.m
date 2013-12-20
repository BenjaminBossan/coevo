% this programm just computes the behavior of the strategies
% without any kind of evolution.
% The goal is to find out the PERFORMANCE of the different strategies
% This particular instance tests PBSL with a proportional choice rule
% 1) no observation is the same as a negative observation
% 2) a negative outcome counts 1/10 of a positive observation
% 3) no observation counts 1/10 of a positive observation

clear all

% PARAMETERS
tmax=100000                    % periods per generation
regime=1                   % how the environment changes. 0->no regression to the mean, 1->medium, 32->high
incr=2/100;                 % increment at which the environment becomes better or worse
pincr=1;                    % probability that environmental quality changes at all after each period
dpA=0;
dpB=0;                    % shift of pA and pB
pA0=1/2;pB0=1/2;
choiceLimit=1/10000        % puts a limit on the min/max proportion, so that they are not 0 or 1

% INITIALIZE
x=zeros(7,tmax);        % percent of A choices
y=x;z=x;
perfmatx=x;              % performance matrix
perfmaty=x;
perfmatz=x;
    
% ROUTINE
    
% ENVIRONMENT
% routine to determine pA and pB
[pA,pB] = randomenvironment4(tmax,regime,incr,pincr,pA0,pB0);
% possible shift dp
pA=min(1,max(0,pA+dpA));pB=min(1,max(pB+dpB,0));


tic

% first choice random
x(:,:,1)=1/2;
y(:,:,1)=1/2;
z(:,:,1)=1/2;

wb=waitbar(0,'progress');

% ROUTINE
for t=1:tmax-1
    
    x(1,t+1)=(1+pB(t)*(-1+x(1,t))+pA(t)*x(1,t))/2;
    x(2,t+1)=(1+2*pB(t)*(-1+x(2,t))+pB(t)^2*(-1+x(2,t))^2+pA(t)*x(2,t)*(2-pA(t)*x(2,t)))/2;
    x(3,t+1)=(1+pB(t)^3*(-1+x(3,t))^3-pB(t)^2*(-1+x(3,t))^2*(-3+pA(t)*x(3,t))-pB(t)*(-1+x(3,t))*(-3+pA(t)^2*x(3,t)^2)+pA(t)*x(3,t)*(3+pA(t)*x(3,t)*(-3+pA(t)*x(3,t))))/2;
    x(4,t+1)=(1+pB(t)^4*(-1+x(4,t))^4-2*pB(t)^3*(-1+x(4,t))^3*(-2+pA(t)*x(4,t))-2*pB(t)^2*(-1+x(4,t))^2*(-3+2*pA(t)*x(4,t))-pA(t)*x(4,t)*(-2+pA(t)*x(4,t))*(2+pA(t)*x(4,t)*(-2+pA(t)*x(4,t)))+2*pB(t)*(-1+x(4,t))*(2+pA(t)^2*x(4,t)^2*(-2+pA(t)*x(4,t))))/2;
    x(5,t+1)=(1+pB(t)^5*(-1+x(5,t))^5-pB(t)^4*(-1+x(5,t))^4*(-5+3*pA(t)*x(5,t))+2*pB(t)^2*(-1+x(5,t))^2*(5-5*pA(t)*x(5,t)+pA(t)^3*x(5,t)^3)+2*pB(t)^3*(-1+x(5,t))^3*(5+pA(t)*x(5,t)*(-5+pA(t)*x(5,t)))-pB(t)*(-1+x(5,t))*(-5+pA(t)^2*x(5,t)^2*(10+pA(t)*x(5,t)*(-10+3*pA(t)*x(5,t))))+pA(t)*x(5,t)*(5+pA(t)*x(5,t)*(-10+pA(t)*x(5,t)*(10+pA(t)*x(5,t)*(-5+pA(t)*x(5,t))))))/2;
    x(6,t+1)=(1+pB(t)^6*(-1+x(6,t))^6-2*pB(t)^5*(-1+x(6,t))^5*(-3+2*pA(t)*x(6,t))-pA(t)*x(6,t)*(-2+pA(t)*x(6,t))*(3+pA(t)*x(6,t)*(-3+pA(t)*x(6,t)))*(1+pA(t)*x(6,t)*(-1+pA(t)*x(6,t)))+2*pB(t)^3*(-1+x(6,t))^3*(10+3*pA(t)*x(6,t)*(-5+2*pA(t)*x(6,t)))+pB(t)^4*(-1+x(6,t))^4*(15+pA(t)*x(6,t)*(-18+5*pA(t)*x(6,t)))+2*pB(t)*(-1+x(6,t))*(3+pA(t)^2*x(6,t)^2*(-2+pA(t)*x(6,t))*(5+pA(t)*x(6,t)*(-5+2*pA(t)*x(6,t))))-pB(t)^2*(-1+x(6,t))^2*(-15+pA(t)*x(6,t)*(20+pA(t)^2*x(6,t)^2*(-12+5*pA(t)*x(6,t)))))/2;
    x(7,t+1)=(1+pB(t)^7*(-1+x(7,t))^7-pB(t)^6*(-1+x(7,t))^6*(-7+5*pA(t)*x(7,t))+pB(t)^5*(-1+x(7,t))^5*(21+pA(t)*x(7,t)*(-28+9*pA(t)*x(7,t)))-pB(t)^3*(-1+x(7,t))^3*(-35+pA(t)*x(7,t)*(70-42*pA(t)*x(7,t)+5*pA(t)^3*x(7,t)^3))-pB(t)^4*(-1+x(7,t))^4*(-35+pA(t)*x(7,t)*(63+5*pA(t)*x(7,t)*(-7+pA(t)*x(7,t))))+pB(t)^2*(-1+x(7,t))^2*(21+pA(t)*x(7,t)*(-35+pA(t)^2*x(7,t)^2*(42+pA(t)*x(7,t)*(-35+9*pA(t)*x(7,t)))))-pB(t)*(-1+x(7,t))*(-7+pA(t)^2*x(7,t)^2*(35+pA(t)*x(7,t)*(-70+pA(t)*x(7,t)*(63+pA(t)*x(7,t)*(-28+5*pA(t)*x(7,t))))))+pA(t)*x(7,t)*(7+pA(t)*x(7,t)*(-21+pA(t)*x(7,t)*(35+pA(t)*x(7,t)*(-35+pA(t)*x(7,t)*(21+pA(t)*x(7,t)*(-7+pA(t)*x(7,t))))))))/2;
    
    y(1,t+1)=y(1,t);
    y(2,t+1)=((11-9*pA(t)*(-1+y(2,t))+9*pB(t)*(-1+y(2,t)))*y(2,t))/11;
    y(3,t+1)=(y(3,t)*(14+3*pB(t)*(7+3*pB(t)*(-1+y(3,t)))*(-1+y(3,t))+9*pA(t)^2*(-1+y(3,t))*y(3,t)-3*pA(t)*(-1+y(3,t))*(7+pB(t)*(-3+6*y(3,t)))))/14;
    y(4,t+1)=(y(4,t)*(4433+27*pB(t)*(341+9*pB(t)*(31+9*pB(t)*(-1+y(4,t)))*(-1+y(4,t)))*(-1+y(4,t))-2187*pA(t)^3*(-1+y(4,t))*y(4,t)^2+243*pA(t)^2*(-1+y(4,t))*y(4,t)*(31+9*pB(t)*(-2+3*y(4,t)))-27*pA(t)*(-1+y(4,t))*(341+9*pB(t)*(-31+62*y(4,t)+9*pB(t)*(-1+y(4,t))*(-1+3*y(4,t))))))/4433;
    y(5,t+1)=(y(5,t)*(52808+9*pB(t)*(15088+27*pB(t)*(656+9*pB(t)*(41+9*pB(t)*(-1+y(5,t)))*(-1+y(5,t)))*(-1+y(5,t)))*(-1+y(5,t))+19683*pA(t)^4*(-1+y(5,t))*y(5,t)^3-2187*pA(t)^3*(-1+y(5,t))*y(5,t)^2*(41+9*pB(t)*(-3+4*y(5,t)))+243*pA(t)^2*(-1+y(5,t))*y(5,t)*(656+9*pB(t)*(27*pB(t)*(-1+y(5,t))*(-1+2*y(5,t))+41*(-2+3*y(5,t))))-9*pA(t)*(-1+y(5,t))*(15088+27*pB(t)*(656*(-1+2*y(5,t))+9*pB(t)*(-1+y(5,t))*(41*(-1+3*y(5,t))+9*pB(t)*(-1+y(5,t))*(-1+4*y(5,t)))))))/52808;
    y(6,t+1)=(y(6,t)*(2618+3*pB(t)*(2618+3*pB(t)*(1309+9*pB(t)*(119+3*pB(t)*(17+3*pB(t)*(-1+y(6,t)))*(-1+y(6,t)))*(-1+y(6,t)))*(-1+y(6,t)))*(-1+y(6,t))-729*pA(t)^5*(-1+y(6,t))*y(6,t)^4+243*pA(t)^4*(-1+y(6,t))*y(6,t)^3*(17+3*pB(t)*(-4+5*y(6,t)))-81*pA(t)^3*(-1+y(6,t))*y(6,t)^2*(119+3*pB(t)*(-51+68*y(6,t)+6*pB(t)*(-1+y(6,t))*(-3+5*y(6,t))))+9*pA(t)^2*(-1+y(6,t))*y(6,t)*(1309+9*pB(t)*(119*(-2+3*y(6,t))+9*pB(t)*(-1+y(6,t))*(-17+34*y(6,t)+2*pB(t)*(-1+y(6,t))*(-2+5*y(6,t)))))-3*pA(t)*(-1+y(6,t))*(2618+3*pB(t)*(1309*(-1+2*y(6,t))+9*pB(t)*(-1+y(6,t))*(119*(-1+3*y(6,t))+3*pB(t)*(-1+y(6,t))*(-17+68*y(6,t)+3*pB(t)*(-1+y(6,t))*(-1+5*y(6,t))))))))/2618;
    y(7,t+1)=(y(7,t)*(23187320+27*pB(t)*(2898415+9*pB(t)*(579683+9*pB(t)*(68198+27*pB(t)*(1586+9*pB(t)*(61+9*pB(t)*(-1+y(7,t)))*(-1+y(7,t)))*(-1+y(7,t)))*(-1+y(7,t)))*(-1+y(7,t)))*(-1+y(7,t))+4782969*pA(t)^6*(-1+y(7,t))*y(7,t)^5-531441*pA(t)^5*(-1+y(7,t))*y(7,t)^4*(61+9*pB(t)*(-5+6*y(7,t)))+59049*pA(t)^4*(-1+y(7,t))*y(7,t)^3*(1586+9*pB(t)*(45*pB(t)*(-1+y(7,t))*(-2+3*y(7,t))+61*(-4+5*y(7,t))))-4374*pA(t)^3*(-1+y(7,t))*y(7,t)^2*(34099+27*pB(t)*(793*(-3+4*y(7,t))+9*pB(t)*(-1+y(7,t))*(45*pB(t)*(-1+y(7,t))*(-1+2*y(7,t))+61*(-3+5*y(7,t)))))+243*pA(t)^2*(-1+y(7,t))*y(7,t)*(579683+9*pB(t)*(68198*(-2+3*y(7,t))+81*pB(t)*(-1+y(7,t))*(1586*(-1+2*y(7,t))+3*pB(t)*(-1+y(7,t))*(45*pB(t)*(-1+y(7,t))*(-1+3*y(7,t))+122*(-2+5*y(7,t))))))-27*pA(t)*(-1+y(7,t))*(2898415+9*pB(t)*(579683*(-1+2*y(7,t))+9*pB(t)*(-1+y(7,t))*(68198*(-1+3*y(7,t))+27*pB(t)*(-1+y(7,t))*(1586*(-1+4*y(7,t))+9*pB(t)*(-1+y(7,t))*(61*(-1+5*y(7,t))+9*pB(t)*(-1+y(7,t))*(-1+6*y(7,t)))))))))/23187320;
    
    z(1,t+1)=(11+10*pB(t)*(-1+z(1,t))+(-11+10*pA(t))*z(1,t))/11;
    z(2,t+1)=(66-55*pB(t)*(-2+z(2,t))*(-1+z(2,t))+50*pB(t)^2*(-1+z(2,t))^2+z(2,t)*(-66-50*pA(t)^2*z(2,t)+55*pA(t)*(1+z(2,t))))/66;
    z(3,t+1)=(3289+2000*pB(t)^3*(-1+z(3,t))^3-200*pB(t)^2*(-1+z(3,t))^2*(-33+(11+10*pA(t))*z(3,t))-10*pB(t)*(-1+z(3,t))*(-759+22*(23+10*pA(t))*z(3,t)+40*pA(t)*(-11+5*pA(t))*z(3,t)^2)+z(3,t)*(-3289+10*pA(t)*(253+2*z(3,t)*(253+100*pA(t)^2*z(3,t)-110*pA(t)*(2+z(3,t))))))/3289;
    z(4,t+1)=(374*(7+(-7+5*pA(t))*z(4,t))+5*(250*pB(t)^4*(-1+z(4,t))^4-25*pB(t)^3*(-1+z(4,t))^3*(-44+(11+20*pA(t))*z(4,t))+55*pB(t)^2*(-1+z(4,t))^2*(34+z(4,t)*(-17+5*pA(t)*(-5+3*z(4,t))))+pA(t)*z(4,t)^2*(1122-5*pA(t)*(187+z(4,t)*(187+50*pA(t)^2*z(4,t)-55*pA(t)*(3+z(4,t)))))+pB(t)*(-1+z(4,t))*(1496+z(4,t)*(-1122+5*pA(t)*(-187+z(4,t)*(374+5*pA(t)*(-22+(-33+20*pA(t))*z(4,t))))))))/2618;
    z(5,t+1)=(1155*(3+(-3+2*pA(t))*z(5,t))+2*(640*pB(t)^5*(-1+z(5,t))^5-64*pB(t)^4*(-1+z(5,t))^4*(-55+(11+30*pA(t))*z(5,t))+16*pB(t)^3*(-1+z(5,t))^3*(495-22*(9+22*pA(t))*z(5,t)+16*pA(t)*(11+5*pA(t))*z(5,t)^2)+8*pB(t)^2*(-1+z(5,t))^2*(1155+z(5,t)*(-693+2*pA(t)*(-693+66*(9+2*pA(t))*z(5,t)+8*pA(t)*(-33+10*pA(t))*z(5,t)^2)))+4*pA(t)*z(5,t)^2*(1155+2*pA(t)*(-462+z(5,t)*(-693+2*pA(t)*(297+2*z(5,t)*(99+20*pA(t)^2*z(5,t)-22*pA(t)*(4+z(5,t)))))))-pB(t)*(-1+z(5,t))*(-5775+4*z(5,t)*(1155+2*pA(t)*(693+2*z(5,t)*(-693+pA(t)*(99+2*z(5,t)*(297+60*pA(t)^2*z(5,t)-22*pA(t)*(7+4*z(5,t))))))))))/3465;
    z(6,t+1)=(69069*(8+(-8+5*pA(t))*z(6,t))+5*(31250*pB(t)^6*(-1+z(6,t))^6-3125*pB(t)^5*(-1+z(6,t))^5*(-66+(11+40*pA(t))*z(6,t))+625*pB(t)^4*(-1+z(6,t))^4*(924+z(6,t)*(-308+5*pA(t)*(-209+(55+50*pA(t))*z(6,t))))-2750*pB(t)^3*(-1+z(6,t))^3*(-322+z(6,t)*(161+5*pA(t)*(98+z(6,t)*(-56+5*pA(t)*(-8+5*z(6,t))))))+pB(t)*(-1+z(6,t))*(414414+5*z(6,t)*(-69069+5*pA(t)*(-21252+z(6,t)*(42504+5*pA(t)*z(6,t)*(-10626+5*pA(t)*(924+154*(8-5*pA(t))*z(6,t)+25*pA(t)*(-11+8*pA(t))*z(6,t)^2))))))-50*pB(t)^2*(-1+z(6,t))^2*(-15939+z(6,t)*(10626+5*pA(t)*(5313+z(6,t)*(-5313+5*pA(t)*(-462+z(6,t)*(924+5*pA(t)*(-33+5*(-11+5*pA(t))*z(6,t))))))))+5*pA(t)*z(6,t)^2*(69069+5*pA(t)*(-10626+z(6,t)*(-21252+5*pA(t)*(3542+z(6,t)*(3542-5*pA(t)*(616+z(6,t)*(308+50*pA(t)^2*z(6,t)-55*pA(t)*(5+z(6,t)))))))))))/552552;
    z(7,t+1)=(219160953*(17+(-17+10*pA(t))*z(7,t))+10*(80000000*pB(t)^7*(-1+z(7,t))^7-8000000*pB(t)^6*(-1+z(7,t))^6*(-77+(11+50*pA(t))*z(7,t))+400000*pB(t)^5*(-1+z(7,t))^5*(5159+2*z(7,t)*(-737+10*pA(t)*(-319+(66+90*pA(t))*z(7,t))))-40000*pB(t)^4*(-1+z(7,t))^4*(-98021+z(7,t)*(42009+10*pA(t)*(16951+10*z(7,t)*(-737+10*pA(t)*(-88+(33+10*pA(t))*z(7,t))))))-1000*pB(t)^3*(-1+z(7,t))^3*(-4606987+4*z(7,t)*(658141+10*pA(t)*(238051+4*z(7,t)*(-42009+5*pA(t)*(-8107+10*z(7,t)*(737+10*pA(t)*(11+(-22+5*pA(t))*z(7,t))))))))+20*pB(t)^2*(-1+z(7,t))^2*(170458519+5*z(7,t)*(-24351217+10*pA(t)*(-7239551+4*z(7,t)*(1974423+10*pA(t)*(126027+2*z(7,t)*(-126027+10*pA(t)*(737+10*z(7,t)*(737+90*pA(t)^2*z(7,t)-55*pA(t)*(5+3*z(7,t))))))))))+2*pA(t)*z(7,t)^2*(657482859+10*pA(t)*(-48702434+5*z(7,t)*(-24351217+10*pA(t)*(1974423+4*z(7,t)*(658141+10*pA(t)*(-56012+z(7,t)*(-42009+10*pA(t)*(3685+2*z(7,t)*(737+100*pA(t)^2*z(7,t)-110*pA(t)*(6+z(7,t)))))))))))-pB(t)*(-1+z(7,t))*(-1534126671+2*z(7,t)*(657482859+50*pA(t)*(24351217+2*z(7,t)*(-24351217+5*pA(t)*(-658141+4*z(7,t)*(1974423+10*pA(t)*(-70015+2*z(7,t)*(-84018+5*pA(t)*(9581+10*z(7,t)*(737+2*pA(t)*(-253+(-66+50*pA(t))*z(7,t))))))))))))))/3725736201;
    
    x(:,t+1)=min(x(:,t+1),1-choiceLimit);
    x(:,t+1)=max(x(:,t+1),choiceLimit);
    
    y(:,t+1)=min(y(:,t+1),1-choiceLimit);
    y(:,t+1)=max(y(:,t+1),choiceLimit);
    
    z(:,t+1)=min(z(:,t+1),1-choiceLimit);
    z(:,t+1)=max(z(:,t+1),choiceLimit);
    
    if mod(t,10)==0
        waitbar(t/tmax);
    end
end

close(wb)

% calculation of PERFORMANCE
ind=find(pA==pB); %indices of when pA and pB are equal (->draw);
delp=pA-pB; %delta p
delp(ind)=[];   %remove draws
delp=sign(delp);    % A or B better
x2=x;   %save complete form
y2=y;   %save complete form
z2=z;   %save complete form
x(:,ind)=[];  %remove draws from choices
y(:,ind)=[];  %remove draws from choices
z(:,ind)=[];  %remove draws from choices
perfmatx2=perfmatx; %save complete form
perfmaty2=perfmaty; %save complete form
perfmatz2=perfmatz; %save complete form
perfmatx(:,ind)=[];
perfmaty(:,ind)=[];
perfmatz(:,ind)=[];
for t=1:tmax-length(ind)
    perfmatx(:,t)=repmat((delp(t)>0),7,1).*x(:,t)+repmat((delp(t)<0),7,1).*(1-x(:,t));
    perfmaty(:,t)=repmat((delp(t)>0),7,1).*y(:,t)+repmat((delp(t)<0),7,1).*(1-y(:,t));
    perfmatz(:,t)=repmat((delp(t)>0),7,1).*z(:,t)+repmat((delp(t)<0),7,1).*(1-z(:,t));
end

'fail=no po -- fail>no po -- fail<no po'
[[1:7]' mean(perfmatx,2) mean(perfmaty,2) mean(perfmatz,2)]

'failue == no outcome'
[max(mean(perfmatx,2)) find(mean(perfmatx,2)==max(mean(perfmatx,2)))]
'failure > no outcome'
[max(mean(perfmaty,2)) find(mean(perfmaty,2)==max(mean(perfmaty,2)))]
'no outcome > failure'
[max(mean(perfmatz,2)) find(mean(perfmatz,2)==max(mean(perfmatz,2)))]



toc