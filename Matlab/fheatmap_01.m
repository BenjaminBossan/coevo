% this function accepts as first input argument the matrix that contains
% the values to be shown in the "heatmap".
% the second argument is the image size in pixels (400 by default)
% the third argument is whether the scale is absolute (0->1) or relative
% (min->max); the default setting is absolute
% the 4th argument is the spacing of the tick marks
% the 5th and 6th argument are the x and y axis as vectors; note that you
% have to adjust them to the spacing
% the 7th and 8th argument are the labels of the x and y axis, respectively

function a=fheatmap_01(inputMat,imageSize,abso,spacing,...
    xticknumbers,yticknumbers,xaxisname,yaxisname);

if nargin<2
    imageSize=400;
    abso=1;
    spacing=1;
elseif nargin<3
    abso=1;
    spacing=1;
elseif nargin<4
    spacing=1;
end

if nargin<7
    xaxisname='weight of losses';
    yaxisname='weight of gains';
elseif nargin<8
    yaxisname='unspecified';
end

% inputMat=squeeze(inputMat);

[s1 s2]=size(inputMat);
binLength=ceil(imageSize/max(s1,s2));

if nargin<5
    xticknumbers=[-s1/2+.5:spacing:s1/2-.5];
    yticknumbers=[-s2/2+.5:spacing:s2/2-.5];
end

a=zeros(s1*binLength,s2*binLength);

for i=1:s1
    for j=1:s2
        a(1+(i-1)*binLength:i*binLength,1+(j-1)*binLength:j*binLength)=...
            inputMat(i,j);
    end
end

% a=fliplr(a);

figure

iptsetpref('ImshowAxesVisible','on')
axis square
set(gca,'fontsize',14)

xax=[0:1/s1/binLength:1]*max(s1)-s1/2;
yax=[0:1/s2/binLength:1]*max(s2)-s2/2;

if abso==1
    imshow(xax,yax,a,[-0.002 1.002]);
elseif abso==0
    imshow(xax,yax,a,[min(min(inputMat)) max(max(inputMat))]);
end

axis xy
axis square

colorbar
colorbar('ytick',[0:.2:1],'fontsize',14)
% colorbar()

colormap([1-[0:.01:1]' 1-[0:.01:1]' 1-[0:.01:1]']);

xlabel(xaxisname)
ylabel(yaxisname)

set(gca,'xtick',[-s1/2+.5:spacing:s1/2-.5])
set(gca,'ytick',[-s2/2+.5:spacing:s2/2-.5])
set(gca,'xticklabel',xticknumbers)
set(gca,'yticklabel',yticknumbers)
set(gca,'ticklength',[0 0])

hold on
if mod(s1,2)==1
    for i=floor(-s1/2):ceil(s1/2)-1
        plot([i+.5 i+.5],[floor(-s2/2)+.5 ceil(s2/2)-.5],'k')
    end
    for j=floor(-s2/2):ceil(s2/2)-1
        plot([floor(-s1/2)+.5 ceil(s1/2)-.5],[j+.5 j+.5],'k')
    end
else
    for i=floor(-s1/2):ceil(s1/2)
        plot([i i],[floor(-s2/2) ceil(s2/2)],'k')
    end
    for j=floor(-s2/2):ceil(s2/2)
        plot([floor(-s1/2) ceil(s1/2)],[j j],'k')
    end
end
