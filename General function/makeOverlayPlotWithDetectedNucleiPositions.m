function makeOverlayPlotWithDetectedNucleiPositions(im,minArea,thresh_factor,maxArea,largeFilt,bounds)

if nargin<2
    minArea=15;
end
if nargin<3
    thresh_factor=25;
end
if nargin<4
    maxArea=Inf;
end
if nargin<5
    largeFilt=0;
end
if nargin<6
    bounds=[]; %for imshow
end

if isstr(im)
    im=imread(im);
    im=imageSubtractBackground(im);
    %im=imageSubtractBG_imopen(im);
end

%coors=getNucleiPositions(im,minArea,thresh_factor,maxArea,method);
coors=getNucleiPositions1(im,minArea,thresh_factor,maxArea,largeFilt);

figure;
imshow(im,bounds); hold on;
%scatter(coors(:,1),coors(:,2),coors(:,4));  %use the measured area as the circle area
scatter(coors(:,1),coors(:,2),11);  %use the measured area as the circle area
hold off;
%showImagesWithLinkedAxes({im,im>thresh_factor});
fprintf('Identified %i objects.\n',size(coors,1));
fprintf('Median area = %.1f pixels\n',median(coors(:,4)));
