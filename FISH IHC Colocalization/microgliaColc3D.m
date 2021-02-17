function [colocTable] = microgliaColc3D(fname, Options)
%% RNA Scope 3D Colocalization with Micrglia IHC
% Written by Vasu Shandar & Kevin Coffey 2021
% 
%   INPUT: RGB Image Stack
%          Red Channel: FISH Signal
%          Green Channel: Microglia IHC
%          Blue Channel: DAPI
%
%  Options:
%% Options Defaults 
if ~isfield(Options, 'fiberSize') 
   Options.fiberSize = 12;
end

if ~isfield(Options, 'fishThreshold') 
   Options.fishThreshold = .185;
end

if ~isfield(Options, 'fishFilter')
   Options.fishFilter = 24;
end

if ~isfield(Options, 'ihcThreshold') 
   Options.ihcThreshold = .185;
end

if ~isfield(Options, 'ihcFilter')
   Options.ihcFilter = 5000;
end

if ~isfield(Options, 'makeFig3D') 
   Options.makeFig3D = 0;
end

if ~isfield(Options, 'figFold') 
   Options.figFold = 0;
end

%% Magic Code Stuff
tiff_info = imfinfo(fname); % return tiff structure, one element per image
im = imread(fname, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(fname, ii);
    im = cat(3 , im, temp_tiff);
end

for i =1 : (size(im,3))/3;
    G(:,:,i) = im(:,:,((i*3)-1));
    %imshow(G(:,:,i));
end 

for i =1 : (size(im,3))/3;
    RC(:,:,i) = im(:,:,((i*3)-2));
    %imshow(R(:,:,i));
end 

for i =1 : (size(im,3))/3;
    B(:,:,i) = im(:,:,((i*3)));
    %imshow(R(:,:,i));
end
% For Example Video
% v = VideoWriter('DNN.avi');
% open(v);
% Improves Fiber Detection
net = denoisingNetwork('DnCNN');   
for i = 1 : size(G,3);
bwG(:,:,i) = denoiseImage(G(:,:,i),net);
% bwG(:,:,i) = imadjust(bwG(:,:,i),[],[],1.25);
% imshowpair(G(:,:,i),bwG(:,:,i), 'montage');
% frame = getframe(gca);
% for ii=1:30
% writeVideo(v,frame);
% end
end
% bwG = fibermetric(bwG,'ObjectPolarity','bright','StructureSensitivity',8);

% close(v);

for i = 1 : size(RC,3);
bwR(:,:,i) = denoiseImage(RC(:,:,i),net);
% bwR(:,:,i) = imadjust(bwR(:,:,i),[],[],1.25);
% imshowpair(RC(:,:,i),bwR(:,:,i), 'montage');
end

for i = 1 : size(B,3);
bwB(:,:,i) = denoiseImage(B(:,:,i),net);
% bwB(:,:,i) = imadjust(bwB(:,:,i),[],[],1.25);
% imshowpair(B(:,:,i),bwB(:,:,i), 'montage');
end

BW_R = imbinarize(bwR,Options.fishThreshold);
BW_R = bwareaopen(BW_R, Options.fishFilter); % Removes all small objects from image

level = graythresh(bwG);
BW_G = imbinarize(bwG,level); % Not sure how to binarize...
se = strel('sphere',4);
BW_G = imclose(BW_G,se);
BW_G = bwareaopen(BW_G, Options.ihcFilter); % Removes all small objects from image 

level = graythresh(bwB);
BW_B = imbinarize(bwB,level);
se = strel('sphere',2);
BW_B = imclose(BW_B,se);
BW_B = bwareaopen(BW_B, 500); % Removes all small objects from image 

ConnectedComponents=bwconncomp(BW_G,26); %returns structure with 4 fields. PixelIdxList contains a 1-by-NumObjects cell array where the k-th element in the cell array is a vector containing the linear indices of the pixels in the k-th object. 26 defines connectivity. This looks at cube of connectivity around pixel.
stats = regionprops3(ConnectedComponents); % Some stats about the cells
numObj = numel(ConnectedComponents.PixelIdxList); %PixelIdxList is field with list of pixels in each connected component. Find how many connected components there are.
s = size(G);

if Options.makeFig3D~=0
% Simple 3D Visualization
    figure('Color','k','Position',[1 1 1024 1024]);
    set(gca,'Color',[0 0 0],'DataAspectRatio',[.5 .5 .5],'XColor','w','YColor','w','XLim',[1 s(1)],'YLim',[1 s(2)]);
    view(0,270); % Lookf at image from top viewpoint instead of side  
    % daspect([1 1 1]); 
    fv=isosurface((BW_R),0);
    patch(fv,'FaceColor',[1 0 0],'FaceAlpha',.95,'EdgeColor','none');
    fv=isosurface((BW_B),0);
    patch(fv,'FaceColor',[0 0 1],'FaceAlpha',.35,'EdgeColor','none');
    fv=isosurface((BW_G),0);%display each object as a surface in 3D. Will automatically add the next object to existing image.
    patch(fv,'FaceColor',[0 1 0],'FaceAlpha',.35,'EdgeColor','none');
    camlight; lighting phong
    f= frame2im(getframe(gca));
end

if Options.figFold~=0
    figure('Color','k','Position',[1 1 2048 2048]);
    threeDMask=cat(3,imadjust(mean(BW_R,3)),imadjust(mean(BW_G,3)),imadjust(mean(BW_B,3)));
    maxProject=cat(3,max(RC,[],3),max(G,[],3),max(B,[],3));
    maxDenoise=cat(3,max(bwR,[],3),max(bwG,[],3),max(bwB,[],3));
    %v=cat(4,bwR,bwG,bwB);
    if Options.makeFig3D~=0
        montage({maxProject,maxDenoise,threeDMask,f});
    else
        montage({maxProject,maxDenoise,threeDMask});
    end
    splt=strsplit(fname,'\');
    export_fig(fullfile(Options.figFold,splt{end}));
end

% volumeViewer(bwG);
%% Cell Specific Colocalizaiton
for m = 1:numObj-1
    ex=zeros(s(1),s(2),s(3));
    ex(ConnectedComponents.PixelIdxList{1,m})=1;%write in only one object to image. Cells are white on black background
    redColoc(m,1)=sum(BW_R(logical(ex)));
    cellSize(m,1)=sum(ex,'all');
end

%% Final Table Creation
splt=strsplit(fname,'\');
tok=strtok(splt{end},'.');
ID=repmat(categorical({[tok(1) tok(5)]}),height(cellSize),1);
colocTable = table(ID,cellSize,redColoc);
close all;
end

