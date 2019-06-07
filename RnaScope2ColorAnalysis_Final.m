
clear all
close all

path='Z:\Neumaier Lab\Pet1 RNAseq\WildType RNA Scope\Images';
files=dir(path);
figpath='Z:\Neumaier Lab\Pet1 RNAseq\WildType RNA Scope';
% Pixels per Square Millimeter = 600000
files=files(3:end,:);

for i=1:length(files)
disp(['Image: ' num2str(i)]);
tmp=strsplit(files(i).name,'-');

% Image Info
Sample(i)=tmp(1,1);
tmp=strsplit(tmp{1,2},'_');
Slide(i)=tmp(1,1);
Slice(i)=tmp(1,2);

% Image Read
img=imread(fullfile(path,files(i).name));

% DAPI Cell Segmentation
imgblue = img(:,:,3) ; % blue, nuclei only signal
background = imopen(imgblue,strel('disk',15));
imgblue = imgblue - background;
bw_blue = imbinarize(imgblue);
bw_blue = bwareaopen(bw_blue, 50);
cc_blue = bwconncomp(bw_blue, 8);

% Create Dapi Centroids & Remove centroids touching the border
radius = 30;
s = regionprops(cc_blue,'centroid');
centroids = cat(1, s.Centroid);
NotAtBorder = centroids(:,1) < size(img,2) - (radius + 1) &...
centroids(:,1) > radius + 1 &...
centroids(:,2) < size(img,1) - (radius + 1) &...
centroids(:,2) > radius + 1;
centroids = centroids(NotAtBorder,:);  
cell_radious(1:length(centroids),1)=30;

% % Show Cell Identification
% % figure;
% % imshow(imgblue)
% % hold on
% % plot(centroids(:,1),centroids(:,2), 'b*')
% % hold on
% % viscircles(centroids,cell_radious,'LineStyle','--','LineWidth',1);

% Red Segmentation
imgred = img(:,:,1) ; % red, neuron only signal
background = imopen(imgred,strel('disk',15));
imgred=imgred-background;
%imgred = imadjust(imgred);
imgred = imbinarize(imgred,.1);
%imgred = bwareaopen(imgred, 50);

% % Show Red Channel Segmentation
% % figure;
% % imshow(imgred);

% Green Channel Segmentation
imggreen = img(:,:,2);
background = imopen(imggreen,strel('disk',15));
imggreen=imggreen-background;
%imggreen = imadjust(imggreen);
imggreen = imbinarize(imggreen,.1);

% % Show Red Channel Segmentation
% % figure;
% % imshow(imggreen);

% Green & Red per Cell
SE = strel('disk',radius,0);
circle = find(SE.Neighborhood);
[I,J] = ind2sub(size(SE.Neighborhood),circle);
I = I - radius - 1;
J = J - radius - 1;
dist = sqrt(J.^2 + I.^2);
[dist,distIdx] = sort(dist,'ascend');
I = I(distIdx);
J = J(distIdx);
FEV_size = zeros(length(centroids),1);
FKBP5_Size = zeros(length(centroids),1);
FEV = zeros(length(centroids),length(I),'uint16');
FKBP5 = zeros(length(centroids),length(I),'uint16');

for j=1:length(centroids)
cellIndex = sub2ind(size(imgred),round(centroids(j,2)) + I,round(centroids(j,1)) + J);
FEV(j,:) = imgred(cellIndex);
FKBP5(j,:)= imggreen(cellIndex);
FEV_size(j) = sum(imgred(cellIndex));
FKBP5_Size(j)=sum(imggreen(cellIndex));
end

% Threshold Cells
RedCellThresh = 100;
GreenCellThresh = 100;
MergeCells=sum(FEV_size>RedCellThresh & FKBP5_Size>GreenCellThresh);
FKBP5Cells=sum(FKBP5_Size>GreenCellThresh & FEV_size<RedCellThresh);
FEVCells=sum(FEV_size>RedCellThresh & FKBP5_Size<GreenCellThresh);

% Final Numbers
PercentMergeCells(i,1)=MergeCells/(MergeCells+FEVCells);
FKBP5perFEVCell(i,1)=mean(FKBP5_Size(FEV_size>RedCellThresh));
FKBP5perOtherCell(i,1)=mean(FKBP5_Size(FEV_size<RedCellThresh));
FEVCellCounts(i,1)=(MergeCells+FEVCells);
FKBP5CellCounts(i,1)=FKBP5Cells;

% % Insert Circles
% I=insertShape(img,'circle',[centroids(FKBP5_Size>GreenCellThresh & FEV_size<RedCellThresh,:) cell_radious(FKBP5_Size>GreenCellThresh & FEV_size<RedCellThresh)],'Color','green','LineWidth',2);
% I=insertShape(I,'circle',[centroids(FEV_size>RedCellThresh & FKBP5_Size<GreenCellThresh,:) cell_radious(FEV_size>RedCellThresh & FKBP5_Size<GreenCellThresh)],'Color','red','LineWidth',2);
% I=insertShape(I,'circle',[centroids(FEV_size>RedCellThresh & FKBP5_Size>GreenCellThresh,:) cell_radious(FEV_size>RedCellThresh & FKBP5_Size>GreenCellThresh)],'Color','yellow','LineWidth',2);
% imwrite(I,fullfile(figpath,['Raw ' files(i).name]));

clear Final img imggreen imgred imgblue;
end

% Save Data Analysis
save('RNASCAOPE');
Sample=Sample';
Slide=Slide';
Slice=Slice';
MasterTable  = table(Sample,Slide,Slice,PercentMergeCells,FKBP5perFEVCell,FKBP5perOtherCell,FEVCellCounts,FKBP5CellCounts);
truth=readtable("Z:\Neumaier Lab\Pet1 RNAseq\WildType RNA Scope\Sample Truth Sheet.xlsx");
for i=1:height(MasterTable)
    IndexC = strfind(truth.WildTypeStressRNAScope3_13_2018, MasterTable.Sample(i));
    Index = find(not(cellfun('isempty', IndexC)));
    Animal(i,1)=truth.Var2(Index(1,1));
    Stress(i,1)=truth.Var3(Index(1,1));
end
Order=([1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3])';
MasterTable  = table(Sample,Slide,Slice,PercentMergeCells,FKBP5perFEVCell,FKBP5perOtherCell,FEVCellCounts,FKBP5CellCounts,Animal,Stress,Order);
writetable(MasterTable,fullfile(figpath,'FKBP5_Analysis.xlsx'));
MasterTable.Sample=categorical(MasterTable.Sample);

% Read Master Table
MasterTable=readtable("Z:\Neumaier Lab\Pet1 RNAseq\WildType RNA Scope\FKBP5_Analysis.xlsx");
sa = grpstats(MasterTable,{'Animal','Stress'},'mean','DataVars',{'PercentMergeCells','FKBP5perFEVCell','FKBP5perOtherCell'})
sa.Stress=categorical(sa.Stress);
MasterTable.Stress=categorical(MasterTable.Stress);
for i=1:height(sa);
if sa.Stress(i)=='FU' |  sa.Stress(i)=='MU'  
Stress1(i)=categorical({'US'});
else
Stress1(i)=categorical({'S'});    
end
if sa.Stress(i)=='FU' |  sa.Stress(i)=='FS'  
Sex(i)=categorical({'F'});
else
Sex(i)=categorical({'M'});    
end
end
sa.Stress1=(Stress1');
sa.Sex=Sex';
for i=1:height(MasterTable);
if MasterTable.Stress(i)=='FU' |  MasterTable.Stress(i)=='MU'  
Stress1(i)=categorical({'US'});
else
Stress1(i)=categorical({'S'});    
end
end
MasterTable.Stress1=(Stress1');

%% Final Plotting

% Fig 5c and 5d
g=gramm('x',sa.Stress1,'y',sa.mean_PercentMergeCells,'color',sa.Stress1);
g.stat_boxplot('width',1);
g.geom_jitter('width',.25);

figure('Position',[100 100 300 300]);
g.set_names('x','','y','% of Pet1 Cells w/ FKBP5');
g.axe_property('YLim',[0 1],'LineWidth',1.5,'FontSize',12);
g.set_order_options('x',{'US','S'});
g.draw();
g.results.geom_jitter_handle(1).MarkerEdgeColor=[0 0 0];
g.results.geom_jitter_handle(2).MarkerEdgeColor=[0 0 0];
g.results.stat_boxplot(1).box_handle.FaceAlpha=.25;
g.results.stat_boxplot(2).box_handle.FaceAlpha=.25;
g.export('file_name','Fig1c','file_type','png');

g=gramm('x',sa.Stress1,'y',sa.mean_FKBP5perFEVCell,'color',sa.Stress1);
g.stat_boxplot('width',1);
g.geom_jitter('width',.25);

figure('Position',[100 100 300 300]);
g.set_names('x','','y','FKBP5 Signal /Pet1 Cell');
g.axe_property('YLim',[0 550],'LineWidth',1.5,'FontSize',12);
g.set_order_options('x',{'US','S'});
g.draw();
g.results.geom_jitter_handle(1).MarkerEdgeColor=[0 0 0];
g.results.geom_jitter_handle(2).MarkerEdgeColor=[0 0 0];
g.results.stat_boxplot(1).box_handle.FaceAlpha=.25;
g.results.stat_boxplot(2).box_handle.FaceAlpha=.25;
g.export('file_name','Fig1d','file_type','png');
try
MU_FKIP=statarray.mean_PercentMergeCells(strcmp(statarray.Stress, 'MU'));
MS_FKIP=statarray.mean_PercentMergeCells(strcmp(statarray.Stress, 'MS'));
FU_FKIP=statarray.mean_PercentMergeCells(strcmp(statarray.Stress, 'FU'));
FS_FKIP=statarray.mean_PercentMergeCells(strcmp(statarray.Stress, 'FS'));

StressMean=mean([MS_FKIP;FS_FKIP]);
UnstressMean=mean([MU_FKIP;FU_FKIP]);
StressSEM=std([MS_FKIP;FS_FKIP])/sqrt(length([MS_FKIP;FS_FKIP]));
UnstressSEM=std([MU_FKIP;FU_FKIP])/sqrt(length([MU_FKIP;FU_FKIP]));

figure1 = figure('Color',[1 1 1],'Position',[0 0 300 300]);
cdat=[140/255 178/255 216/255;254/255 69/255 0/255];
% Create axes
axes1 = axes('Parent',figure1,...
    'TickDir','out',...
    'FontSize',12,...
    'LineWidth',1,...
    'Box','Off');
hold(axes1,'on');
% % vargout = barwitherrCat(errors,varargin,Cat)
[b2 eb2]=barwitherr([UnstressSEM StressSEM],[UnstressMean StressMean])
xlim([0 3]);
xticks([1 2]);
xticklabels({'US' 'S'});
ylabel('% of FEV Cells w/ FKBP5');
set(b2,'FaceColor','flat','EdgeAlpha',0);
set(eb2,'LineWidth',2);
b2.CData(1,:) = cdat(1,:);
b2.CData(2,:) = cdat(2,:);
export_fig(fullfile(figpath,'Percent of FEV Cells with FKBP5.png'),'-m3');
close all

[h,p,ci,stats] = ttest2([MU_FKIP;FU_FKIP],[MS_FKIP;FS_FKIP])
save('Stats_Fig5c','h','p','ci','stats')

MU_FKIP=statarray.mean_FKBP5perFEVCell(strcmp(statarray.Stress, 'MU'));
MS_FKIP=statarray.mean_FKBP5perFEVCell(strcmp(statarray.Stress, 'MS'));
FU_FKIP=statarray.mean_FKBP5perFEVCell(strcmp(statarray.Stress, 'FU'));
FS_FKIP=statarray.mean_FKBP5perFEVCell(strcmp(statarray.Stress, 'FS'));

StressMean=mean([MS_FKIP;FS_FKIP]);
UnstressMean=mean([MU_FKIP;FU_FKIP]);
StressSEM=std([MS_FKIP;FS_FKIP])/sqrt(length([MS_FKIP;FS_FKIP]));
UnstressSEM=std([MU_FKIP;FU_FKIP])/sqrt(length([MU_FKIP;FU_FKIP]));

figure1 = figure('Color',[1 1 1],'Position',[0 0 300 300]);
cdat=[140/255 178/255 216/255;254/255 69/255 0/255];
% Create axes
axes1 = axes('Parent',figure1,...
    'TickDir','out',...
    'FontSize',12,...
    'LineWidth',1,...
    'Box','Off');
hold(axes1,'on');
% % vargout = barwitherrCat(errors,varargin,Cat)
[b2 eb2]=barwitherr([UnstressSEM StressSEM],[UnstressMean StressMean])
xlim([0 3]);
xticks([1 2]);
xticklabels({'US' 'S'});
ylabel('FKBP5 Grains / Fev Cell');
set(b2,'FaceColor','flat','EdgeAlpha',0);
set(eb2,'LineWidth',2);
b2.CData(1,:) = cdat(1,:);
b2.CData(2,:) = cdat(2,:);
export_fig(fullfile(figpath,'FKBP5 Grains per Fev Cell.png'),'-m3');
close all

[h,p,ci,stats] = ttest2([MU_FKIP;FU_FKIP],[MS_FKIP;FS_FKIP])
save('Stats_Fig5d','h','p','ci','stats')
end

% Fig 5e
g=gramm('x',sa.Stress1,'y',sa.mean_FKBP5perOtherCell,'color',sa.Stress1);
g.stat_boxplot('width',1);
g.geom_jitter('width',.25);

figure('Position',[100 100 300 300]);
g.set_names('x','','y','FKBP5 Signal /FKBP5 Cell');
g.axe_property('YLim',[0 175],'LineWidth',1.5,'FontSize',12);
g.set_order_options('x',{'US','S'});
g.draw();
g.results.geom_jitter_handle(1).MarkerEdgeColor=[0 0 0];
g.results.geom_jitter_handle(2).MarkerEdgeColor=[0 0 0];
g.results.stat_boxplot(1).box_handle.FaceAlpha=.25;
g.results.stat_boxplot(2).box_handle.FaceAlpha=.25;
g.export('file_name','Fig1e','file_type','png');
try
MU_FKIP=statarray.mean_FKBP5perOtherCell(strcmp(statarray.Stress, 'MU'));
MS_FKIP=statarray.mean_FKBP5perOtherCell(strcmp(statarray.Stress, 'MS'));
FU_FKIP=statarray.mean_FKBP5perOtherCell(strcmp(statarray.Stress, 'FU'));
FS_FKIP=statarray.mean_FKBP5perOtherCell(strcmp(statarray.Stress, 'FS'));

StressMean=mean([MS_FKIP;FS_FKIP]);
UnstressMean=mean([MU_FKIP;FU_FKIP]);
StressSEM=std([MS_FKIP;FS_FKIP])/sqrt(length([MS_FKIP;FS_FKIP]));
UnstressSEM=std([MU_FKIP;FU_FKIP])/sqrt(length([MU_FKIP;FU_FKIP]));

figure1 = figure('Color',[1 1 1],'Position',[0 0 300 300]);
cdat=[140/255 178/255 216/255;254/255 69/255 0/255];
% Create axes
axes1 = axes('Parent',figure1,...
    'TickDir','out',...
    'FontSize',12,...
    'LineWidth',1,...
    'Box','Off');
hold(axes1,'on');
% % vargout = barwitherrCat(errors,varargin,Cat)
[b2 eb2]=barwitherr([UnstressSEM StressSEM],[UnstressMean StressMean])
xlim([0 3]);
xticks([1 2]);
xticklabels({'US' 'S'});
ylabel('FKBP5 Grains / Fev Cell');
set(b2,'FaceColor','flat','EdgeAlpha',0);
set(eb2,'LineWidth',2);
b2.CData(1,:) = cdat(1,:);
b2.CData(2,:) = cdat(2,:);
export_fig(fullfile(figpath,'FKBP5 Grains per Fev Cell.png'),'-m3');
close all

[h,p,ci,stats] = ttest2([MU_FKIP;FU_FKIP],[MS_FKIP;FS_FKIP])
save('Stats_Fig5e','h','p','ci','stats')
end

% Fig 5f 
g=gramm('x',sa.Sex,'y',sa.mean_PercentMergeCells,'color',sa.Stress1);
g.stat_boxplot('width',.65,'dodge',.65);
g.geom_jitter('width',.25,'dodge',.65);

figure('Position',[100 100 300 300]);
g.set_names('x','','y','% of Pet1 Cells w/ FKBP5');
g.axe_property('YLim',[0 1],'LineWidth',1.5,'FontSize',12);
g.set_order_options('x',{'F','M'},'color',{'US','S'});
g.draw();
g.results.geom_jitter_handle(1).MarkerEdgeColor=[0 0 0];
g.results.geom_jitter_handle(2).MarkerEdgeColor=[0 0 0];
g.results.stat_boxplot(1).box_handle.FaceAlpha=.25;
g.results.stat_boxplot(2).box_handle.FaceAlpha=.25;
g.export('file_name','Fig1f','file_type','png');

%fig 5g
g=gramm('x',sa.Sex,'y',sa.mean_FKBP5perFEVCell,'color',sa.Stress1);
g.stat_boxplot('width',.65,'dodge',.65);
g.geom_jitter('width',.25,'dodge',.65);

figure('Position',[100 100 300 300]);
g.set_names('x','','y','FKBP5 Signal /Pet1 Cell');
g.axe_property('YLim',[0 550],'LineWidth',1.5,'FontSize',12);
g.set_order_options('x',{'F','M'},'color',{'US','S'});
g.draw();
g.results.geom_jitter_handle(1).MarkerEdgeColor=[0 0 0];
g.results.geom_jitter_handle(2).MarkerEdgeColor=[0 0 0];
g.results.stat_boxplot(1).box_handle.FaceAlpha=.25;
g.results.stat_boxplot(2).box_handle.FaceAlpha=.25;
g.export('file_name','Fig1g','file_type','png');
try
MU_FKIP=statarray.mean_PercentMergeCells(strcmp(statarray.Stress, 'MU'));
MS_FKIP=statarray.mean_PercentMergeCells(strcmp(statarray.Stress, 'MS'));
FU_FKIP=statarray.mean_PercentMergeCells(strcmp(statarray.Stress, 'FU'));
FS_FKIP=statarray.mean_PercentMergeCells(strcmp(statarray.Stress, 'FS'));

mStressMean=mean([MS_FKIP]);
mUnstressMean=mean([MU_FKIP]);
mStressSEM=std([MS_FKIP])/sqrt(length(MS_FKIP));
mUnstressSEM=std([MU_FKIP])/sqrt(length(MU_FKIP));

fStressMean=mean([FS_FKIP]);
fUnstressMean=mean([FU_FKIP]);
fStressSEM=std([FS_FKIP])/sqrt(length(FS_FKIP));
fUnstressSEM=std([FU_FKIP])/sqrt(length(FU_FKIP));

figure1 = figure('Color',[1 1 1],'Position',[0 0 300 300]);
cdat=[140/255 178/255 216/255;254/255 69/255 0/255];
% Create axes
axes1 = axes('Parent',figure1,...
    'TickDir','out',...
    'FontSize',12,...
    'LineWidth',1,...
    'Box','Off');
% % vargout = barwitherrCat(errors,varargin,Cat)
[b2 eb2]=barwitherr([fUnstressSEM fStressSEM;mUnstressSEM mStressSEM],[fUnstressMean fStressMean;mUnstressMean mStressMean])
box off;
xticks([1 2]);
xticklabels({'FUS   FS ' 'MUS   MS '});
ylabel('% of FEV Cells w/ FKBP5');
set(b2(1),'FaceColor',cdat(1,:),'EdgeAlpha',0);
set(b2(2),'FaceColor',cdat(2,:),'EdgeAlpha',0);
set(eb2,'LineWidth',2);
export_fig(fullfile(figpath,'Sex Percent of FEV Cells with FKBP5.png'),'-m3');

[h,p,ci,stats] = ttest2([MU_FKIP],[MS_FKIP])
save('Stats_Fig5f_Male','h','p','ci','stats')
[h,p,ci,stats] = ttest2([FU_FKIP],[FS_FKIP])
save('Stats_Fig5f_female','h','p','ci','stats')
close all

% % Male v Female Stress Effects
% Overall Stress Effects
MU_FKIP=statarray.mean_FKBP5perFEVCell(strcmp(statarray.Stress, 'MU'));
MS_FKIP=statarray.mean_FKBP5perFEVCell(strcmp(statarray.Stress, 'MS'));
FU_FKIP=statarray.mean_FKBP5perFEVCell(strcmp(statarray.Stress, 'FU'));
FS_FKIP=statarray.mean_FKBP5perFEVCell(strcmp(statarray.Stress, 'FS'));

mStressMean=mean([MS_FKIP]);
mUnstressMean=mean([MU_FKIP]);
mStressSEM=std([MS_FKIP])/sqrt(length(MS_FKIP));
mUnstressSEM=std([MU_FKIP])/sqrt(length(MU_FKIP));

fStressMean=mean([FS_FKIP]);
fUnstressMean=mean([FU_FKIP]);
fStressSEM=std([FS_FKIP])/sqrt(length(FS_FKIP));
fUnstressSEM=std([FU_FKIP])/sqrt(length(FU_FKIP));

figure1 = figure('Color',[1 1 1],'Position',[0 0 300 300]);
cdat=[140/255 178/255 216/255;254/255 69/255 0/255];
% Create axes
axes1 = axes('Parent',figure1,...
    'TickDir','out',...
    'FontSize',12,...
    'LineWidth',1,...
    'Box','Off');
% % vargout = barwitherrCat(errors,varargin,Cat)
[b2 eb2]=barwitherr([fUnstressSEM fStressSEM;mUnstressSEM mStressSEM],[fUnstressMean fStressMean;mUnstressMean mStressMean])
box off;
xticks([1 2]);
xticklabels({'FUS   FS ' 'MUS   MS '});
ylabel('FKBP5 Grains / Fev Cell');
set(b2(1),'FaceColor',cdat(1,:),'EdgeAlpha',0);
set(b2(2),'FaceColor',cdat(2,:),'EdgeAlpha',0);
set(eb2,'LineWidth',2);
export_fig(fullfile(figpath,'Sex FKBP5 Grains per Fev Cell.png'),'-m3');
[h,p,ci,stats] = ttest2([MS_FKIP],[MU_FKIP])
save('Stats_Fig5g_Male','h','p','ci','stats')
[h,p,ci,stats] = ttest2([FS_FKIP],[FU_FKIP])
save('Stats_Fig5g_female','h','p','ci','stats')

close all
end

% Positional Effects Fig 5h
g=gramm('x',MasterTable.Order,'y',MasterTable.FKBP5perFEVCell,'color',MasterTable.Stress1);
g.stat_boxplot('width',.35,'dodge',.35);
g.geom_jitter('width',.15,'dodge',.35);

figure('Position',[100 100 300 300]);
g.set_names('x','','y','FKBP5 Signal /Pet1 Cell');
g.axe_property('YLim',[0 550],'LineWidth',1.5,'FontSize',12);
g.set_order_options('x',{1 2 3},'color',{'US','S'});
g.draw();
g.results.geom_jitter_handle(1).MarkerEdgeColor=[0 0 0];
g.results.geom_jitter_handle(2).MarkerEdgeColor=[0 0 0];
g.results.stat_boxplot(1).box_handle.FaceAlpha=.25;
g.results.stat_boxplot(2).box_handle.FaceAlpha=.25;
g.export('file_name','Fig1h','file_type','png');
try
% Positional Stress Effects

U1_FKIP=MasterTable.FKBP5perFEVCell((strcmp(MasterTable.Stress, 'MU') | strcmp(MasterTable.Stress, 'FU')) & MasterTable.Order==1);
S1_FKIP=MasterTable.FKBP5perFEVCell((strcmp(MasterTable.Stress, 'MS') | strcmp(MasterTable.Stress, 'FS')) & MasterTable.Order==1);
U2_FKIP=MasterTable.FKBP5perFEVCell((strcmp(MasterTable.Stress, 'MU') | strcmp(MasterTable.Stress, 'FU')) & MasterTable.Order==2);
S2_FKIP=MasterTable.FKBP5perFEVCell((strcmp(MasterTable.Stress, 'MS') | strcmp(MasterTable.Stress, 'FS')) & MasterTable.Order==2);
U3_FKIP=MasterTable.FKBP5perFEVCell((strcmp(MasterTable.Stress, 'MU') | strcmp(MasterTable.Stress, 'FU')) & MasterTable.Order==3);
S3_FKIP=MasterTable.FKBP5perFEVCell((strcmp(MasterTable.Stress, 'MS') | strcmp(MasterTable.Stress, 'FS')) & MasterTable.Order==3);

StressMean1=mean(S1_FKIP);
UnstressMean1=mean(U1_FKIP);
StressSEM1=std(S1_FKIP)/sqrt(length(S1_FKIP));
UnstressSEM1=std(U1_FKIP)/sqrt(length(U1_FKIP));

StressMean2=mean(S2_FKIP);
UnstressMean2=mean(U2_FKIP);
StressSEM2=std(S2_FKIP)/sqrt(length(S2_FKIP));
UnstressSEM2=std(U2_FKIP)/sqrt(length(U2_FKIP));

StressMean3=mean(S3_FKIP);
UnstressMean3=mean(U3_FKIP);
StressSEM3=std(S3_FKIP)/sqrt(length(S3_FKIP));
UnstressSEM3=std(U3_FKIP)/sqrt(length(U3_FKIP));

figure1 = figure('Color',[1 1 1],'Position',[0 0 400 300]);
cdat=[140/255 178/255 216/255;254/255 69/255 0/255];
% Create axes
axes1 = axes('Parent',figure1,...
    'TickDir','out',...
    'FontSize',12,...
    'LineWidth',1,...
    'Box','Off');
% % vargout = barwitherrCat(errors,varargin,Cat)
[b2 eb2]=barwitherr([UnstressSEM1 StressSEM1;UnstressSEM2 StressSEM2;UnstressSEM3 StressSEM3],[UnstressMean1 StressMean1;UnstressMean2 StressMean2;UnstressMean3 StressMean3])
box off;
xticks([1 2 3]);
xticklabels({'-4.25' '-4.5' '-4.75'});
ylabel('FKBP5 Grains / Fev Cell');
set(b2(1),'FaceColor',cdat(1,:),'EdgeAlpha',0);
set(b2(2),'FaceColor',cdat(2,:),'EdgeAlpha',0);
set(eb2,'LineWidth',2);
%l=legend(b2,'Location','northeastoutside','US','S');
set(l,'Box','off');
export_fig(fullfile(figpath,'Position FKBP5 Grains per Fev Cell5.png'),'-m3');

[h,p,ci,stats] = ttest2([S1_FKIP],[U1_FKIP])
save('Stats_Fig5h1','h','p','ci','stats')
[h,p,ci,stats] = ttest2([S2_FKIP],[U2_FKIP])
save('Stats_Fig5h2','h','p','ci','stats')
[h,p,ci,stats] = ttest2([S3_FKIP],[U3_FKIP])
save('Stats_Fig5h3','h','p','ci','stats')

close all
end



