% Sets Directory For Data & Figures
directory=uigetdir('D:\Neumaier Lab\Pet1 RNAseq\Figure Creation Folders\For Combined'); % Raw data
fig_fold='D:\Neumaier Lab\Pet1 RNAseq\Figure Creation Folders\For Combined'; % Figure folder
cd(directory);
files=dir('*DESeq2*'); %find FST data files
%load('DRGenes.mat');


% Generate cmap
% lightness, chroma, hue range
lightness = [65, 65];
chroma = [75, 75];
hue = [205 385];

% colormap resolution
n = 100;
LHC = [
    linspace(lightness(1),lightness(2),n)
    linspace(chroma(1),chroma(2),n)
    linspace(hue(1),hue(2),n)
    ]';
cmap = pa_LCH2RGB(LHC);


% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 1000 600]);
ha = tight_subplot(2,3,[.05 .05],[.12 .12],[.12 .12]) ;

for l=1:6
tmp=importdata(files(l).name); % Import seq data
filename=strtok(files(l).name,'.')    
C = strsplit(filename,'_');
tmp.textdata=tmp.textdata(2:end,1);
% G1_name=C(2);
% G2_name=C(3);
% tmp.data(isnan(tmp.data(:,6)),6)=1;


axes(ha(l));
set(ha,'FontSize',12,'LineWidth',1,'TickDir','out');
hold on;
scatter(log2(tmp.data(:,1)),tmp.data(:,2),6,[.5 .5 .5],'filled','o');
idx=(tmp.data(:,6))<.20;
scatter(log2(tmp.data(idx,1)),tmp.data(idx,2),12,tmp.data(idx,6),'filled','o');

colormap(flipud(cmap))
caxis([0 .2])
ylim([-1.5 1.5]);
xlim([-5 15]);

if l==1
yticks([-1 0 1 ])
yticklabels({'-1','0','1'})
ylabel('log_2(Fold Change)');
text(-4,1.5,'Male IP','FontSize',12)
% Create xlabel
%xlabel('log2(Mean TPM)','FontWeight','bold');
end

if l==2
    yticks([-1 0 1 ])
text(-4,1.5,'Female IP','FontSize',12)
end

if l==3
    yticks([-1 0 1 ])
text(-4,1.5,'All IP','FontSize',12)
end


if l==4
yticks([-1 0 1])
yticklabels({'-1','0','1'})
xticks([-5 0 5 10 15])
xticklabels({'-5','0','5','10','15'});
ylabel('log_2(Fold Change)');
% Create xlabel
xlabel('log_2(TPM)');
text(-4,1.5,'Male Input','FontSize',12)
end

if l==5
xticks([-5 0 5 10 15])
yticks([-1 0 1 ])
xticklabels({'-5','0','5','10','15'});
xlabel('log_2(TPM)');
text(-4,1.5,'Female Input','FontSize',12)
end

if l==6
xticks([-5 0 5 10 15])
yticks([-1 0 1 ])
xticklabels({'-5','0','5','10','15'});
xlabel('log_2(TPM)');
text(-4,1.5,'All Input','FontSize',12)
end

clear G1_name G2_name
end

pngFileName = 'Stress DifferentialExpression.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure

% Plotting colorbar Expression 
f2=figure('color','w','position',[100 100 300 600]);
set(gca,'FontSize',12,'LineWidth',1,'TickDir','out');
h = colorbar('Box','off');
caxis([0 .2]);
colormap(flipud(cmap))
set(get(h,'title'),'string','FDR');pngFileName = 'Colorbar.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure

% Piechart
f1=figure('color','w','position',[100 100 1000 600]);
set(gca,'FontSize',16,'FontWeight','bold');
for l=1:6
tmp=importdata(files(l).name); % Import seq data
filename=strtok(files(l).name,'.');    
C = strsplit(filename,'_');
tmp.textdata=tmp.textdata(2:end,1);  
if l==1
   SFG=tmp.textdata(tmp.data(:,6)<=.2 & tmp.data(:,2)>0,1);
elseif l==2
   SMG=tmp.textdata(tmp.data(:,6)<=.2 & tmp.data(:,2)>0,1); 
elseif l==3
   UFG=tmp.textdata(tmp.data(:,6)<=.2 & tmp.data(:,2)>0,1);  
elseif l==4
   UMG=tmp.textdata(tmp.data(:,6)<=.2 & tmp.data(:,2)>0,1);   
end
end

AllGenes=[SFG;SMG;UFG;UFG];
UniqGenes=unique(AllGenes);
SFG_i=intersect(unique([SMG;UFG;UFG]),SFG);
SFG_un=length(SFG)-length(SFG_i);
SMG_i=intersect(unique([SFG;UFG;UFG]),SMG);
SMG_un=length(SMG)-length(SMG_i);
UFG_i=intersect(unique([SMG;SFG;UMG]),UFG);
UFG_un=length(UFG)-length(UFG_i);
UMG_i=intersect(unique([SMG;SFG;UFG]),UMG);
UMG_un=length(UMG)-length(UMG_i);
sumUnGenes=sum([SFG_un SMG_un UFG_un UMG_un]);
X=[(length(UniqGenes)-sumUnGenes), SFG_un, SMG_un, UFG_un, UMG_un];
h = pie(X);
colormap(jet(5));
hText = findobj(h,'Type','text'); % text object handles
percentValues = get(hText,'String'); % percent values
txt = {'Overlap (1761): ';'Stressed F (252): ';'Stressed M (245): ';'Unstressed F (46): ';'Unstressed M (307): '}; % strings
combinedtxt = strcat(txt,percentValues); % strings and percent values
oldExtents_cell = get(hText,'Extent');
oldExtents = cell2mat(oldExtents_cell); % numeric array
hText(1).String = combinedtxt(1);
hText(2).String = combinedtxt(2);
hText(3).String = combinedtxt(3);
hText(4).String = combinedtxt(4);
hText(5).String = combinedtxt(5);
newExtents_cell = get(hText,'Extent'); % cell array
newExtents = cell2mat(newExtents_cell); % numeric array 
width_change = newExtents(:,3)-oldExtents(:,3);
signValues = sign(oldExtents(:,1));
offset = signValues.*(width_change/2);
textPositions_cell = get(hText,{'Position'}); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array
textPositions(:,1) = textPositions(:,1) + offset; % add offset 
hText(1).Position = textPositions(1,:);
hText(2).Position = textPositions(2,:);
hText(3).Position = textPositions(3,:);
hText(4).Position = textPositions(4,:);
hText(5).Position = textPositions(5,:);
set(hText,'FontWeight','bold','FontSize',13);
set(h,'LineWidth',2);
pngFileName = 'ControlPie.png'; % Set the File name 
fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
export_fig(fullFileName, '-m5'); % Save the Figure

xlswrite('SignalGenes.xlsx',UniqGenes);

%% Signal Genes Only Analysis

signal=importdata('SignalGenes.xlsx');

for l=3

tmp=importdata(files(l).name); % Import seq data
filename=strtok(files(l).name,'.')    
C = strsplit(filename,'_');
tmp.textdata=tmp.textdata(2:end,1);

[genes ia ib]=intersect(tmp.textdata,signal.Sheet1);

% Signal Only Waldt
SigGene=tmp.textdata(ia,1);
Waldt=tmp.data(ia,4);

xlswrite('GSEAsignalonly.xlsx',SigGenes);

% P Value Distrobutions

AllP=tmp.data(:,5);

SigP=tmp.data(ia,5);

figure
x=hist(AllP,100)
bar(x)

figure
y=hist(SigP,100)
bar(y)

end