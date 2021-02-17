%% Script to batch process IHC FISH Overlap Images for Morphine Microglia Project 2016-2021........ 
close all
clear all

%% Options
Options.fishThreshold = .175;
Options.fishFilter = 20;
Options.ihcThreshold = .15;
Options.ihcFilter = 10000;
Options.makeFig3D = 0;
Options.figFold = 'Images';

d=dir('Images');
d=d(3:end,:);

for i=1:height(d);
    fname=fullfile(d(i).folder,d(i).name);
    [colocTable] = microgliaColc3D(fname, Options);
    if i==1;
        mT=colocTable;
    else
        mT=[mT; colocTable];
    end
end