% this file is meant to go through the steps of the analysis of the movies

% 1) segment the images on the HSF1 channel - Movie_Segmentation_v2.m
% 2) perform drift correction - imageregistration.m 
% 3) track cells in time - automatictracking.m 
% 4) extract the data

%% images segmentation for t_0
clear all

basefolder = 'Y:\IN_Cell_Analyzer_6000\Giorgio\2019-01-10_U2Os_H2030_LR7_TMRE_MG132_Sta9090\RawData';
folderSeg = 'SegmentedImages\';
folderRaw = 'time_0\';
 
% channel = ' wv Blue - FITC - time '; 
%imageloc = [basefolder folderRaw wells1{i1} ' - ' wells2{i2}  '(fld ' num2str(flds(i3),'%02d') channel  num2str(timept(i4),'%02d') ' - ' num2str(timems(i4)) ' ms).tif'];
%For preMG: imageloc = [basefolder folderRaw wells1{i1} ' - ' wells2{i2}  '(fld ' flds{i3} ' time ' num2str(timept(i4),'%01d') ' - ' num2str(timems(i4)) 'ms).tif'];
%For imageloc
% wells1 = {'C';'D';'E';'F'};
% wells2 = {'02','03','04','05','06','07','08','09','10','11'};
wells = {'B - 03', 'B - 04', 'B - 05', 'B - 06', 'B - 07', 'B-08', 'B - 09', 'B - 10', ...
        'C - 03', 'C - 04', 'C - 05', 'C - 06', 'C - 07', 'C - 08', 'C - 09', 'C - 10', ...
        'D - 03', 'D - 04', 'D - 05', 'D - 06', 'D - 07', 'D - 08', 'D - 09', 'D - 10', ...
        'E - 03', 'E - 04', 'E - 05', 'E - 06', 'E - 07', 'E - 08', 'E - 09', 'E - 10', ...
        'F - 07', 'F - 08', 'F - 09', 'F - 10', ...
        'G - 07', 'G - 08', 'G - 09', 'G - 10'}; 
wells2 = {'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09', 'B10', ...
        'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C10', ...
        'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09', 'D10', ...
        'E03', 'E04', 'E05', 'E06', 'E07', 'E08', 'E09', 'E10', ...
        'F07', 'F08', 'F09', 'F10', ...
        'G07', 'G08', 'G09', 'G10'};   
flds = linspace(1,25,25);
channel = {' wv Blue - FITC.tif', ' wv Green dsRed).tif'}; 
channels = {'FITC', 'dsRed'}; 
timept = linspace(1,12,12);
timems = [linspace(0,30,31) 30+linspace(1,4,4)*3]*1200000;% 18+12*3+linspace(1,9,9)*6]*1200000;

cellsize = 50; 

addpath(basefolder)
mkdir([basefolder folderSeg])
tic

%for i2 = 1:length(wells2)
    for i1 = 1:length(wells)
        for i3 = 1:length(flds)
            for i4 = 1:length(channel) 
                toc
           disp([i1 i2 i3 i4])  
            imageloc = [basefolder folderRaw wells{i1}  '(fld '  num2str(flds(i3),'%02d') channel{i4}];
            outputfilename = [basefolder folderSeg wells2{i1} 'fld' num2str(flds(i3),'%02d') '_time00_' channels{i4} '_Seg'];
            Frame_Segmentation_v1(imageloc,outputfilename,cellsize)        

            end
        end
    end
%end

%% images segmentation for timecourse 
clear all

basefolder = 'W:\IN_Cell_Analyzer_6000\Giorgio\2018-12-11 U2OS LR7 cardglyco MG132\RawData\';
folderSeg = 'SegmentedImages\';
folderRaw = 'postMG_timecourse\';
 
% channel = ' wv Blue - FITC - time '; 
%imageloc = [basefolder folderRaw wells1{i1} ' - ' wells2{i2}  '(fld ' num2str(flds(i3),'%02d') channel  num2str(timept(i4),'%02d') ' - ' num2str(timems(i4)) ' ms).tif'];
%For preMG: imageloc = [basefolder folderRaw wells1{i1} ' - ' wells2{i2}  '(fld ' flds{i3} ' time ' num2str(timept(i4),'%01d') ' - ' num2str(timems(i4)) 'ms).tif'];
%For imageloc
wells = {'B - 03', 'B - 04', 'B - 05', 'B - 06', 'B - 07', 'B-08', 'B - 09', 'B - 10', ...
        'C - 03', 'C - 04', 'C - 05', 'C - 06', 'C - 07', 'C - 08', 'C - 09', 'C - 10', ...
        'D - 03', 'D - 04', 'D - 05', 'D - 06', 'D - 07', 'D - 08', 'D - 09', 'D - 10', ...
        'E - 03', 'E - 04', 'E - 05', 'E - 06', 'E - 07', 'E - 08', 'E - 09', 'E - 10', ...
        'F - 07', 'F - 08', 'F - 09', 'F - 10', ...
        'G - 07', 'G - 08', 'G - 09', 'G - 10'}; 
wells2 = {'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09', 'B10', ...
        'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C10', ...
        'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09', 'D10', ...
        'E03', 'E04', 'E05', 'E06', 'E07', 'E08', 'E09', 'E10', ...
        'F07', 'F08', 'F09', 'F10', ...
        'G07', 'G08', 'G09', 'G10'};  
timept = linspace(1,21,21);
timems = [linspace(0,20,21)]*1200000;% 18+12*3+linspace(1,9,9)*6]*1200000;

cellsize = 50; 

addpath(basefolder)
mkdir([basefolder folderRaw folderSeg])
tic

for i2 = 1:length(wells2)
    for i1 = 1:length(wells1)
        for i3 = 1:length(flds)
            for i4 = 1:length(timept)
                toc
           disp([i1 i2 i3 i4])  
            imageloc = [basefolder folderRaw wells{i1}  '(fld '  num2str(flds(i3),'%02d') ' time ' num2str(timept(i4),'%02d') ' - ' num2str(timems(i4)) 'ms).tif'];
            outputfilename = [basefolder folderRaw folderSeg wells1{i1} wells2{i2} 'fld'    flds{i3} '_time' num2str(timept(i4),'%02d') '_Seg'];
            Frame_Segmentation_v1(imageloc,outputfilename,cellsize)        

            end
        end
    end
end


%% drift correction - image registration in time
clear all

basefolder = 'W:\IN_Cell_Analyzer_6000\Giorgio\2018-12-11 U2OS LR7 cardglyco MG132\RawData\postMG_timecourse\';
folderRaw = 'postMG_timecourse\';


posnum  = linspace(1,12,12); %fields 

Well = {'C02','C03','C04','C05','C06','C07','C08','C09','C10','C11',...
        'D02','D03','D04','D05','D06','D07','D08','D09','D10','D11',...
        'E02','E03','E04','E05','E06','E07','E08','E09','E10','E11'...
        'F02','F03','F04','F05','F06','F07','F08','F09','F10','F11'};
for i = 1:36
    
Prefix = [Well{i} 'fld'];
corrfile = [basefolder Well{i} '_drift.mat'];
skip = 0;

time = linspace(1,12,12); %time points s
d_t = [2]; %digits of timept
d_p = [2]; %digits of fields 
shift = 0.15;              % percentage of shift on each side
basesize = [2048 2048];
imagesize = basesize.*[1 1];

save(corrfile,'imagesize','time')

imageregistration_noprint_incell(basefolder,Prefix,posnum,d_p,imagesize,shift,corrfile,time,d_t)

saveregisteredimages_incell(basefolder,Prefix,posnum,d_p,corrfile,time,d_t)

end

%% automated tracking of cells
clear all
basefolder = 'W:\IN_Cell_Analyzer_6000\Giorgio\2018-12-11 U2OS LR7 cardglyco MG132\RawData\postMG_timecourse\';
folderRaw = 'postMG_timecourse\';
Well = {'C02','C03','C04','C05','C06','C07','C08','C09','C10','C11',...
        'D02','D03','D04','D05','D06','D07','D08','D09','D10','D11',...
        'E02','E03','E04','E05','E06','E07','E08','E09','E10','E11'...
        'F02','F03','F04','F05','F06','F07','F08','F09','F10','F11'};

posnum  = linspace(1,25,25);
    
for posnum = 1:12
    posnum
    for i = 1:36
        filename = [basefolder 'AutomaticTracks\' Well{i} 'fld' num2str(posnum,'%02d') '.tiff'];
        if ~exist(filename,'file')

        Prefix = [Well{i} 'fld'];
        corrfile = [basefolder Well{i} '_drift.mat'];
        skip = 0;

        time = linspace(1,12,12);
        d_t = 2;
        d_p = 2;

        automatictracking_skip_incell(basefolder,Prefix,posnum,d_p,corrfile,time,skip)
        end
    end
end
    

%%
clear all

rawdir  = 'Z:\sorger\data\IN Cell Analyzer 6000\Giorgio\2018-03-23 U2OS LR3 live cell TMRE HS RHT\RawImages\timecourse';
segdir  = 'Z:\sorger\data\IN Cell Analyzer 6000\Giorgio\2018-03-23 U2OS LR3 live cell TMRE HS RHT\RawImages\SegmentedImages';
trackdir= 'Z:\sorger\data\IN Cell Analyzer 6000\Giorgio\2018-03-23 U2OS LR3 live cell TMRE HS RHT\RawImages\AutomaticTracks';
results_base = 'Z:\sorger\data\IN Cell Analyzer 6000\Giorgio\2018-03-23 U2OS LR3 live cell TMRE HS RHT\TrackedResults_';

Well = {{'F03','F - 03'},{'F04','F - 04'},{'F05','F - 05'},{'F06','F - 06'},{'F07','F - 07'},...
        {'E03','E - 03'},{'E04','E - 04'},{'E05','E - 05'},{'E06','E - 06'},{'E07','E - 07'},...
        {'D03','D - 03'},{'D04','D - 04'},{'D05','D - 05'},{'D06','D - 06'},{'D07','D - 07'},...
        {'C03','C - 03'},{'C04','C - 04'},{'C05','C - 05'},{'C06','C - 06'},{'C07','C - 07'}};


poscode = linspace(1,22,22);
xy_digits = '%02d';
t_digits = '%02d';
timept = linspace(1,33,33);
realtime = [linspace(0,30,31) 30+linspace(1,4,4)*3]*1200000;
cellsize = 50;
thresh_spatialvar = 10;

for i = 1:20
    
    trackbasepos = [Well{i}{1} 'fld'];
    rawbasepos =  [Well{i}{2} '(fld'];
    results_file = [results_base Well{i}{1} '.mat'];
    save(results_file)
    
    MeasureImages_threecolor_foci_cyto_v5(results_file,rawdir,segdir,trackdir,poscode,rawbasepos,trackbasepos,xy_digits,timept,t_digits,realtime,cellsize,thresh_spatialvar)
end



