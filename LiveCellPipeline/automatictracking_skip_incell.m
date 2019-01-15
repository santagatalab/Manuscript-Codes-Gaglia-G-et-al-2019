function automatictracking_skip_incell(inDir,Prefix,posnum,dp,corrfile,timeNEW,skip)

% posnum  = linspace(66,66,1);
% 
% inDir = '/Volumes/lindquist_lab/linlogin-dropbox/Giorgio/Lindquist Lab/LAB/Microscopy Data/ANALYSIS/2017-03-02 U2OS LR3 HSPA6rep clones test2';
% dataDir = '/Volumes/lindquist_imaging/Giorgio/2017-03-02 U2OS LR3 pHSPA6p-mCer-PEST-NLS clones test2/Tiffs';
% Prefix = '2017-03-2_u2oslr3_hspa6p-mcer_clone_test2';
% corrfile = strcat(inDir,'/',Prefix,'_drift.mat');
% skip = 0;
% 
% time = linspace(2,22,21);
% d_t = [2];
% Channel = [2 3];

p_digits = ['%0' num2str(dp) 'd'];

SegDir = strcat(inDir,'ShiftCorrected\');
TrackDir = strcat(inDir,'AutomaticTracks\');
mkdir(TrackDir);
trackname = strcat(TrackDir,Prefix);

imagebase = strcat(SegDir,Prefix);


load(corrfile);

time = timeNEW;

% imagesize = [1040 1392];
% shift = 0.15;

AreaThresh = 0.2;
SolidThresh = [0.05 0.92];
OverlapThres = 0.5;




for pos = 1:length(posnum)


    
filename = [imagebase num2str(posnum(pos),p_digits) '.tiff'];
image1 = double(imread(filename,1));
 
if skip > 0
    imwrite(zeros(size(image1)), [trackname num2str(posnum(pos),p_digits) '.tiff'] , 'WriteMode','append','Compression','none'); 
end
if skip > 1
    for fillblank = 2:skip
        imwrite(zeros(size(image1)), [trackname num2str(posnum(pos),p_digits) '.tiff'] , 'WriteMode', 'append','Compression','none'); 
    end
end


% clear up cells from the first image that are not the right shape --------
% ie that has solidity less than at threshold specified above -------------

lbimage1 = bwlabel(image1);
stats_t1 = regionprops(lbimage1,'Solidity','PixelIdxList');


for i = 1:max(lbimage1(:))
    if stats_t1(i).Solidity < SolidThresh(2)
        lbimage1(stats_t1(i).PixelIdxList) = 0;
    else
        lbimage1(stats_t1(i).PixelIdxList) = 1;
    end
end

% -------------------------------------------------------------------------

lbimage1 = bwlabel(lbimage1);


image2 = double(imread(filename,2));
image3 = double(imread(filename,3));

lbimage2 = bwlabel(image2);
lbimage3 = bwlabel(image3);

% back calculate the drift

    cumdriftX = cumsum(Drift{1,posnum(pos)}.X);
    cumdriftY = cumsum(Drift{1,posnum(pos)}.Y);
    
    limit_driftX = [-min(cumdriftX) max(cumdriftX)];
    limit_driftY = [-min(cumdriftY) max(cumdriftY)];
    
    row_start = limit_driftY(1)+1;
    row_end   = row_start+imagesize(1)-1;

    col_start = limit_driftX(1)+1;
    col_end   = col_start+imagesize(2)-1;

lbimage1_shiftcorrect = lbimage1(row_start:row_end,col_start:col_end);
lbimage1_shiftcorrect = uint16(lbimage1_shiftcorrect);

imwrite(lbimage1_shiftcorrect, [trackname num2str(posnum(pos),p_digits) '.tiff'] , 'WriteMode','append','Compression','none');

SkippedCells = zeros(size(image1));

% start looping in time ---------------------------------------------------

% figure(1)
% imshow(lbimage1,[])


for tt = 2:(length(time)-2-skip)

stats_image1 = regionprops(lbimage1,'Area','Solidity','Centroid');
stats_image2 = regionprops(lbimage2,'Area','Solidity','Centroid');
stats_image3 = regionprops(lbimage3,'Area','Solidity','Centroid');
 
A_tot = zeros(size(image1));
B_tot = zeros(size(image1));

if max(lbimage1(:)) > 0
    
    stillthere = unique(lbimage1(:)) ;
    stillthere = sort(stillthere);
    
for k1 = 2:length(stillthere)
    
    A = zeros(size(image1));
    B = zeros(size(image1));
    k = stillthere(k1);
if ~isempty(find(lbimage1==k,1))
    [r,c] = find(lbimage1==k);
%     cell1 = bwselect(lbimage1,c,r,8).* lbimage1;
%     cell2 = bwselect(lbimage2,c,r,8).* lbimage2;
    overlap = bwselect(lbimage1,c,r,8).* lbimage2;
    overlap2 = bwselect(lbimage1,c,r,8).* lbimage3;

% first try to get an overlap with the next image by comparing adjacent
% images
    % initialize A

    
    if length(unique(overlap(:))) < 3 && length(unique(overlap(:))) > 1
        AreaDiff = abs( stats_image2(max(overlap(:))).Area/stats_image1(k).Area -1);
        SolidDiff = abs( stats_image2(max(overlap(:))).Solidity/stats_image1(k).Solidity -1);  
        if AreaDiff < AreaThresh && SolidDiff < SolidThresh(1) && ...
                stats_image2(max(overlap(:))).Solidity > SolidThresh(2)
            A = bwselect(lbimage2,stats_image2(max(overlap(:))).Centroid(1),stats_image2(max(overlap(:))).Centroid(2),8)*k; %(:,:,k)
        end
    elseif length(unique(overlap(:))) > 1
        poss_cells = unique(overlap(:),'sorted');
        poss_cells_area = zeros(length(poss_cells),1);
        poss_cells_solidity = zeros(length(poss_cells),1);
        poss_cells_areaoverlap = zeros(length(poss_cells),1);
        for kk = 2:length(unique(overlap(:),'sorted')) 
            poss_cells_area(kk)     =  abs( stats_image2(poss_cells(kk)).Area/stats_image1(k).Area - 1);
            poss_cells_solidity(kk) =  abs( stats_image2(poss_cells(kk)).Solidity/stats_image1(k).Solidity - 1);
            [r2,c2] = find(lbimage2==poss_cells(kk));
            poss_cells_areaoverlap(kk)  = sum(sum(bwselect(lbimage2,c2,r2,8).* bwselect(lbimage1,c,r,8))) /stats_image1(k).Area;
        end
        if max(poss_cells_areaoverlap) > OverlapThres
            [maxoverlap, rightcell] = max(poss_cells_areaoverlap);
        else
            MetricDiff = poss_cells_area + poss_cells_solidity;
            MetricDiff(1) = 100; %fool the code to avodi the first 0
            [minmetric, rightcell] = min(MetricDiff);
        end
        
        A = bwselect(lbimage2,stats_image2(poss_cells(rightcell)).Centroid(1),stats_image2(poss_cells(rightcell)).Centroid(2),8)*k; %(:,:,k)
    end
    
%     if length(A(1,1,:)) < k
%         A(:,:,k) = zeros(size(image1));
%     end
    
% if finding the cell in the next image has failed we try to find it in the
% next frame

    if max(max(A)) == 0 %(:,:,k)
        
        if length(unique(overlap2(:))) < 3 && length(unique(overlap2(:))) > 1
        
            AreaDiff = abs( stats_image3(max(overlap2(:))).Area/stats_image1(k).Area -1);
            SolidDiff = abs( stats_image3(max(overlap2(:))).Solidity/stats_image1(k).Solidity -1);  
            if AreaDiff < AreaThresh && SolidDiff < SolidThresh(1) && ...
                   stats_image3(max(overlap2(:))).Solidity > SolidThresh(2)
                B = bwselect(lbimage3,stats_image3(max(overlap2(:))).Centroid(1),stats_image3(max(overlap2(:))).Centroid(2),8)*k; %(:,:,k)
            end
        elseif length(unique(overlap2(:))) > 1
            poss_cells = unique(overlap2(:),'sorted');
            poss_cells_area = zeros(length(poss_cells),1);
            poss_cells_solidity = zeros(length(poss_cells),1);
            poss_cells_areaoverlap = zeros(length(poss_cells),1);
            for kk = 2:length(unique(overlap2(:),'sorted')) 
                poss_cells_area(kk)     =  abs( stats_image3(poss_cells(kk)).Area/stats_image1(k).Area - 1);
                poss_cells_solidity(kk) =  abs( stats_image3(poss_cells(kk)).Solidity/stats_image1(k).Solidity - 1);
                [r2,c2] = find(lbimage3==poss_cells(kk));
                poss_cells_areaoverlap(kk)  = sum(sum(bwselect(lbimage3,c2,r2,8).* bwselect(lbimage1,c,r,8))) /stats_image1(k).Area;
            end
            if max(poss_cells_areaoverlap) > OverlapThres
            [maxoverlap, rightcell] = max(poss_cells_areaoverlap);
            else
            MetricDiff = poss_cells_area + poss_cells_solidity;
            MetricDiff(1) = 100; %fool the code to avodi the first 0
            [minmetric, rightcell] = min(MetricDiff);
            end
        
            B = bwselect(lbimage3,stats_image3(poss_cells(rightcell)).Centroid(1),stats_image3(poss_cells(rightcell)).Centroid(2),8)*k; %(:,:,k)
        
        end        
        
        
    end
    
end



A_tot = A_tot + A;
B_tot = B_tot + B;

% figure(1)
% imshow(A_tot,[])

end
end

TrackedCells = A_tot + SkippedCells;
SkippedCells = B_tot;
clear A_tot
clear B_tot


row_start = limit_driftY(1)+cumdriftY(time(tt))+1;
row_end   = row_start+imagesize(1)-1;

col_start = limit_driftX(1)+cumdriftX(time(tt))+1;
col_end   = col_start+imagesize(2)-1;

TrackedCells_shiftcorrect = TrackedCells(row_start:row_end,col_start:col_end);
TrackedCells_shiftcorrect = uint16(TrackedCells_shiftcorrect);
% imwrite(TrackedCells(:,:,pos)/255, strcat(outbasename,num2str(position,'%03d'),'.tiff') , 'WriteMode', 'append','Compression','none');
imwrite(TrackedCells_shiftcorrect, [trackname num2str(posnum(pos),p_digits) '.tiff'] , 'WriteMode', 'append','Compression','none');

lbimage1 = TrackedCells;

image2 = double(imread(filename,tt+1));
image3 = double(imread(filename,tt+2));
lbimage2 = bwlabel(image2);

lbimage3 = bwlabel(image3);

% figure(tt)
% imshow(lbimage1,[])

[posnum(pos),tt]

end

end








