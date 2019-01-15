function MeasureImages_threecolor_foci_cyto_v5(results_file,rawdir,segdir,trackdir,poscode,rawbasepos,trackbasepos,xy_digits,timept,t_digits,realtime,cellsize,thresh_spatialvar)

results_key = ['1 - CellAnnotation Time-------';...
               '2 - CellAnnotation Index------';...
               '3 - centroid row y pos--------';...
               '4 - centroid col x pos--------';... 
               '5 - Nuclear Segmented Solidity';... 
               '6 - Nuclear Segmented Area----';...
               '7 - Nuclear FITC-488 Fluor----';...
               '8 - Nuclear dsRed-555 Fluor---';...
               '9 - Nuclear Cy5-647 Fluor-----';...
               '10 - Cytoplasm Segmented Area-';...
               '11 - Cytoplasm FITC-488 Fluor-';...
               '12 - Cytoplasm dsRed-555 Fluor';...
               '13 - Cytoplasm Cy5-647 Fluor--';...
               '14 - Nuclear FITC-488 Total---';...
               '15 - Nuclear FITC-488 Foci ---'...
               ];
% Include crispness??? 
    %if max(Foci_Image(:)) > 0
       % lb_Borders = Borders.*(lb_Nuc_Image+lb_CytImage);
       % Borders_stats = regionprops(lb_Borders,SpatVar_Raw,'PixelValues');
        % Crispness of cells
       % OutputMatrix(cell_start:cell_end,5) = cellfun(@mean,{Borders_stats.PixelValues});
   % end

save(results_file,'-append')
tic
for i = 1:length(poscode)

tracks  = [trackdir '\' trackbasepos num2str(poscode(i),xy_digits) '.tiff'];
Track_1 = uint8(imread(tracks,1));
totcells = max(Track_1(:));
    
    for t = 1:length(timept)
        Results.Field{poscode(i)}.TimePoint{t}.Matrix = zeros(totcells,15);
        
%%%%%% THIS IS WHERE THE ACTION BEGINS! %%%%%%%%%%%%   
        
        % load the channels you are interested in and the tracking data
        fileA488    = [rawdir '\' rawbasepos ' ' num2str(poscode(i),xy_digits) ' wv Blue - FITC - time ' num2str(timept(t),t_digits) ' - ' num2str(realtime(t)) ' ms).tif'];
        fileA555    = [rawdir '\' rawbasepos ' ' num2str(poscode(i),xy_digits) ' wv Green - dsRed - time ' num2str(timept(t),t_digits) ' - ' num2str(realtime(t)) ' ms).tif'];
        fileA647    = [rawdir '\' rawbasepos ' ' num2str(poscode(i),xy_digits) ' wv Red - Cy5 - time ' num2str(timept(t),t_digits) ' - ' num2str(realtime(t)) ' ms).tif'];
        
        Track_t = uint16(imread(tracks,t));
        
        % load images and do BK subtraction as soon as the come in
        A488_Image = uint16(imread(fileA488));
        A488_BK = A488_Image-imopen(A488_Image, strel('disk',cellsize*2));
        % ADD IN BASIC CORRECTION!?!?!?!?
        flagA555 = 0;
        flagA647 = 0;
        
        if exist(fileA555,'file')
            flagA555 = 1;
            A555_Image = uint16(imread(fileA555));
            A555_BK = A555_Image-imopen(A555_Image, strel('disk',cellsize*2)); 
        end
        if exist(fileA647,'file')
            flagA647 = 1;
            A647_Image = uint16(imread(fileA647));
            A647_BK = A647_Image-imopen(A647_Image, strel('disk',cellsize*2)); 
        end
        
        % we already have the labelled nuclear mask
        lb_Nuc_Image = Track_t;
        lb_Nuc_Image(lb_Nuc_Image(:)>totcells) = 0;
        SegImage = lb_Nuc_Image > 0;
        
        % create cytoplasmic mask
        lb_CytImage = imdilate(lb_Nuc_Image,offsetstrel('ball',20,0));
        lb_CytImage(lb_Nuc_Image>0)=0;
        CytImage = lb_CytImage > 0;
        imwrite(uint16(CytImage),[segdir filesep trackbasepos num2str(poscode(i),xy_digits) '_time' num2str(timept(t),t_digits) '_CytoSeg.tif'])
        
        stats_NucImage = regionprops(lb_Nuc_Image,'Area','Solidity','Centroid');
        stats_CytImage = regionprops(lb_CytImage,'Area');

        % segment the foci and non foci parts of the nuclei
        % segement HSF1         
        HSF1Foci_BK = A488_Image-imopen(A488_Image, strel('disk',round(cellsize/3))); 
        SpatVar_Seg = stdfilt(stdfilt(Track_t));
        SpatVar_Raw = stdfilt(HSF1Foci_BK)./double(sqrt(double(A488_Image)));
        % FUNDAMENTAL PARAMETER HERE!!!!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     figure(11)
        %     imshow(A488_Image,[0 20000])
    
        FociSeg = SpatVar_Raw>thresh_spatialvar;       
        %     figure(1)
        %     imshow(SpatVar_Raw,[])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        FociSeg = imfill(FociSeg,'holes');
        %     figure(2)
        %     imshow(FociSeg,[])
        FociSeg = imclose(FociSeg,strel('disk',1));
        %     figure(2)
        %     imshow(FociSeg,[])
        FociSeg = imfill(FociSeg,'holes');
        %     figure(3)
        %     imshow(FociSeg,[])
        FociSeg = bwmorph(FociSeg,'thin',1);
        %     figure(4)
        %     imshow(FociSeg,[])
        %     FociSeg = imfill(FociSeg,'holes');
        % make binary foci mask
        FociSeg = uint16(FociSeg).*uint16(SegImage) > 0;
        %     figure(6)
        %     imshow(FociSeg,[])
    
        imwrite(uint16(FociSeg),[segdir filesep trackbasepos num2str(poscode(i),xy_digits) '_time' num2str(timept(t),t_digits) '_FociSeg.tif'])
            
        lb_FociImage = uint16(FociSeg).*lb_Nuc_Image;
        Foci_stats = regionprops(lb_FociImage,A488_BK,'PixelValues');
        
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_FociImage(:))),15) = cellfun(@sum,{Foci_stats.PixelValues});
        
        Area = {stats_NucImage.Area};
        Solidity = {stats_NucImage.Solidity};
        Centroid = {stats_NucImage.Centroid};
        CentroidVect = cell2mat(Centroid);
        CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
        Area_Cyt = {stats_CytImage.Area};
        Nuclei_stats = regionprops(lb_Nuc_Image,lb_Nuc_Image,'PixelValues');
        
        
        % totcells is undefined - need to make sure only to capture cells
        % present from the first timepoint
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(:,1) = zeros(totcells,1)+realtime(t);    
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),2) = cellfun(@mean,{Nuclei_stats.PixelValues});
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),3) = CentroidMat(1,:);
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),4) = CentroidMat(2,:); 
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),5) = cell2mat(Solidity);
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),6) = cell2mat(Area);
        
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),10)= cell2mat(Area_Cyt);
        
        Nuclei_A488_stats = regionprops(lb_Nuc_Image,A488_BK,'PixelValues');
        Cytopl_A488_stats = regionprops(lb_CytImage,A488_BK,'PixelValues');

        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),7) = cellfun(@mean,{Nuclei_A488_stats.PixelValues});
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),11) = cellfun(@mean,{Cytopl_A488_stats.PixelValues});
        Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),14) = cellfun(@sum,{Nuclei_A488_stats.PixelValues});
        
        if flagA555 > 0
            Nuclei_A555_stats = regionprops(lb_Nuc_Image,A555_BK,'PixelValues');
            Cytopl_A555_stats = regionprops(lb_CytImage,A555_BK,'PixelValues');
            Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),8) = cellfun(@mean,{Nuclei_A555_stats.PixelValues});
            Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),12) = cellfun(@mean,{Cytopl_A555_stats.PixelValues});
        end
        if flagA647 > 0
            Nuclei_A647_stats = regionprops(lb_Nuc_Image,A647_BK,'PixelValues');
            Cytopl_A647_stats = regionprops(lb_CytImage,A647_BK,'PixelValues');
            Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),9) = cellfun(@mean,{Nuclei_A647_stats.PixelValues});
            Results.Field{poscode(i)}.TimePoint{t}.Matrix(1:(max(lb_Nuc_Image(:))),13) = cellfun(@mean,{Cytopl_A647_stats.PixelValues});
        end
        
        
        
        
    disp([num2str(poscode(i)) ' ' num2str(t)])
    save(results_file,'Results','-append')
    end
    toc
    
end



