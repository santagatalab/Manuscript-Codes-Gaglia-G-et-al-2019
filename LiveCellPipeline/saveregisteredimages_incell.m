function saveregisteredimages_incell(inDir,Prefix,posnum,dp,corrfile,time_vect,dt)

ShiftDir = strcat(inDir,'ShiftCorrected\');
mkdir(ShiftDir)

segm_shift = strcat(ShiftDir,Prefix);
t_digits = strcat('%0',num2str(dt),'d');
p_digits = strcat('%0',num2str(dp),'d');
load(corrfile);



for pos = 1:length(posnum)

    % ---- use the correlation data to create stacks -----
    
    % calculate the running sum of the correlation
    Drift{1,posnum(pos)}.X(1:1) = 0;  %set the skipped corrs to 0
    Drift{1,posnum(pos)}.Y(1:1) = 0;
    cumdriftX = cumsum(Drift{1,posnum(pos)}.X);
    cumdriftY = cumsum(Drift{1,posnum(pos)}.Y);

    limit_driftX = [-min(cumdriftX) max(cumdriftX)];
    limit_driftY = [-min(cumdriftY) max(cumdriftY)];

    % create padding images for this positiong
    
    padding_img = zeros(imagesize(1)+sum(limit_driftY),...
                        imagesize(2)+sum(limit_driftX));
    
    segm_name1 = [inDir 'SegmentedImages\' Prefix num2str(posnum(pos),p_digits) '_time' num2str(time_vect(1),t_digits) '_Seg.tif'];
    
    segm_img1 = imread(segm_name1);
    
    row_start = limit_driftY(1)+1;
    row_end   = row_start+imagesize(1)-1;

    col_start = limit_driftX(1)+1;
    col_end   = col_start+imagesize(2)-1;

    shift_segm_img1 = padding_img;
    shift_segm_img1(row_start:row_end,col_start:col_end) = segm_img1;
   
    imwrite(shift_segm_img1, [segm_shift num2str(posnum(pos),p_digits) '.tiff'] , 'WriteMode', 'append','Compression','none');

    
    for tt = 2:length(time_vect)
        
        segm_name_t = [inDir 'SegmentedImages\' Prefix num2str(posnum(pos),p_digits) '_time' num2str(time_vect(tt),t_digits) '_Seg.tif'];
        segm_img_t = imread(segm_name_t);
 
        row_start = limit_driftY(1)+cumdriftY(time_vect(tt))+1;
        row_end   = row_start+imagesize(1)-1;

        col_start = limit_driftX(1)+cumdriftX(time_vect(tt))+1;
        col_end   = col_start+imagesize(2)-1;
    
    shift_segm_img = padding_img;

    shift_segm_img(row_start:row_end,col_start:col_end) = segm_img_t;

    imwrite(shift_segm_img, strcat(segm_shift,num2str(posnum(pos),'%02d'),'.tiff') , 'WriteMode', 'append','Compression','none');
    
       [posnum(pos), time_vect(tt)]
        
    end

end


