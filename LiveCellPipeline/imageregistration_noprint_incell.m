function imageregistration_noprint_incell(inDir,Prefix,posnum,dp,imagesize,shift,corrfile,time_vect,dt)

% OutDir = strcat(inDir,'/ShiftCorrected');
% mkdir(OutDir);
% OutDir2 = strcat(inDir,'/SegmentedImagesShiftCorrected');
% mkdir(OutDir2);

t_digits = strcat('%0',num2str(dt),'d');
p_digits = strcat('%0',num2str(dp),'d');

load(corrfile);

for pos = 1:length(posnum)

imagebase = strcat(inDir,'SegmentedImages\',Prefix);

% define the size of the smaller image to be correlated using the parameter ?shift?
Xl = [round(imagesize(2)*shift) round(imagesize(2)*(1-shift))]
Yl = [round(imagesize(1)*shift) round(imagesize(1)*(1-shift))]

% create a matrix of correlations to keep track of the shifting images
corr_offset = zeros(1,2);

 % start with loading the first two timepoints of the stack we are going to align
filename1 = [imagebase num2str(posnum(pos),p_digits) '_time' num2str(time_vect(1),t_digits) '_Seg.tif'];
filename2 = [imagebase num2str(posnum(pos),p_digits) '_time',num2str(time_vect(2),t_digits),'_Seg.tif'];
 

% now loop starting from the second position onwards, recording the offset and adjusting the position of the image

for tt = 2:length(time_vect)
    
    referimg = imread(filename1);
    sampleimg = imread(filename2);
% %     eval(['saveimg = imread(filenameout);'])
    
    if sum(sampleimg(:)) == 0
        corr_offset(time_vect(tt),1) = 0;
        corr_offset(time_vect(tt),2) = 0;
    else
        c = normxcorr2(sampleimg(Yl(1):Yl(2),Xl(1): Xl(2)),referimg);
                 
        [max_c, imax] = max(c(:));
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
    
        limit_test(1) = abs(xpeak-Xl(2)) - imagesize(1)*shift;
        limit_test(2) = abs(xpeak-Xl(2)) - imagesize(2)*shift;
    
        if max_c > 0.2 && max(limit_test) < 0
            corr_offset(time_vect(tt),1) = xpeak - Xl(2);
            corr_offset(time_vect(tt),2) = ypeak - Yl(2);
        else
            corr_offset(time_vect(tt),1) = 0;
            corr_offset(time_vect(tt),2) = 0;
        end
    end

    filename1 = [imagebase num2str(posnum(pos),p_digits) '_time' num2str(time_vect(tt),t_digits) '_Seg.tif'];
    filename2 = [imagebase num2str(posnum(pos),p_digits) '_time',num2str(time_vect(tt)+1,t_digits),'_Seg.tif'];

clear limit_test

[posnum(pos),time_vect(tt)]
end
posnum(pos)
Drift{posnum(pos)}.X = corr_offset(:,1);
Drift{posnum(pos)}.Y = corr_offset(:,2);
size_shift(posnum(pos))   = shift;

% save drifting parameters for tracking
save(corrfile,'Drift','size_shift','-append');
end





