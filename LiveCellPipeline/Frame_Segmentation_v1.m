function Frame_Segmentation_v1(input,output,cellsize)

% cellsize = 25;
% input = 'timecourse\B - 02(fld 01 wv Blue - FITC - time 01 - 0 ms).tif';
% MagPar = 1;
            %       *IMAGE PROCESSING GOES HERE.*
            
% --------------------------- FINDING NUCLEI ------------------------------
% takes the DNA image and gives back a "Nuclei" image

RawImage = double(imread(input));
coImage=imopen(imclose(RawImage,strel('disk',cellsize/5)),strel('disk',cellsize));
sImage=RawImage-coImage;%imopen(imclose(RawImage-coImage,strel('disk',3)),strel('disk',15*MagPar));

% figure(1)
% imshow(RawImage,[])
% figure(2)
% imshow(coImage,[])
% figure(3)
% imshow(sImage,[])
% figure(4)
% ksdensity(log2(sImage(sImage>8)))

% new way of determining the threshold
% old_thr = prctile(sImage(~imclose(sImage<0,strel('disk',5))),95)/8

% ThresImage = imopen(RawImage-coImage,strel('disk',10));
ThresImage = imopen(RawImage-coImage,strel('disk',cellsize/5));
% figure(5)
% imshow(ThresImage,[])
% thImage = stdfilt(ThresImage)./ThresImage>0.3;
% thImage = thImage.*ThresImage;
% figure(6)
% imshow(thImage,[0 1500])
ThreshImage2 = sImage;
% ThreshImage2(thImage<1) = NaN;

temp_v = ThreshImage2(ThreshImage2>0);
thr = prctile(temp_v,50);
clear temp_v

seeds=imregionalmax(imfilter(sImage,fspecial('gaussian',round(cellsize*0.9),round(cellsize/4))));
seeds=seeds&sImage>thr;
seeds=imdilate(seeds,strel('disk',4));

% figure(1)
% imshow(ThresImage+seeds*2000,[])

ThresImage=sImage>thr;

L=watershed(imimposemin(-sImage,seeds));

Nuclei=ThresImage&L;
Nuclei=bwareaopen(Nuclei,cellsize*2);


% FIXING BROKEN CELLS
% LineSegments=ThresImage&~L;
% LineSegments=bwareaopen(LineSegments,cellsize/10);
% LineSegments=bwlabel(LineSegments,8);
% % imshow(LineSegments,[])
% 
% 
% 
% for Seg_i=1:max(LineSegments(:))
%     
%     FullCell=imreconstruct(LineSegments==Seg_i,ThresImage);
% %     imshow(ThresImage,[])
%     BrokenCells=bwlabel(FullCell&~LineSegments>0);
%     BrokenStats=regionprops(BrokenCells,'solidity','area');
%     FullStats=regionprops(FullCell,'solidity','area');
%     
%     for Seg_k=1:max(BrokenCells(:))
%         BrokenArea(Seg_k)=BrokenStats(Seg_k).Area;
%         BrokenSolidity(Seg_k)=BrokenStats(Seg_k).Solidity;
%     end
%     %
%     FullArea=FullStats(1).Area;
%     FullSolidity=FullStats(1).Solidity;
%     if sum(FullSolidity>BrokenSolidity) && FullArea<10000
%         
%         Nuclei=Nuclei|LineSegments==Seg_i;
%     end
%     
%     Nuclei=imopen(Nuclei,strel('disk',2));
% end
 
Nuclei = imfill(Nuclei,'holes');
SegmentedImage = bwlabel(Nuclei,8);
 
imwrite(uint16(SegmentedImage),[output '.tif'])

% overlay = cat(3,uint16(sImage),uint16(Nuclei));
% overlay = cat(3,overlay,(zeros(size(RawImage))));
% % figure(18)
% % imshow(overlay,[])
% 
% imwrite(overlay,[output '_overlay.tif'])

% figure(5)
% imshow(double(SegmentedImage),[])
end 
