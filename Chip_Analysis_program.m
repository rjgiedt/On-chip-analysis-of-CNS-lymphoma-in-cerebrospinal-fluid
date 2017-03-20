%% Initial setup
clear all; close all;

%% CORRECTION FACTORS: ENTER FILTERS HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Filter_disk = 2; % Probably don't screw with this... just removes small dots
minimum_size_cutoff = 100; % Minimum size for structure to be included
maximum_size_cutoff = 2500; % Maximum size for structure to be included
Dilation_disk = 12; %Controls How much of the area around regions you see... the bigger you pick the more size and circular you get!
nth_percent = 0.25; % Top nth% of pixels for each region controlled

%% Input Video (TRITC)
input_video_fname = 'Input.avi';
              
%Prompt for input video file name (use above name as default)
[input_video_fname,input_video_path] = uigetfile('*.avi',...
    'Select the TRITC input video',input_video_fname);
if isequal(input_video_fname,0) || isequal(input_video_path,0)
   disp('User selected Cancel')
   return;
else
   fprintf('User selected "%s" as input video\n',input_video_fname);
end

%% Input Video (DAPI)
input_video_fname2 = 'Input2.avi';
              
%Prompt for input video file name (use above name as default)
[input_video_fname2,input_video_path2] = uigetfile('*.avi',...
    'Select the DAPI video',input_video_fname2);
if isequal(input_video_fname2,0) || isequal(input_video_path2,0)
   disp('User selected Cancel')
   return;
else
   fprintf('User selected "%s" as input video\n',input_video_fname2);
end

%%  Input Video (CY5)
input_video_fname3 = 'Input3.avi';
              
%Prompt for input video file name (use above name as default)
[input_video_fname3,input_video_path3] = uigetfile('*.avi',...
    'Select the CY5 video',input_video_fname3);
if isequal(input_video_fname3,0) || isequal(input_video_path3,0)
   disp('User selected Cancel')
   return;
else
   fprintf('User selected "%s" as input video\n',input_video_fname3);
end

%% Read in selected videos 
% Read in information for TRITC Video
movie=VideoReader(input_video_fname);  %Read in TRITC video
nframes = movie.NumberOfFrames;  %Determine the total number of frames

% Read in information of DAPI Video
movie2=VideoReader(input_video_fname2); %Read in Dapi video
nframes2 = movie2.NumberOfFrames; %Determine total number of frames for Dapi video

% Read in information CY5 Video
movie3 = VideoReader(input_video_fname3); % Read in CY5 Video
nframes3 = movie3.NumberOfFrames; % Determine total number of frames for CY5 video 


%% Loop to Find information on individual Cells
% Close previous windows
close all;

parfor n=1:nframes
n
    % Segment frame from TRITC Channel
I = read(movie, n); % Read in movie image  
dapi = read(movie2,n); % Read in dapi image
cy5 = read(movie3, n); % Read in Cy5 image

%[t, threshold] = huangentropy2(I); % Segment via Huang Entropy Method
level = graythresh(I);
if level>=0.100000001
level_correct(n) = level-0.1;
else
    level_correct(n) = 0.1;
end

threshold_image = im2bw(I,level_correct(n));

% Filter Results from Huang to eliminate potential noise
% Standard disk filtering
se = strel('disk',Filter_disk); % Set disk size
B = threshold_image;
Io = imopen(B,se);
Ie = imerode(B, se);
Iobr = imreconstruct(Ie, B);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
Iobrcbr = double(Iobrcbr);

Iobrcbr_filled = imfill(Iobrcbr, 'holes');
%Iobrcbr_filled = rgb2gray(Iobrcbr_filled);
Iobrcbr_filled_reverse = Iobrcbr_filled;
%Iobrcbr_filled_reverse = imadjust(Iobrcbr_filled,[0 1], [1 0]);

% Filter Objects based on sizes
min_cutoff = minimum_size_cutoff; % Original 75... prof 45
max_cutoff = maximum_size_cutoff; %Original 1000...prof 1500
size_filtering = bwlabel(Iobrcbr_filled_reverse);
s = regionprops(size_filtering, 'All');

% Filter based on total area sizes
area_values = [s.Area];
idx = find((min_cutoff <= area_values) & (max_cutoff >= area_values));
size_filtered_images = ismember(size_filtering, idx);

% Find properties of filtered image
segment = bwlabel(size_filtered_images); 
properties_size_regions = regionprops(segment, 'All');

% Dilate for each centroid found in segmention
large_disk = strel('disk', Dilation_disk);
dilate_image = imdilate(size_filtered_images, large_disk);
dilate_image_label = bwlabel(dilate_image);
 
% Calculate Background Values
% Create Background Masks
background_mask_segment = imadjust(im2double(size_filtered_images),...
    [0 1], [1 0]);  
background_mask_dilate = imadjust(im2double(dilate_image),...
    [0 1], [1 0]);

% Find Intensity of Background for each mask with each channel
backgroundprops_segment_dapi = regionprops(logical(background_mask_segment),...
    rgb2gray(dapi),'MeanIntensity');
backgroundprops_segment_cy5 = regionprops(logical(background_mask_segment),...
    rgb2gray(cy5),'MeanIntensity');
backgroundprops_dilate_dapi = regionprops(logical(background_mask_dilate),...
    rgb2gray(dapi),'MeanIntensity');
backgroundprops_dilate_cy5 = regionprops(logical(background_mask_dilate),...
    rgb2gray(cy5),'MeanIntensity'); 

% Take largest area mean intensity as background signal for each image
background_segment_dapi(n) = backgroundprops_segment_dapi(1).MeanIntensity(1,1);
background_segment_cy5(n) = backgroundprops_segment_cy5(1).MeanIntensity(1,1);
background_dilate_dapi(n) = backgroundprops_dilate_dapi(1).MeanIntensity(1,1);
background_dilate_cy5(n) = backgroundprops_dilate_cy5(1).MeanIntensity(1,1);

%% Find properties of images for each segmented cell with NO BACKGROUND SUBTRACTION
% Properties for Dapi
regionprops_segment_dapi{n} = regionprops(segment, rgb2gray(dapi),...
    'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Area');
regionprops_dilate_dapi{n} = regionprops(dilate_image_label, rgb2gray(dapi),...
    'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Area');

% Properties for Cy5
regionprops_segment_cy5{n} = regionprops(segment, rgb2gray(cy5),...
    'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Area');
regionprops_dilate_cy5{n} = regionprops(dilate_image_label, rgb2gray(cy5),...
    'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Area');

% Find Intensity of each individual pixel in each region
% Properties for Dapi
regionprops_pixel_segment_dapi{n} = regionprops(segment, rgb2gray(dapi),...
    'PixelValues');
regionprops_pixel_dilate_dapi{n} = regionprops(dilate_image_label, rgb2gray(dapi),...
    'PixelValues');

% Properties for Cy5
regionprops_pixel_segment_cy5{n} = regionprops(segment, rgb2gray(cy5),...
    'PixelValues');
regionprops_pixel_dilate_cy5{n} = regionprops(dilate_image_label, rgb2gray(cy5),...
    'PixelValues');


%% FIND PROPERTIES OF IMAGES WITH BACKGROUND SUBTRACTION
regionprops_segment_dapi_b{n} = regionprops(segment, rgb2gray(dapi)...
    -background_segment_dapi(n),...
    'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Area');
regionprops_dilate_dapi_b{n} = regionprops(dilate_image_label, rgb2gray(dapi)...
    -background_dilate_dapi(n),...
    'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Area');

% Properties for Cy5
regionprops_segment_cy5_b{n} = regionprops(segment, rgb2gray(cy5)...
    -background_segment_cy5(n),...
    'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Area');
regionprops_dilate_cy5_b{n} = regionprops(dilate_image_label, rgb2gray(cy5)...
    -background_dilate_cy5(n),...
    'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Area');

% Find Individual Pixel Values for each region identified
regionprops_pixel_segment_dapi_b{n} = regionprops(segment, rgb2gray(dapi)...
    -background_segment_dapi(n),'PixelValues');
regionprops_pixel_dilate_dapi_b{n} = regionprops(dilate_image_label, rgb2gray(dapi)...
    -background_dilate_dapi(n),'PixelValues');

% Properties for Cy5
regionprops_pixel_segment_cy5_b{n} = regionprops(segment, rgb2gray(cy5)...
    -background_segment_cy5(n),'PixelValues');
regionprops_pixel_dilate_cy5_b{n} = regionprops(dilate_image_label, rgb2gray(cy5)...
    -background_dilate_cy5(n),'PixelValues');

% Save Images for compatibility with parallel processing for loop
segmented_final{n} = segment;
dilated_final{n} = dilate_image_label;

% Create overlay output images
% Create overlay of borders
overlay_segmented{n} = imoverlay(I, bwperim(segment), [.3 1 .3]);
overlay_dilated{n} = imoverlay(I, bwperim(dilate_image_label), [.3 1 .3]);

overlay_segmented_dapi{n} = imoverlay(dapi, bwperim(segment), [.3 1 .3]);
overlay_segmented_cy5{n} = imoverlay(cy5, bwperim(segment), [.3 1 .3]);
overlay_dilated_dapi{n} = imoverlay(dapi, bwperim(dilate_image_label), [.3 1 .3]);
overlay_dilated_cy5{n} = imoverlay(cy5, bwperim(dilate_image_label), [.3 1 .3]);

end

%% Scatter Plot Analysis
% Create variables for all analyzed regions with NO BACKGROUND SUBTRACTION
% Create Variables
conc_segment_dapi = regionprops_segment_dapi{1,1};
conc_dilate_dapi = regionprops_dilate_dapi{1,1};
conc_segment_cy5 = regionprops_segment_cy5{1,1};
conc_dilate_cy5 = regionprops_dilate_cy5{1,1};

for n = 2:nframes
conc_segment_dapi = [conc_segment_dapi; regionprops_segment_dapi{1,n}];
conc_dilate_dapi = [conc_dilate_dapi; regionprops_dilate_dapi{1,n}];
conc_segment_cy5 = [conc_segment_cy5; regionprops_segment_cy5{1,n}];
conc_dilate_cy5 = [conc_dilate_cy5; regionprops_dilate_cy5{1,n}];
end

% Create variables for all analyzed regions with BACKGROUND SUBTRACTION
% Create Variables
conc_segment_dapi_b = regionprops_segment_dapi_b{1,1};
conc_dilate_dapi_b = regionprops_dilate_dapi_b{1,1};
conc_segment_cy5_b = regionprops_segment_cy5_b{1,1};
conc_dilate_cy5_b = regionprops_dilate_cy5_b{1,1};

for n = 2:nframes
conc_segment_dapi_b = [conc_segment_dapi_b; regionprops_segment_dapi_b{1,n}];
conc_dilate_dapi_b = [conc_dilate_dapi_b; regionprops_dilate_dapi_b{1,n}];
conc_segment_cy5_b = [conc_segment_cy5_b; regionprops_segment_cy5_b{1,n}];
conc_dilate_cy5_b = [conc_dilate_cy5_b; regionprops_dilate_cy5_b{1,n}];
end

%% Scatter Plots
% Create Variables outside of structure loops
% Segmented Variables
segment_dapi_av_p = struct2cell(conc_segment_dapi);
segment_cy5_av_p = struct2cell(conc_segment_cy5);
segment_dapi_av_b_p = struct2cell(conc_segment_dapi_b);
segment_cy5_av_b_p = struct2cell(conc_segment_cy5_b);

% Dilated Variables 
dilate_dapi_av_p = struct2cell(conc_dilate_dapi);
dilate_cy5_av_p = struct2cell(conc_dilate_cy5);
dilate_dapi_av_b_p = struct2cell(conc_dilate_dapi_b);
dilate_cy5_av_b_p = struct2cell(conc_dilate_cy5_b);

% Create single variables
% Segmented Variables
segement_dapi_av = segment_dapi_av_p(2,:);
segement_cy5_av = segment_cy5_av_p(2,:);
segment_dapi_av_b = segment_dapi_av_b_p(2,:);
segment_cy5_av_b = segment_cy5_av_b_p(2,:);

% Dilated Variables 
dilate_dapi_av = dilate_dapi_av_p(2,:);
dilate_cy5_av = dilate_cy5_av_p(2,:);
dilate_dapi_av_b = dilate_dapi_av_b_p(2,:);
dilate_cy5_av_b = dilate_cy5_av_b_p(2,:);

% Segmented Graphs
figure; scatter(cell2mat(segement_dapi_av), cell2mat(segement_cy5_av));
title('Plot of segmented data with no background reduction')
figure; scatter(cell2mat(segment_dapi_av_b), cell2mat(segment_cy5_av_b));
title('Plot of segmented data WITH background reduction')

% Dilated Graphs
figure; scatter(cell2mat(dilate_dapi_av), cell2mat(dilate_cy5_av));
title('Plot of Dilated data with no background reduction')
figure; scatter(cell2mat(dilate_dapi_av_b), cell2mat(dilate_cy5_av_b));
title('Plot of Dilated data WITH background reduction')


%% Create Output Videos
% Write Segment Video
writerObj = VideoWriter('segment_video.avi');
writerObj.FrameRate = 1;  
open(writerObj);

for n=1:nframes   
 inputframe = im2double(overlay_segmented{n});
 writeVideo(writerObj,inputframe); 
end

close(writerObj);

% Write Dilate Video
writerObj2 = VideoWriter('dilate_video.avi');
writerObj2.FrameRate = 1;  
open(writerObj2);

for n=1:nframes   
 inputframe = im2double(overlay_dilated{n});
 writeVideo(writerObj2,inputframe); 
end

close(writerObj2);

% Write Segment Video for Dapi
writerObj3 = VideoWriter('segmented_dapi_video.avi');
writerObj3.FrameRate = 1;  
open(writerObj3);

for n=1:nframes   
 inputframe = im2double(overlay_segmented_dapi{n});
 writeVideo(writerObj3,inputframe); 
end

close(writerObj3);

% Write Segment Video for Cy5
writerObj4 = VideoWriter('segmented_cy5_video.avi');
writerObj4.FrameRate = 1;  
open(writerObj4);

for n=1:nframes   
 inputframe = im2double(overlay_segmented_cy5{n});
 writeVideo(writerObj4,inputframe); 
end

close(writerObj4);

% Write Dilate Video for Dapi
writerObj5 = VideoWriter('dilated_dapi_video.avi');
writerObj5.FrameRate = 1;  
open(writerObj5);

for n=1:nframes   
 inputframe = im2double(overlay_dilated_dapi{n});
 writeVideo(writerObj5,inputframe); 
end

close(writerObj5);

% Write Dilate Video for Cy5
writerObj6 = VideoWriter('dilated_cy5_video.avi');
writerObj6.FrameRate = 1;  
open(writerObj6);

for n=1:nframes   
 inputframe = im2double(overlay_dilated_cy5{n});
 writeVideo(writerObj6,inputframe); 
end

close(writerObj6);

%% Additional Analysis of pixel regions 
% For segmented, no background subtraction
for n=1:nframes % For each FRAME of video
    for j = 1:length(regionprops_pixel_segment_dapi{1,n}) % For Each Region Found in the frame
    sort_segment_dapi{n,j} = sort(regionprops_pixel_segment_dapi{1,n}(j).PixelValues,'descend');  
    sort_segment_cy5{n,j} = sort(regionprops_pixel_segment_cy5{1,n}(j).PixelValues, 'descend');
    % Find Average of predefined intensities
    t = sort_segment_dapi{n,j};
    sintensity_segment_dapi{n,j} = mean(t(1:floor(nth_percent*length(t))));
    q = sort_segment_cy5{n,j};
    sintensity_segment_cy5{n,j} = mean(q(1:floor(nth_percent*length(q))));
    clearvars q t;
    end
end 

% For dilated, no background subtraction
for n=1:nframes % For each FRAME of video
    for j = 1:length(regionprops_pixel_dilate_dapi{1,n}) % For Each Region Found in the frame
    sort_dilate_dapi{n,j} = sort(regionprops_pixel_dilate_dapi{1,n}(j).PixelValues,'descend');  
    sort_dilate_cy5{n,j} = sort(regionprops_pixel_dilate_cy5{1,n}(j).PixelValues, 'descend');
    % Find Average of predefined intensities
    t = sort_dilate_dapi{n,j};
    sintensity_dilate_dapi{n,j} = mean(t(1:floor(nth_percent*length(t))));
    q = sort_dilate_cy5{n,j};
    sintensity_dilate_cy5{n,j} = mean(q(1:floor(nth_percent*length(q))));
    clearvars q t;
    end
end

% For segmented, background subtracted
for n=1:nframes % For each FRAME of video
    for j = 1:length(regionprops_pixel_segment_dapi_b{1,n}) % For Each Region Found in the frame
    sort_segment_dapi_b{n,j} = sort(regionprops_pixel_segment_dapi_b{1,n}(j).PixelValues,'descend');  
    sort_segment_cy5_b{n,j} = sort(regionprops_pixel_segment_cy5_b{1,n}(j).PixelValues, 'descend');
    % Find Average of predefined intensities
    t = sort_segment_dapi_b{n,j};
    sintensity_segment_dapi_b{n,j} = mean(t(1:floor(nth_percent*length(t))));
    q = sort_segment_cy5_b{n,j};
    sintensity_segment_cy5_b{n,j} = mean(q(1:floor(nth_percent*length(q))));
    clearvars q t;
    end
end 

% For dilated, background subtracted 
for n=1:nframes % For each FRAME of video
    for j = 1:length(regionprops_pixel_dilate_dapi_b{1,n}) % For Each Region Found in the frame
    sort_dilate_dapi_b{n,j} = sort(regionprops_pixel_dilate_dapi_b{1,n}(j).PixelValues,'descend');  
    sort_dilate_cy5_b{n,j} = sort(regionprops_pixel_dilate_cy5_b{1,n}(j).PixelValues, 'descend');
    % Find Average of predefined intensities
    t = sort_dilate_dapi_b{n,j};
    sintensity_dilate_dapi_b{n,j} = mean(t(1:floor(nth_percent*length(t))));
    q = sort_dilate_cy5_b{n,j};
    sintensity_dilate_cy5_b{n,j} = mean(q(1:floor(nth_percent*length(q))));
    clearvars q t;
    end
end

%% Linearize outputs from frames format
% Set variables for n=1 
% For segmented, no background subtraction
for n = 1:length(regionprops_pixel_segment_dapi{1,1})
intensity_adj_segment_dapi{n} =  sintensity_segment_dapi{1,n};
intensity_adj_segment_cy5{n} = sintensity_segment_cy5{1,n}; 
end

% For dilated, no background subtraction
for n = 1:length(regionprops_pixel_dilate_dapi{1,1})
intensity_adj_dilate_dapi{n} =  sintensity_dilate_dapi{1,n};
intensity_adj_dilate_cy5{n} = sintensity_dilate_cy5{1,n}; 
end

% For segmented, with background subtraction
for n = 1:length(regionprops_pixel_segment_dapi_b{1,1})
intensity_adj_segment_dapi_b{n} =  sintensity_segment_dapi_b{1,n};
intensity_adj_segment_cy5_b{n} = sintensity_segment_cy5_b{1,n}; 
end

% For dilated, with background subtraction
for n = 1:length(regionprops_pixel_dilate_dapi_b{1,1})
intensity_adj_dilate_dapi_b{n} =  sintensity_dilate_dapi_b{1,n};
intensity_adj_dilate_cy5_b{n} = sintensity_dilate_cy5_b{1,n}; 
end

%% Create Linearized Variables 
% For segmented, no background subtraction
for n = 2:nframes
     j = length(regionprops_pixel_segment_dapi{1,n});
intensity_adj_segment_dapi = [intensity_adj_segment_dapi sintensity_segment_dapi{n,1:j}];
intensity_adj_segment_cy5 = [intensity_adj_segment_cy5 sintensity_segment_cy5{n,1:j}];
clearvars j;    
end

% For dilated, no background subtraction
for n = 2:nframes
     j = length(regionprops_pixel_dilate_dapi{1,n});
intensity_adj_dilate_dapi = [intensity_adj_dilate_dapi sintensity_dilate_dapi{n,1:j}];
intensity_adj_dilate_cy5 = [intensity_adj_dilate_cy5 sintensity_dilate_cy5{n,1:j}];
clearvars j;    
end

% For segmented, with background subtraction
for n = 2:nframes
     j = length(regionprops_pixel_segment_dapi_b{1,n});
intensity_adj_segment_dapi_b = [intensity_adj_segment_dapi_b sintensity_segment_dapi_b{n,1:j}];
intensity_adj_segment_cy5_b = [intensity_adj_segment_cy5_b sintensity_segment_cy5_b{n,1:j}];
clearvars j;    
end

% For dilated, with background subtraction
for n = 2:nframes
     j = length(regionprops_pixel_dilate_dapi_b{1,n});
intensity_adj_dilate_dapi_b = [intensity_adj_dilate_dapi_b sintensity_dilate_dapi_b{n,1:j}];
intensity_adj_dilate_cy5_b = [intensity_adj_dilate_cy5_b sintensity_dilate_cy5_b{n,1:j}];
clearvars j;    
end

%% Create Scatter Plots for X% pixel Intensity Averages

% Segmented Graphs
figure; scatter(cell2mat(intensity_adj_segment_dapi), cell2mat(intensity_adj_segment_cy5));
title('Plot of segmented data with no background reduction AND pixel Cor')
figure; scatter(cell2mat(intensity_adj_segment_dapi_b), cell2mat(intensity_adj_segment_cy5_b));
title('Plot of segmented data WITH background reduction AND pixel Cor')

% Dilated Graphs
figure; scatter(cell2mat(intensity_adj_dilate_dapi), cell2mat(intensity_adj_dilate_cy5));
title('Plot of Dilated data with no background reduction AND pixel Cor')
figure; scatter(cell2mat(intensity_adj_dilate_dapi_b), cell2mat(intensity_adj_dilate_cy5_b));
title('Plot of Dilated data WITH background reduction AND pixel Cor')