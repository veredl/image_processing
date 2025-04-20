%% Reading experiment details of the series of images in the current folder

details = readlines('details.txt');

cell_line = details(1);
WT_or_KO = details(2);

% First chromosome
first_ch = details(3); 
first_ch_details = split(first_ch);
first_ch_number = first_ch_details(1);
first_ch_color = first_ch_details(2);

% Second chromosome
second_ch = details(4);
second_ch_details = split(second_ch);
second_ch_number = second_ch_details(1);
second_ch_color = second_ch_details(2);

series_name = convertStringsToChars(details(5));

% The total number of images in the series (from 3 channels)
tot_number_of_images = str2num(details(6));

% The total number of RGB images (3 channels combined)
number_of_images = tot_number_of_images/3;

% The threshold for blue channel chosen for this cell_line images' series
threshold_for_dapi = str2num(details(7));

% Resolution of all images
resolution = [1002 1004];

% The smallest size chosen for a nucleus to include in segmentation
smallest_size_nu = 1000;

% The smallest size chosen for a chromosome to include in segmentation
smallest_size_ch = 10;

% The smallest size chosen for a chromosome to include in measuring
% its circularity
smallest_size_ch_circ = 150;

%% Reading the raw data

cd('Image processing/1 raw data');

i=1;
j=1;

while i<=tot_number_of_images

    blue_images(:,:,j) = imread([series_name, sprintf('%03d', i),'.tif']);
    i = i+1;
    green_images(:,:,j) = imread([series_name, sprintf('%03d', i),'.tif']);
    i = i+1;
    red_images(:,:,j) = imread([series_name, sprintf('%03d', i),'.tif']);
    i = i+1;

    j = j+1;

end

cd('../..');

%% Printing the raw data to folder, as normalized to min - max values 
% 
% cd('Image processing/2 normalized data');
% 
% blue_images_normalized = mat2gray(blue_images);
% green_images_normalized = mat2gray(green_images);
% red_images_normalized = mat2gray(red_images);
% 
% i=1;
% j=1;
% 
% while i<=tot_number_of_images
% 
%     imwrite(blue_images_normalized(:,:,j),['normalized_', series_name, sprintf('%03d', i),'.tif']);
%     i = i+1;
%     imwrite(green_images_normalized(:,:,j),['normalized_', series_name, sprintf('%03d', i),'.tif']);
%     i = i+1;
%     imwrite(red_images_normalized(:,:,j),['normalized_', series_name, sprintf('%03d', i),'.tif']);
%     i = i+1;
% 
%     j = j+1;
% 
% end
% 
% cd('../..');

%% Segmentation

for image=1:number_of_images
    
    % Blue channel
    blue_image = blue_images(:,:,image);
    blue_image_as_bw = im2bw(blue_image, threshold_for_dapi);
    blue_image_smooth = medfilt2(blue_image_as_bw);
    blue_image_noborder = imclearborder(blue_image_smooth);
    segmented_B = bwconncomp(blue_image_noborder);

    i=1; 
    while i<=length(segmented_B.PixelIdxList)

        if length(segmented_B.PixelIdxList{1,i})<smallest_size_nu
            segmented_B.PixelIdxList(:,i) = [];
        else i=i+1;
        end

    end 
    
    blue_nuclei(image,1:length(segmented_B.PixelIdxList)) = segmented_B.PixelIdxList;

    how_many_nuclei(image) = length(segmented_B.PixelIdxList);

    % Green channel
    green_image = green_images(:,:,image);
     
    for nu=1:length(segmented_B.PixelIdxList)

        green_ch_in_nu = uint16(false(resolution));
        green_ch_in_nu(blue_nuclei{image,nu}) = green_image(blue_nuclei{image,nu});
        
        green_thresh_for_nu = graythresh(green_image(blue_nuclei{image,nu}));

        green_image_as_bw = im2bw(green_ch_in_nu, green_thresh_for_nu);
        green_image_smooth = medfilt2(green_image_as_bw);
        segmented_G = bwconncomp(green_image_smooth);
        
        i=1;
        while i<=length(segmented_G.PixelIdxList)

            if length(segmented_G.PixelIdxList{1,i})<smallest_size_ch
                segmented_G.PixelIdxList(:,i) = [];
            else i=i+1;
            end

        end 
        
        green_chromosomes{image,nu} = segmented_G.PixelIdxList;
    end
    
    % Red channel
    red_image = red_images(:,:,image);
     
    for nu=1:length(segmented_B.PixelIdxList)

        red_ch_in_nu = uint16(false(resolution));
        red_ch_in_nu(blue_nuclei{image,nu}) = red_image(blue_nuclei{image,nu});
        
        red_thresh_for_nu = graythresh(red_image(blue_nuclei{image,nu}));
        
        red_image_as_bw = im2bw(red_ch_in_nu, red_thresh_for_nu);
        red_image_smooth = medfilt2(red_image_as_bw);
        segmented_R = bwconncomp(red_image_smooth);
        
        i=1;
        while i<=length(segmented_R.PixelIdxList)

            if length(segmented_R.PixelIdxList{1,i})<smallest_size_ch
                segmented_R.PixelIdxList(:,i) = [];
            else i=i+1;
            end

        end 
        
        red_chromosomes{image,nu} = segmented_R.PixelIdxList;
    end
    
end

%% Printing the segmented images to "segmented data" folder
%
% cd('Image processing/3 segmented data');
% 
% for image=1:number_of_images
% 
%     nucleus = false(resolution);        
%     for nu=1:how_many_nuclei(image)    
%         nucleus(blue_nuclei{image,nu}) = true;
%     end
% 
%     ch_green = false(resolution);
%     for nu=1:how_many_nuclei(image)       
%         for ch=1:length(green_chromosomes{image,nu})
%             ch_green(green_chromosomes{image,nu}{1,ch}) = true;
%         end
%     end
% 
%     ch_red = false(resolution);
%     for nu=1:how_many_nuclei(image)       
%         for chr=1:length(red_chromosomes{image,nu})
%             ch_red(red_chromosomes{image,nu}{1,chr}) = true;
%         end
%     end
% 
%     rgbImage = double(cat(3, ch_red, ch_green, nucleus));
%     imwrite(rgbImage,['Segmented_', series_name ,sprintf('%03d', image),'.tif']);
%
% end
%
% cd('../..');

%% Removing artifacts (reading from a txt file prepared in advance)

blue_nuclei_cleaned = blue_nuclei;
green_chromosomes_cleaned = green_chromosomes;
red_chromosomes_cleaned = red_chromosomes;

[images,max_nuclei] = size(blue_nuclei_cleaned);

file = strcat('artifacts_', cell_line, '_', WT_or_KO, '.txt')

artifacts = readcell(file);
[nr,nc] = size(artifacts);

for i=1:nr

    image = artifacts{i, 1};
    clear_nu = artifacts{i, 2};

    if clear_nu<max_nuclei
        blue_nuclei_cleaned(image,clear_nu:(end-1)) = blue_nuclei_cleaned(image,(clear_nu+1):end);
        green_chromosomes_cleaned(image,clear_nu:(end-1)) = green_chromosomes_cleaned(image,(clear_nu+1):end);
        red_chromosomes_cleaned(image,clear_nu:(end-1)) = red_chromosomes_cleaned(image,(clear_nu+1):end);
    end
    
    blue_nuclei_cleaned(image,end) = {[]};
    green_chromosomes_cleaned(image,end) = {[]};
    red_chromosomes_cleaned(image,end) = {[]};

    how_many_nuclei(image) = how_many_nuclei(image)-1;

end

%% Printing the images to "segmented data witout artifacts" folder
% 
% cd('Image processing/4 segmented data without artifacts');
% 
% for image=1:number_of_images
%    
%     nucleus = false(resolution);        
%     for nu=1:how_many_nuclei(image)    
%         nucleus(blue_nuclei_cleaned{image,nu}) = true;
%     end
%    
%     ch_green = false(resolution);
%     for nu=1:how_many_nuclei(image)       
%         for ch=1:length(green_chromosomes_cleaned{image,nu})
%             ch_green(green_chromosomes_cleaned{image,nu}{1,ch}) = true;
%         end
%     end
%    
%     ch_red = false(resolution);
%     for nu=1:how_many_nuclei(image)       
%         for ch=1:length(red_chromosomes_cleaned{image,nu})
%             ch_red(red_chromosomes_cleaned{image,nu}{1,ch}) = true;
%         end
%     end
%    
%     rgbImage = double(cat(3, ch_red, ch_green, nucleus));
%     imwrite(rgbImage,['processed_segmented_', series_name ,sprintf('%03d', image),'.tif']);
%
% end
% 
% cd('../..');

%% Building a struct of chromosomes' properties 

for image=1:number_of_images
    for nu=1:how_many_nuclei(image) 
    
    % Area 
    sum_g = 0;
    sum_r = 0;
    
    for ch=1:length(green_chromosomes_cleaned{image,nu})
        sum_g = sum_g+length(green_chromosomes_cleaned{image,nu}{1,ch});
    end
    
    for ch=1:length(red_chromosomes_cleaned{image,nu})
        sum_r = sum_r+length(red_chromosomes_cleaned{image,nu}{1,ch});
    end
    
    props.area.blue_area(image,nu) = length(blue_nuclei_cleaned{image,nu});
    
    props.area.green_area(image,nu) = sum_g;
    props.area.green_area_devided_by_nuclei(image,nu) = sum_g/length(blue_nuclei_cleaned{image,nu});
    
    props.area.red_area(image,nu) = sum_r;
    props.area.red_area_devided_by_nuclei(image,nu) = sum_r/length(blue_nuclei_cleaned{image,nu});
    
    % Circularity
    for chg=1:length(green_chromosomes_cleaned{image,nu})
        peri_frame = false(resolution); 
        peri_frame(green_chromosomes_cleaned{image,nu}{1,chg}) = true;
        
        props.Circularity.green_perimeter.idx{image,nu}{1,chg} = find(bwperim(peri_frame)==1);
        props.Circularity.green_perimeter.length{image,nu}{1,chg} = length(props.Circularity.green_perimeter.idx{image,nu}{1,chg});
        perimeter = props.Circularity.green_perimeter.length{image,nu}{1,chg};
        
        props.Circularity.green_area{image,nu}{1,chg} = length(green_chromosomes_cleaned{image,nu}{1,chg});
        area = props.Circularity.green_area{image,nu}{1,chg};
        
        props.Circularity.green_circ{image,nu}{1,chg} = (perimeter.^ 2) ./ (4 * pi * area);
    end
    
    for chr=1:length(red_chromosomes_cleaned{image,nu})
        peri_frame = false(resolution); 
        peri_frame(red_chromosomes_cleaned{image,nu}{1,chr}) = true;
        
        props.Circularity.red_perimeter.idx{image,nu}{1,chr} = find(bwperim(peri_frame)==1);
        props.Circularity.red_perimeter.length{image,nu}{1,chr} = length(props.Circularity.red_perimeter.idx{image,nu}{1,chr});
        perimeter = props.Circularity.red_perimeter.length{image,nu}{1,chr};
        
        props.Circularity.red_area{image,nu}{1,chr} = length(red_chromosomes_cleaned{image,nu}{1,chr});
        area = props.Circularity.red_area{image,nu}{1,chr};
        
        props.Circularity.red_circ{image,nu}{1,chr} = (perimeter.^ 2) ./ (4 * pi * area);
    end
    
    % Co-localization 
    A = [];
    for chr=1:length(green_chromosomes_cleaned{image,nu})
        A = [A,green_chromosomes_cleaned{image,nu}{1,chr}'];
    end
    
    B = [];
    for chr=1:length(red_chromosomes_cleaned{image,nu})
        B = [B,red_chromosomes_cleaned{image,nu}{1,chr}'];
    end
    
    props.Coloc.tot(image,nu) = length(intersect(A,B));
    props.Coloc.idx{image,nu} = intersect(A,B);
    props.Coloc.normalized(image,nu) = props.Coloc.tot(image,nu)/props.area.blue_area(image,nu);
    
    end
end

%% Printing the overlap area of chromosomes to folder
% 
% cd('Image processing/5 circularity and colocalization');
%
% for image=1:number_of_images
% 
%     nucleus = false(resolution);        
%     for nu=1:how_many_nuclei(image)    
%         nucleus(blue_nuclei_cleaned{image,nu}) = true;
%     end
% 
%     ch_green = false(resolution);
%     for nu=1:how_many_nuclei(image)       
%         ch_green(props.Coloc.idx{image,nu}) = true;
%     end
% 
%     ch_red = false(resolution);
%     for nu=1:how_many_nuclei(image)       
%         ch_red(props.Coloc.idx{image,nu}) = true;
%     end
% 
%     rgbImage = double(cat(3, ch_red, ch_green, nucleus));
%     imwrite(rgbImage,['coloc_', series_name, sprintf('%03d', image),'.tif']);
% 
% end
% 
% cd('../..');

%% Printing the perimeters of chromosomes to folder
% 
% cd('Image processing/5 circularity and colocalization');
%
% for image=1:number_of_images
% 
%     nucleus = false(resolution);        
%     for nu=1:how_many_nuclei(image)    
%         nucleus(blue_nuclei_cleaned{image,nu}) = true;
%     end
% 
%     ch_green = false(resolution);
%     for nu=1:how_many_nuclei(image) 
%         for chg=1:length(green_chromosomes_cleaned{image,nu})    
%             ch_green(props.Circularity.green_perimeter.idx{image,nu}{1,chg}) = true;
%         end
%     end
% 
%     ch_red = false(resolution);
%     for nu=1:how_many_nuclei(image) 
%         for chr=1:length(red_chromosomes_cleaned{image,nu})
%             ch_red(props.Circularity.red_perimeter.idx{image,nu}{1,chr}) = true;
%         end
%     end
% 
%     rgbImage = double(cat(3, ch_red, ch_green, nucleus));
%     imwrite(rgbImage,['perimeter_', series_name, sprintf('%03d', image),'.tif']);
% 
% end
% 
% cd('../..');

%% Export data of chromosomes' properties to txt files 

% area

% green
green_area_final = props.area.green_area_devided_by_nuclei;
idx_plus = green_area_final > 0;

% red
red_area_final = props.area.red_area_devided_by_nuclei;

% circularity

% green
circ_green_final = [];
for image=1:number_of_images        
    for nu=1:how_many_nuclei(image)  
        for chg=1:length(green_chromosomes_cleaned{image,nu})
            if length(green_chromosomes_cleaned{image,nu}{1,chg})>smallest_size_ch_circ
                circ_green_final = [circ_green_final, props.Circularity.green_circ{image,nu}{1,chg}];
            end
        end
    end
end

% red
circ_red_final = [];
for image=1:number_of_images        
    for nu=1:how_many_nuclei(image)  
        for chg=1:length(red_chromosomes_cleaned{image,nu})
            if length(red_chromosomes_cleaned{image,nu}{1,chg})>smallest_size_ch_circ
                circ_red_final = [circ_red_final, props.Circularity.red_circ{image,nu}{1,chg}];
            end
        end
    end
end

% coloc

coloc_final = [];
for image=1:number_of_images
    coloc_final = [coloc_final, props.Coloc.normalized(image,1:how_many_nuclei(image))];
end


d = dictionary([first_ch_color, second_ch_color], [first_ch_number, second_ch_number])

chr_number_green = convertStringsToChars( d('green') )
chr_number_red = convertStringsToChars( d('red') )

% area
writematrix( num2str(green_area_final(idx_plus), '%.5f '), strcat('area_', chr_number_green, '_', WT_or_KO,'.txt') )
writematrix( num2str(red_area_final(idx_plus), '%.5f '), strcat('area_', chr_number_red, '_', WT_or_KO,'.txt') )

% circ
writematrix( num2str(circ_green_final', '%.5f '), strcat('circ_', chr_number_green, '_', WT_or_KO,'.txt') )
writematrix( num2str(circ_red_final', '%.5f '), strcat('circ_', chr_number_red, '_', WT_or_KO,'.txt') )

% coloc
writematrix( num2str(coloc_final', '%.5f '), strcat('coloc_', chr_number_green, '_', chr_number_red, '_', WT_or_KO,'.txt') )
