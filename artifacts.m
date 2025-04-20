%% Creating artifacts.txt file (Removing artifacts manually)

file = strcat('artifacts_', cell_line, '_', WT_or_KO, '.txt');

blue_nuclei_cleaned = blue_nuclei;
green_chromosomes_cleaned = green_chromosomes;
red_chromosomes_cleaned = red_chromosomes;

%% Reaching an image 

image = 100; % choose which image

nucleus = false(resolution);        
for nu=1:How_many_nuclei(image,1)    
    nucleus(blue_nuclei_cleaned{image,nu}) = true;
end

ch_green = false(resolution);
for nu=1:How_many_nuclei(image,1)       
    for chg=1:length(green_chromosomes_cleaned{image,nu})
        ch_green(green_chromosomes_cleaned{image,nu}{1,chg}) = true;
    end
end

ch_red = false(resolution);
for nu=1:How_many_nuclei(image,1)       
    for chr=1:length(red_chromosomes_cleaned{image,nu})
        ch_red(red_chromosomes_cleaned{image,nu}{1,chr}) = true;
    end
end

figure
rgbImage = double(cat(3, ch_red , ch_green , nucleus));
imshow(rgbImage)

%% Reaching a nucleus in the previous image

nu = 8; % choose which nucleus

nucleus = false(resolution);
nucleus(blue_nuclei_cleaned{image,nu}) = true;

ch_green = false(resolution);
for ch=1:length(green_chromosomes_cleaned{image,nu})
    ch_green(green_chromosomes_cleaned{image,nu}{1,ch}) = true;
end

ch_red = false(resolution);
for ch=1:length(red_chromosomes_cleaned{image,nu})
    ch_red(red_chromosomes_cleaned{image,nu}{1,ch}) = true;
end

figure
rgbImage = double(cat(3, ch_red , ch_green , nucleus));
imshow(rgbImage)

%% Clear this nucleus + all its chromosomes from the data 

clear_nu = nu; % Remove this nucleus (artifact)

if clear_nu<length(blue_nuclei_cleaned(image,:))
    blue_nuclei_cleaned(image,clear_nu:(end-1)) = blue_nuclei_cleaned(image,(clear_nu+1):end);
    green_chromosomes_cleaned(image,clear_nu:(end-1)) = green_chromosomes_cleaned(image,(clear_nu+1):end);
    red_chromosomes_cleaned(image,clear_nu:(end-1)) = red_chromosomes_cleaned(image,(clear_nu+1):end);
end

blue_nuclei_cleaned(image,end) = {[]};
green_chromosomes_cleaned(image,end) = {[]};
red_chromosomes_cleaned(image,end) = {[]};

frame_and_nucleus = [image, nu]

if isfile(file) & dir(file).bytes ~= 0
     writematrix(frame_and_nucleus, file, 'WriteMode', 'append');
else
     writematrix(frame_and_nucleus, file);
end


% how many nuclei in each image
for i=1:length(blue_nuclei_cleaned(:,1))
    if isempty(blue_nuclei_cleaned{i,end})==0
        How_many_nuclei(i,1) = length(blue_nuclei_cleaned(i,:));
    else
        How_many_nuclei(i,1) = find(cellfun(@isempty, blue_nuclei_cleaned(i,:)), 1)-1;
    end
end

close all;

%% Correcting the image to its original form

image = 1; % choose which image

blue_nuclei_cleaned(image,:) = blue_nuclei(image,:);
green_chromosomes_cleaned(image,:) = green_chromosomes(image,:);
red_chromosomes_cleaned(image,:) = red_chromosomes(image,:);

fix_frame = ['frame number:  ', num2str(image), '  returned to its original. ' ...
    'pls clear manually the frame number from artifacts.txt file'];
disp(fix_frame)

% how many nuclei in each image
for i=1:length(blue_nuclei_cleaned(:,1))
    if isempty(blue_nuclei_cleaned{i,end})==0
        How_many_nuclei(i,1) = length(blue_nuclei_cleaned(i,:));
    else
        How_many_nuclei(i,1) = find(cellfun(@isempty, blue_nuclei_cleaned(i,:)), 1)-1;
    end
end
