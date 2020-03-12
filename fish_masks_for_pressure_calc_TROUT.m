%February 22, 2017

%This code is a modified version of John Dabiri's
%autokinematicsandpressureintegrator_nopressure.m script.


%This version will read in digitized outlines of a fish, plus the vector fields,
%and will generate blanking boundaries for queen2.

%Two modes: Mode 1 allows for iteratively changing the boundary for the
%fish, plotting the boundaries alongside the vector coordinates, and
%determining if the boundaries actually enclose vectors.  Once boundary
%dilation settings are finalized, an animation of the plots can be saved
%out as proof of working.  Mode 2 uses the settings chosen in Mode 1,
%generates the interface files for use with queen2, and saves an animation
%of these interfaces plotted over the movie frame, to show that the
%interfaces match the video.



clear all
close all

j=0;                                                         % index for movie storage

%Set metadata: filepath information, individual fish to analyze, scale factors,
%location of axis origins in DaVis, numberformat in the filenames,
%delimiters, and increment between frames used for pressure code.

vecbase = 'H:/fish PIV vector stacks_v2/trout/';
individuals = {'trout4'};

movbase = 'H:/fish vids for pressure fields/trout/clips for steady analysis_long_BCF_v2/';

interfacebase = 'H:/fish_outlines_v2/';
maskbase = 'F:/fish masks_v2/';

pixscale = 0.252564;
DaVis_x0 = 198;
DaVis_y0 = 952;
xmin = -DaVis_x0*pixscale;
xmax = xmin+(1024*pixscale);
ymin = (-1024+DaVis_y0)*pixscale;
ymax = ymin+(1024*pixscale);

numformat = '%05d';                                         % format of numbers in file name. '%05d' is five digit format with leading zeros
delimiter = '\t';                                             % delimiter between columns in vector file
headerlines = 1;                                            % number of header lines in vector file
increment = 10;                                              % increment between file numbers

%loop through each individual which will be analyzed.
for i = 1:1:size(individuals,2),
    
    %Find all sequences for that individual
%     d = dir([vecbase individuals{1,i}]);
    d = dir([interfacebase individuals{1,i}]);
    isub = [d(:).isdir];
    folders = {d(isub).name};
    folders(ismember(folders,{'.','..'}))=[];
    
    %For each sequence (saved to unique folders),
    for f = 1:1:size(folders,2),
        
        %Report which sequence we're on - for tracking code progress
        folders{f}
        clear pelvics
        
        %Extract the first and last frame number from the filename
        frameStr = strsplit(folders{f},'frm');
        frameStr = strsplit(frameStr{2},'-');
        first = str2num(frameStr{1});
        last = str2num(frameStr{2});
        
%         if strcmp(folders{f},'trout7_2.5BLs_1_ventral_frm43-240') || strcmp(folders{f},'trout9_2.5BLs_3_ventral_frm819-1048'),
%             last = last + 18;
%         end
        
        %Set filepath and name base for the corresponding cropped vector fields
        cropstub = [vecbase individuals{i} '/' folders{f} '/cropped_stack/B'];
        
        %Find the movie for the sequence, set first frame in movie = 1, and
        %increment for reading frames = 10
        if strcmp(folders{f},'trout7_2.5BLs_1_ventral_frm43-240'),
            mov = VideoReader('H:/fish vids for pressure fields/trout/trout7/trout7_2.5BLs_1_ventral.avi');
            fr = first;
        elseif strcmp(folders{f},'trout9_2.5BLs_3_ventral_frm819-1048'),
            mov = VideoReader('H:/fish vids for pressure fields/trout/trout9/trout9_2.5BLs_3_ventral.avi');
            fr = first;
        else
            mov = VideoReader([movbase individuals{i} '/' folders{f} '.avi']);
            fr = 1;
        end
        frstep = increment;
        
        %Set filepath and name base for manually-digitized fish outlines
        interfacestub = [interfacebase individuals{i} '/' folders{f} '/' folders{f} '_'];
        
        %Make sure folder where will save final masks exists, and if not,
        %create it
        maskfolder = [maskbase individuals{i} '/' folders{f}];
        if ~exist(maskfolder,'dir'),
            mkdir(maskfolder);
        end
        
        %For every 10th frame in the sequence,
        for filenum = first:increment:last-1                          % loop through files                                                 
                clf                                                     % clear current figure
                clear pelvicFin1
                clear pelvicFin2
                filenum                                                 % display current file number

            im=flipud(read(mov,fr));    %Get current image, flip upside-down to move origin from upper left to lower left

            cropped = dlmread([cropstub num2str(filenum,numformat) '.dat'],delimiter,headerlines,0); % read in vector file
            interface = xlsread([interfacestub num2str(filenum,numformat) '.xlsx'],'Outline');         %read in the fish outlines

            
            %Split fish outline data to x and y lists, move origin to match
            %origin in vec files
            xdata = interface(:,2)+xmin;
            ydata = interface(:,4)+ymin;

            %Add the first point to the end of the list to make a closed shape
            xdata(end+1) = xdata(1,1);
            ydata(end+1) = ydata(1,1);



            figure(1)
            imshow(im,'XData',[xmin xmax],'YData',[ymin ymax])    % display current movie frame, scaled to coordinates in pressure file
            %axis equal
            axis on
            set(gca,'YDir','normal')                                                    %set axes to correct direction (because they are flipped in images)
            hold on    %should see correct axes here
            %axis([xmin xmax ymin ymax])
            axis square
%             return




            %Manually adjust selected areas of the boundaries - ex: tail is too thin,
            %but rest of body came out okay.  Generally, either this or previous will
            %be ran.  However, if both are needed, leave the first 7 lines here
            %commented.

            % %Get the fish boundaries to adjust
            % BWsfinal = A;
            % BWsfinal(end-10:end,:) = 0;                                                 %Set edges to black
            % BWsfinal(:,end-10:end) = 0;
            % BWsfinal(:,1:10) = 0;
            % BWsfinal(1:10,:) = 0;
            % BWsfinal = flipud(BWsfinal);
            % B=bwboundaries(BWsfinal);

            %Set the x-value above which adjustments will occur by changing the
            %*[decimal] term.  Ex: (max-min)*0.7 means that adjustments will occur
            %in the last 30% of the fish's length (after 70%).
            adjx = round(range(xdata)*0.8) + min(xdata);

            %Optionally also set a leading edge adjustment - ex: if fish's nose is too
            %thin
            adjxLE = round(range(xdata)*0.05) + min(xdata);

            %Find maximum x value on the fish, and a trigger that will go off when
            %we've hit that value
            xvalmax = max(xdata);
            trigger = -1;


            %Scroll through all of the x-values in the boundary points list.  Note that the
            %list is formatted so that you move clockwise from the fish's nose.  For each
            %one,
            for k=1:length(xdata),
                %Save current value
                xval = xdata(k);
                %If current value is equal to the maximum, trigger.
                if xval == xvalmax,
                    trigger = 0;
                    xdata(k) = xdata(k) + 1;
                    if abs(xdata(k)-xdata(k-1)) < abs(xdata(k+1)-xdata(k)),
                        ydata(k) = ydata(k) - 0.75;
                    else
                        ydata(k) = ydata(k) + 0.75;
                    end
                elseif trigger == 0,
                    trigger = 1;
                end
                %If current x value is bigger than our limit for adjustment, and we're
                %still increasing x values (haven't hit the maximum yet), we're on the
                %top side of the fish.  Adjust the y-coordinate for this point upward.
                %If we're past the maximum (triggered), we're on the bottom side of the
                %fish.  Adjust the y-coordinate downward.  Adjustments are made by some
                %number of pixels.
                if xval > adjx,
                    if trigger < 0,
                        ydata(k)=ydata(k)+0.75;
                    elseif trigger > 0,
                        ydata(k)=ydata(k)-0.75;
                    end
                end
                %If we also need to adjust at the leading edge, turn this part on.
                %This adjusts the x-values at the leading edge leftward, extending the
                %nose of the fish.
                if xval < adjxLE,
                    xdata(k)=xdata(k)-0.5;
                end


            end


            %Turn on for Mode 1 - seeing if the boundary contains vectors.  Red dots
            %are the locations of the vectors.
            plot(xdata,ydata)
            hold on
            plot(cropped(:,1),cropped(:,2),'.r','markers',4)
            %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.6 1]);
            %pause(1.0)


%             return
            


            % % %Turn on to convert x and y data from mm to m, if needed.
            % boundx=boundx/1000;
            % boundy=boundy/1000;

            
            %Turn on for Mode 2 - saving the final interfaces.
            print([maskfolder '/mask_' num2str(filenum,numformat)], '-dtiffn','-r300');
            dlmwrite([maskfolder '/mask_' num2str(filenum,numformat) '.dat'], cat(1,zeros(1,2),cat(2,xdata,ydata)), '\t');
%             return

            fr=fr+frstep;                                                          %Increment movie frame number
            %pause(0.5)                                                              % pause to see output
            
            
        end
        
        %Turn on to save out animation
        
        d2 = dir(maskfolder);
        isub2 = [d2(:).isdir];
        isub2 = abs(1-isub2);
        numFrames = sum(isub2)/2;
        
        currFr = first;
        v = VideoWriter([maskfolder '/mask_incr10_vectorproof.avi'],'Uncompressed AVI');
        open(v);
        for frame = 1:1:numFrames,
            im = imread([maskfolder '/mask_' num2str(currFr,numformat) '.tif']);
            writeVideo(v,im);
            currFr = currFr + 10;
        end
        close(v);
        
%         return
        
    end

end

