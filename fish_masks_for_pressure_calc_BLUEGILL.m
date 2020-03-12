%April 12, 2017

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



% clear all
% close all

j=0;                                                         % index for movie storage

%Set metadata: filepath information, individual fish to analyze, scale factors,
%location of axis origins in DaVis, numberformat in the filenames,
%delimiters, and increment between frames used for pressure code.

vecbase = 'H:/fish PIV vector stacks_v2/bluegill/';
individuals = {'klbg2','klbg3','klbg4','klbg5'};

movbase = 'H:/fish vids for pressure fields/bluegill/clips for steady analysis_long_BCF_v2/';

interfacebase = 'H:/fish_outlines_v2/';
maskbase = 'F:/fish masks_v2/';

pixscale = 0.254432;
DaVis_x0 = 110;
DaVis_y0 = 881;
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
    %d = dir([vecbase individuals{1,i}]);
    d = dir([interfacebase individuals{1,i}]);
    isub = [d(:).isdir];
    folders = {d(isub).name};
    folders(ismember(folders,{'.','..'}))=[];
    
    %Get corresponding flow speed (useful for determining klbg3 source
    %folders)
    if strcmp(individuals{i},'klbg2'),
        flowBase = 0.1;
    elseif strcmp(individuals{i},'klbg4'),
        flowBase = 0.097;
    elseif strcmp(individuals{i},'klbg5'),
        flowBase = 0.096;
    elseif strcmp(individuals{i},'klbg7'),
        flowBase = 0.093;
    elseif strcmp(individuals{i},'klbg3'),
        flowBase = 0.115;
    end
    
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
        
        %Extract flow speed from filename and convert to meters/sec
        flowBL = strsplit(folders{f},'BLs');
        flowBL = strsplit(flowBL{1},'_');
        flow = str2num(flowBL{2})*flowBase;
        
        %Set filepath and name base for the corresponding cropped vector fields
        if strcmp(individuals{i},'klbg3'),
            teststr = strsplit(folders{f},'_');
            if ismember(str2num(flowBL{2}),[2.0,2.25,2.5]),
                if str2num(flowBL{2}) == 2.25 && ~strcmp(teststr{3},'1'),
                    cropstub = [vecbase individuals{i} '_10.03/' folders{f} '/cropped_stack/B'];
                else
                    cropstub = [vecbase individuals{i} '_09.15/' folders{f} '/cropped_stack/B'];
                end
            else
                cropstub = [vecbase individuals{i} '_10.03/' folders{f} '/cropped_stack/B'];
            end
        else
            cropstub = [vecbase individuals{i} '/' folders{f} '/cropped_stack/B'];
        end
        
        %Find the movie for the sequence, set first frame in movie = 1, and
        %increment for reading frames = 10
        mov = VideoReader([movbase individuals{i} '/' folders{f} '.avi']);
        fr = 1;
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

            
            %Check for pelvic fins
            [status, sheets] = xlsfinfo([interfacestub num2str(filenum,numformat) '.xlsx']);
            check = find(strcmp(sheets, 'Pelvic Fins'));
            if not(isempty(check)),
                pelvics = xlsread([interfacestub num2str(filenum,numformat) '.xlsx'],'Pelvic Fins');
            end

            if exist('pelvics','var'),
%                 %Use these plot features to verify that the pectoral/pelvic
%                 %overlap has been handled correctly.  Turn off for final
%                 %running of code.
%                 plot(xdata,ydata)
%                 hold on

                %Find the coordinates for upper(1) and lower(2) pelvic fins
                if pelvics(1) >0,
                    pelvicFin1 = interface(pelvics(1):pelvics(2),:);
                    pelvicFin1 = [pelvicFin1(:,2)+xmin pelvicFin1(:,4)+ymin];
                end
                if pelvics(3) > 0,
                    pelvicFin2 = interface(pelvics(3):pelvics(4),:);
                    pelvicFin2 = [pelvicFin2(:,2)+xmin pelvicFin2(:,4)+ymin];
                end
                
                %Check to see if the pelvic fins overlap with pectoral fins
                %(by seeing if the arc of the pelvic fin crosses the rest
                %of the body/fins outline).  Repeat for upper and lower
                %fins.
                if exist('pelvicFin1','var'),
                    p1 = InterX(pelvicFin1', [xdata(1:pelvics(1)) ydata(1:pelvics(1))]');
                    %Make sure the anteriormost point of the pelvic fin isn't marked
                    %as an overlap
                    rep = find(abs(p1(1,:) - pelvicFin1(1,1))<0.0001);
                    if not(isempty(rep)),
                        p1(:,rep) = [];
                    end
                end

                if exist('pelvicFin2','var'),
                    p2 = InterX(pelvicFin2', [xdata(pelvics(4):end) ydata(pelvics(4):end)]');
                    rep = find(abs(p2(1,:) - pelvicFin2(end,1))<0.0001);
                    %Make sure the anteriormost point of the pelvic fin isn't marked
                    %as an overlap
                    if not(isempty(rep)),
                        p2(:,rep) = [];
                    end
                end

                %find index of tail tip (used to split fish outline into
                %upper and lower halves later)
                [~, idx] = max(xdata);
                
                if exist('pelvicFin1','var'),
                    %If the upper pelvic fin overlaps with the upper pectoral
                    %fin,
                    if not(isempty(p1)),
                        %Find the points of the pelvic fin that are inside the
                        %body outline of the fish
                        internalP1 = inpolygon(pelvicFin1(:,1),pelvicFin1(:,2),[xdata(1:pelvics(1)-1);xdata((find(xdata(pelvics(2):end)<xdata(pelvics(1)),1,'first')+pelvics(2)-1):end)],[ydata(1:pelvics(1)-1);ydata((find(xdata(pelvics(2):end)<xdata(pelvics(1)),1,'first')+pelvics(2)-1):end)]);
                        %Find the points of the body outline that are inside
                        %the pelvic fin outline
                        internalB1 = inpolygon(xdata(1:pelvics(1)-1),ydata(1:pelvics(1)-1),xdata(pelvics(1):(find(xdata(pelvics(2):end)<xdata(pelvics(1)),1,'first')+pelvics(2)-1)),ydata(pelvics(1):(find(xdata(pelvics(2):end)<xdata(pelvics(1)),1,'first')+pelvics(2)-1)));
                        %Get all the points marked "true" (internal points)
                        internalP1 = find(internalP1 == 1);
                        %Add the number of points in the body outline prior to
                        %pelvic fin start so as to make indices match those in
                        %body outline
                        internalP1 = internalP1 + pelvics(1)-1;
                        %Get all the points marked "true" (internal points)
                        internalB1 = find(internalB1 == 1);
                        %Find the posteriormost point where the pelvic fin arc and
                        %body outline cross
                        [~, pidx] = max(p1(1,:));

                        %If there are pelvic fin points inside the body
                        %outline (pectoral fin passes behind pelvic fin),
                        if not(isempty(internalP1)),
                            if length(internalP1) == 2 && abs(internalP1(1)-pelvics(1))<0.001 && abs(internalP1(2)-pelvics(2))<0.001 && size(p1,2)==2,
                                upperx = [xdata(1:find(xdata(1:pelvics(1)-1)<p1(1,1),1,'last')); p1(1,1); xdata(pelvics(1)+1:pelvics(2)-1); p1(1,2); xdata(find(xdata(1:pelvics(1)-1)>p1(1,2),1,'first'):pelvics(1)-1); xdata(pelvics(2)+1:idx)];
                                uppery = [ydata(1:find(xdata(1:pelvics(1)-1)<p1(1,1),1,'last')); p1(2,1); ydata(pelvics(1)+1:pelvics(2)-1); p1(2,2); ydata(find(xdata(1:pelvics(1)-1)>p1(1,2),1,'first'):pelvics(1)-1); ydata(pelvics(2)+1:idx)];
                            %And if there are NOT body points inside the pelvic
                            %fin,
                            elseif isempty(internalB1),
                                if xdata(pelvics(1)-1) < xdata(pelvics(1)),
                                    firstval = pelvics(1)-1;
                                    val = firstval;
                                    while val > 1,
                                        %Check to see if the line segment between
                                        %current outline point and prev outline point
                                        %contains the posteriormost cross-over point.
                                        %Test is distance between segment and point.
                                        test = abs(det([[xdata(val);ydata(val)]-[xdata(val-1);ydata(val-1)],p1(:,pidx)-[xdata(val-1);ydata(val-1)]]))/norm([xdata(val);ydata(val)]-[xdata(val-1);ydata(val-1)]);
                                        %If segment doesn't contain 3rd
                                        %crossover, keep searching
                                        if test > 0.001,
                                            val = val-1;
                                        %If segment does contain crossover,
                                        %mark that we've hit the crossover and
                                        %stop searching
                                        else
                                            idxstop = val-1;
                                            break
                                        end
                                    end
                                    upperx = [xdata(1:idxstop); p1(1,pidx); xdata(internalP1(end)+1:idx)];
                                    uppery = [ydata(1:idxstop); p1(2,pidx); ydata(internalP1(end)+1:idx)];
                                else
                                    %Drop the pelvic fin points internal to the
                                    %body, replace with cross-over point
                                    upperx = [xdata(1:pelvics(1)-1); p1(1,pidx); xdata(internalP1(end)+1:idx)];
                                    uppery = [ydata(1:pelvics(1)-1); p1(2,pidx); ydata(internalP1(end)+1:idx)];
                                end
                            %Otherwise, if there ARE body points inside the
                            %pelvic fin,
                            else
                                %Drop the internal points and replace with the
                                %cross-over point
                                %If we have a simple cross-over,
                                if size(p1,2) < 3 && xdata(pelvics(2))>xdata(pelvics(1)-1),
                                    upperx = [xdata(1:internalB1(1)-1); p1(1,pidx); xdata(internalP1(end)+1:idx)];
                                    uppery = [ydata(1:internalB1(1)-1); p1(2,pidx); ydata(internalP1(end)+1:idx)];
                                elseif size(p1,2) == 4,
                                    lastbody = find(xdata(1:pelvics(1)-1) > p1(1,2),1,'first')-1;
                                    
                                    abovelowerlim = xdata(pelvics(1)+1:pelvics(2)) > p1(1,2);
                                    belowupperlim = xdata(pelvics(1)+1:pelvics(2)) < p1(1,3);
                                    both = pelvics(1)+find(abovelowerlim+belowupperlim == 2);
                                    
                                    firstval = lastbody+1;
                                    for val = lastbody:1:pelvics(2)-1,
                                        %Check to see if current point is an
                                        %internal point.  If so, skip.
                                        if any(abs(internalP1-val)<0.01) || any(abs(xdata(internalB1)-xdata(firstval))<0.01),
                                            firstval = val+1;
                                            continue
                                        end
                                        %Check to see if the line segment between
                                        %current outline point and next outline point
                                        %contains the posteriormost cross-over point.
                                        %Test is distance between segment and point.
                                        test = abs(det([[xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)],p1(:,pidx)-[xdata(val);ydata(val)]]))/norm([xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)]);
                                        %If segment doesn't contain 4th
                                        %crossover, keep searching
                                        if test > 0.001,
                                            continue
                                        %If segment does contain crossover,
                                        %mark that we've hit the crossover and
                                        %stop searching
                                        else
                                            lastonpeak = val;
                                            break
                                        end
                                    
                                    end

                                    %Get the index where we start following the
                                    %outline normally again
                                    idxcontoutline = find(xdata(pelvics(2):idx)>p1(1,pidx),1,'first') + pelvics(2)-1;

                                    for val=pelvics(1):1:pelvics(2)-1,
                                        if xdata(val) < xdata(idxcontoutline) && xdata(val) > p1(1,pidx),
                                            idxcontoutline = val;
                                            break
                                        end
                                    end

                                    upperx = [xdata(1:lastbody); p1(1,2); xdata(both); p1(1,3); xdata(firstval:lastonpeak); p1(1,pidx); xdata(idxcontoutline:idx)];
                                    uppery = [ydata(1:lastbody); p1(2,2); ydata(both); p1(2,3); ydata(firstval:lastonpeak); p1(2,pidx); ydata(idxcontoutline:idx)];
%                                     plot(upperx,uppery,'ro')
                                %If we have any part of pelvic fin appearing
                                %anterior to pectoral fin (assuming two peaks -
                                %peak 1 is pelvic, peak 2 is pect, and the
                                %peaks don't overlap)
                                else
                                    %Find last body point between first
                                    %crossover
                                    [~,minidx] = min(p1(1,:));
                                    lastbody = find(xdata(1:pelvics(1)-1) < p1(1,minidx),1,'last');

                                    abovelowerlim = xdata(pelvics(1):pelvics(2)) > p1(1,minidx);
                                    belowupperlim = xdata(pelvics(1):pelvics(2)) < p1(1,2);
                                    both = pelvics(1)-1 + find(abovelowerlim+belowupperlim == 2);

                                    firstval = lastbody+1;
                                    for val = lastbody:1:pelvics(2)-1,
                                        %Check to see if current point is an
                                        %internal point.  If so, skip.
                                        if any(abs(internalP1-val)<0.01) || any(abs(xdata(internalB1)-xdata(firstval))<0.01),
                                            firstval = val+1;
                                            continue
                                        end
                                        %Check to see if the line segment between
                                        %current outline point and next outline point
                                        %contains the posteriormost cross-over point.
                                        %Test is distance between segment and point.
                                        test = abs(det([[xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)],p1(:,pidx)-[xdata(val);ydata(val)]]))/norm([xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)]);
                                        %If segment doesn't contain 3rd
                                        %crossover, keep searching
                                        if test > 0.001,
                                            continue
                                        %If segment does contain crossover,
                                        %mark that we've hit the crossover and
                                        %stop searching
                                        else
                                            lastonpeak = val;
                                            break
                                        end
                                    
                                    end

                                    %Get the index where we start following the
                                    %outline normally again
                                    idxcontoutline = find(xdata(pelvics(2):idx)>p1(1,pidx),1,'first') + pelvics(2)-1;

                                    for val=pelvics(1):1:pelvics(2)-1,
                                        if xdata(val) < xdata(idxcontoutline) && xdata(val) > p1(1,pidx),
                                            idxcontoutline = val;
                                            break
                                        end
                                    end

                                    upperx = [xdata(1:lastbody); p1(1,minidx); xdata(both); p1(1,2); xdata(firstval:lastonpeak); p1(1,pidx); xdata(idxcontoutline:idx)];
                                    uppery = [ydata(1:lastbody); p1(2,minidx); ydata(both); p1(2,2); ydata(firstval:lastonpeak); p1(2,pidx); ydata(idxcontoutline:idx)];
%                                     plot(upperx,uppery,'ro')
                                end
%                                 return
                            end
                        %If there are NOT pelvic fin points inside the body outline
                        %(potentially the pectoral fin arcs back & crosses the
                        %pelvic fin, leaving a "hole" between the two fins),
                        else
                            %Set a trigger
                            trigger = 0;
                            %Go through the outline points one-by-one until we
                            %reach the end of the pelvic fin
                            for val = 1:1:length(xdata(1:pelvics(2))),
                                %Check to see if the line segment between
                                %current outline point and next outline point
                                %contains the posteriormost cross-over point.
                                %Test is distance between segment and point.
                                test = abs(det([[xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)],p1(:,pidx)-[xdata(val);ydata(val)]]))/norm([xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)]);

                                %If the segment contains the cross-over point
                                %and we haven't seen a matching segment before
                                %(trigger=0; this means we're on the pectoral
                                %fin and crossed into the pelvic fin),
                                if test <0.001 && trigger == 0,
                                    %Indicate that current outline point's
                                    %index is the start of where we have to
                                    %delete internal points
                                    start = val;
                                    %Mark that we've started finding internal
                                    %points
                                    trigger = 1;
                                %If the segment contains the cross-over point
                                %and we've seen a segment do this before 
                                %(trigger=1; this means that we're on the pelvic
                                %fin and have finished crossing the pectoral fin),
                                elseif test<0.001 && trigger == 1,
                                    %Note that current point is the last one to
                                    %be deleted
                                    stop = val+1;
                                    %Stop searching
                                    break
                                end
                            end
                            %Drop all internal points, replace with cross-over
                            %point
                            upperx = [xdata(1:start); p1(1,pidx); xdata(stop:idx)];
                            uppery = [ydata(1:start); p1(2,pidx); ydata(stop:idx)];
                        end
                    
                    %If no cross-over points were found
                    else
                        %If no cross-over points were found but internal body
                        %points exist (pectoral is completely behind pelvic),
                        internalB1 = inpolygon(xdata(1:pelvics(1)-1),ydata(1:pelvics(1)-1),xdata(pelvics(1):(find(xdata(pelvics(2):end)<xdata(pelvics(1)),1,'first')+pelvics(2)-1)),ydata(pelvics(1):(find(xdata(pelvics(2):end)<xdata(pelvics(1)),1,'first')+pelvics(2)-1)));
                        %Get all the points marked "true" (internal points)
                        internalB1 = find(internalB1 == 1);
                        if not(isempty(internalB1)),
                            upperx = [xdata(1:internalB1(1)-1); xdata(pelvics(1):idx)];
                            uppery = [ydata(1:internalB1(1)-1); ydata(pelvics(1):idx)];
                        else
                            upperx = xdata(1:idx);
                            uppery = ydata(1:idx);
                        end

                    end
                else
                    upperx = xdata(1:idx);
                    uppery = ydata(1:idx);
                end

%                 return
                if exist('pelvicFin2','var'),
                    %Repeat process for second pelvic fin.  Note that outline
                    %here runs from posterior to anterior, so have to flip
                    %search directions
                    if not(isempty(p2)),
                        internalP2 = inpolygon(pelvicFin2(:,1),pelvicFin2(:,2),[xdata(1:find(xdata(1:pelvics(3))<xdata(pelvics(4)),1,'last'));xdata(pelvics(4)+1:end)],[ydata(1:find(xdata(1:pelvics(3))<xdata(pelvics(4)),1,'last'));ydata(pelvics(4)+1:end)]);
                        internalB2 = inpolygon(xdata(pelvics(4)+1:end),ydata(pelvics(4)+1:end),xdata(find(xdata(1:pelvics(3))<xdata(pelvics(4)),1,'last')+1:pelvics(4)),ydata(find(xdata(1:pelvics(3))<xdata(pelvics(4)),1,'last')+1:pelvics(4)));
                        internalP2 = find(internalP2 == 1);
                        internalP2 = internalP2 + pelvics(3)-1;
                        internalB2 = find(internalB2 == 1);
                        [~, pidx] = max(p2(1,:));

                        if not(isempty(internalP2)),
                            if length(internalP2) == 2 && abs(internalP2(1)-pelvics(3))<0.001 && abs(internalP2(2)-pelvics(4))<0.001 && size(p2,2)==2,
                                dummyx = flipud(xdata(idx+1:length(xdata)));
                                dummyy = flipud(ydata(idx+1:length(ydata)));
                                idxPelvic2start = abs(pelvics(4)-length(xdata)-1);
                                idxPelvic2end = abs(pelvics(3)-length(xdata)-1);

                                lowerx = [dummyx(1:find(dummyx(1:idxPelvic2start-1)<p2(1,1),1,'last')); p2(1,1); dummyx(idxPelvic2start+1:idxPelvic2end-1); p2(1,2); dummyx(find(dummyx(1:idxPelvic2start-1)>p2(1,2),1,'first'):idxPelvic2start-1); dummyx(idxPelvic2end+1:end)];
                                lowery = [dummyy(1:find(dummyx(1:idxPelvic2start-1)<p2(1,1),1,'last')); p2(2,1); dummyy(idxPelvic2start+1:idxPelvic2end-1); p2(2,2); dummyy(find(dummyx(1:idxPelvic2start-1)>p2(1,2),1,'first'):idxPelvic2start-1); dummyy(idxPelvic2end+1:end)];
                                
                                lowerx = flipud(lowerx);
                                lowery = flipud(lowery);
                                
                                xdata = [upperx; lowerx];
                                ydata = [uppery; lowery];
                            elseif isempty(internalB2),
                                if xdata(pelvics(4)+1) < xdata(pelvics(4)),
                                    firstval = pelvics(4)+1;
                                    for val = firstval:1:length(xdata),
                                        %Check to see if the line segment between
                                        %current outline point and next outline point
                                        %contains the posteriormost cross-over point.
                                        %Test is distance between segment and point.
                                        test = abs(det([[xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)],p2(:,pidx)-[xdata(val);ydata(val)]]))/norm([xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)]);
                                        %If segment doesn't contain 3rd
                                        %crossover, keep searching
                                        if test > 0.001,
                                            continue
                                        %If segment does contain crossover,
                                        %mark that we've hit the crossover and
                                        %stop searching
                                        else
                                            idxcontoutline = val+1;
                                            break
                                        end
                                    end
                                    xdata = [upperx; xdata(idx+1:internalP2(1)-1); p2(1,pidx); xdata(idxcontoutline:end)];
                                    ydata = [uppery; ydata(idx+1:internalP2(1)-1); p2(2,pidx); ydata(idxcontoutline:end)];    
                                else
                                    xdata = [upperx; xdata(idx+1:internalP2(1)-1); p2(1,pidx); xdata(pelvics(4)+1:end)];
                                    ydata = [uppery; ydata(idx+1:internalP2(1)-1); p2(2,pidx); ydata(pelvics(4)+1:end)];
                                end
                            else
                                if size(p2,2) < 3 && xdata(pelvics(3))>xdata(pelvics(4)+1),
                                    internalB2 = internalB2 + pelvics(4);

                                    xdata = [upperx; xdata(idx+1:internalP2(1)-1); p2(1,pidx); xdata(internalB2(end)+1:end)];
                                    ydata = [uppery; ydata(idx+1:internalP2(1)-1); p2(2,pidx); ydata(internalB2(end)+1:end)];
                                elseif size(p2,2) == 4,
                                    %Find last body point before first
                                    %crossover
                                    dummyx = flipud(xdata(idx+1:length(xdata)));
                                    dummyy = flipud(ydata(idx+1:length(ydata)));
                                    idxPelvic2start = abs(pelvics(4)-length(xdata)-1);
                                    idxPelvic2end = abs(pelvics(3)-length(xdata)-1);

                                    lastbody = find(dummyx(1:idxPelvic2start-1) > p2(1,2),1,'first')-1;
                                    
                                    abovelowerlim = dummyx(idxPelvic2start+1:idxPelvic2end) > p2(1,2);
                                    belowupperlim = dummyx(idxPelvic2start+1:idxPelvic2end) < p2(1,3);
                                    both = idxPelvic2start+find(abovelowerlim+belowupperlim == 2);
                                    
                                    firstval = lastbody+1;
                                    for val = lastbody:1:idxPelvic2end-1,
                                        %Check to see if current point is an
                                        %internal point.  If so, skip.
                                        internalB2scaled=internalB2+pelvics(4);
                                        if any(abs((abs(internalP2-length(xdata)-1)-val)<0.01)) || any(abs(xdata(internalB2scaled)-dummyx(firstval))<0.01),
                                            firstval = val+1;
                                            continue
                                        end
                                        %Check to see if the line segment between
                                        %current outline point and next outline point
                                        %contains the posteriormost cross-over point.
                                        %Test is distance between segment and point.
                                        test = abs(det([[dummyx(val+1);dummyy(val+1)]-[dummyx(val);dummyy(val)],p2(:,pidx)-[dummyx(val);dummyy(val)]]))/norm([dummyx(val+1);dummyy(val+1)]-[dummyx(val);dummyy(val)]);
                                        %If segment doesn't contain 3rd
                                        %crossover, keep searching
                                        if test > 0.001,
                                            continue
                                        %If segment does contain crossover,
                                        %mark that we've hit the crossover and
                                        %stop searching
                                        else
                                            lastonpeak = val;
                                            break
                                        end
                                    end

                                    %Get the index where we start following the
                                    %outline normally again
                                    idxcontoutline = find(dummyx(idxPelvic2end:length(dummyx))>p2(1,pidx),1,'first') + idxPelvic2end-1;

                                    for val=idxPelvic2start:1:idxcontoutline-1,
                                        if dummyx(val) < dummyx(idxcontoutline) && dummyx(val) > p2(1,pidx),
                                            idxcontoutline = val;
                                            break
                                        end
                                    end


                                    lowerx = [dummyx(1:lastbody); p2(1,2); dummyx(both); p2(1,3); dummyx(firstval:lastonpeak); p2(1,pidx); dummyx(idxcontoutline:length(dummyx))];
                                    lowery = [dummyy(1:lastbody); p2(2,2); dummyy(both); p2(2,3); dummyy(firstval:lastonpeak); p2(2,pidx); dummyy(idxcontoutline:length(dummyy))];
    %                                 plot(lowerx,lowery,'ro')
                                    xdata = [upperx; flipud(lowerx)];
                                    ydata = [uppery; flipud(lowery)];
                                    
                                    
                                else

                                    %Find last body point before first
                                    %crossover
                                    dummyx = flipud(xdata(idx+1:length(xdata)));
                                    dummyy = flipud(ydata(idx+1:length(ydata)));
                                    idxPelvic2start = abs(pelvics(4)-length(xdata)-1);
                                    idxPelvic2end = abs(pelvics(3)-length(xdata)-1);

                                    [~,minidx] = min(p2(1,:));
                                    lastbody = find(dummyx(1:idxPelvic2start-1) < p2(1,minidx),1,'last');

                                    abovelowerlim = dummyx(idxPelvic2start:idxPelvic2end) > p2(1,minidx);
                                    belowupperlim = dummyx(idxPelvic2start:idxPelvic2end) < p2(1,2);
                                    both = idxPelvic2start-1 + find(abovelowerlim+belowupperlim == 2);

                                    firstval = lastbody+1;
                                    for val = lastbody:1:idxPelvic2end-1,
                                        %Check to see if current point is an
                                        %internal point.  If so, skip.
                                        internalB2scaled=internalB2+pelvics(4);
                                        if any(abs((abs(internalP2-length(xdata)-1)-val)<0.01)) || any(abs(xdata(internalB2scaled)-dummyx(firstval))<0.01),
                                            firstval = val+1;
                                            continue
                                        end
                                        %Check to see if the line segment between
                                        %current outline point and next outline point
                                        %contains the posteriormost cross-over point.
                                        %Test is distance between segment and point.
                                        test = abs(det([[dummyx(val+1);dummyy(val+1)]-[dummyx(val);dummyy(val)],p2(:,pidx)-[dummyx(val);dummyy(val)]]))/norm([dummyx(val+1);dummyy(val+1)]-[dummyx(val);dummyy(val)]);
                                        %If segment doesn't contain 3rd
                                        %crossover, keep searching
                                        if test > 0.001,
                                            continue
                                        %If segment does contain crossover,
                                        %mark that we've hit the crossover and
                                        %stop searching
                                        else
                                            lastonpeak = val;
                                            break
                                        end
                                    end

                                    %Get the index where we start following the
                                    %outline normally again
                                    idxcontoutline = find(dummyx(idxPelvic2end:length(dummyx))>p2(1,pidx),1,'first') + idxPelvic2end-1;

                                    for val=idxPelvic2start:1:idxcontoutline-1,
                                        if dummyx(val) < dummyx(idxcontoutline) && dummyx(val) > p2(1,pidx),
                                            idxcontoutline = val;
                                            break
                                        end
                                    end


                                    lowerx = [dummyx(1:lastbody); p2(1,minidx); dummyx(both); p2(1,2); dummyx(firstval:lastonpeak); p2(1,pidx); dummyx(idxcontoutline:length(dummyx))];
                                    lowery = [dummyy(1:lastbody); p2(2,minidx); dummyy(both); p2(2,2); dummyy(firstval:lastonpeak); p2(2,pidx); dummyy(idxcontoutline:length(dummyy))];
%                                     plot(lowerx,lowery,'ro')
                                    xdata = [upperx; flipud(lowerx)];
                                    ydata = [uppery; flipud(lowery)];

                                end
                            end

                        else
                            trigger = 0;
                            for val = 1:1:length(xdata(pelvics(3):end-1)),
                                val = val + pelvics(3)-1;
                                test = abs(det([[xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)],p2(:,pidx)-[xdata(val);ydata(val)]]))/norm([xdata(val+1);ydata(val+1)]-[xdata(val);ydata(val)]);
                                if test <0.001 && trigger == 0,
                                    start = val;
                                    trigger = 1;
                                elseif test<0.001 && trigger == 1,
                                    stop = val+1;
                                    break
                                end
                            end

                            xdata = [upperx; xdata(idx+1:start); p2(1,pidx); xdata(stop:end)];
                            ydata = [uppery; ydata(idx+1:start); p2(2,pidx); ydata(stop:end)];

                        end
                    else
                        internalB2 = inpolygon(xdata(pelvics(4)+1:end),ydata(pelvics(4)+1:end),xdata(find(xdata(1:pelvics(3))<xdata(pelvics(4)),1,'last')+1:pelvics(4)),ydata(find(xdata(1:pelvics(3))<xdata(pelvics(4)),1,'last')+1:pelvics(4)));
                        internalB2 = find(internalB2 == 1);
                        if not(isempty(internalB2)),
                            xdata = [upperx; xdata(idx+1:pelvics(4)); xdata(internalB2(length(internalB2)):end)];
                            ydata = [uppery; ydata(idx+1:pelvics(4)); ydata(internalB2(length(internalB2)):end)];
                        else
                            xdata = [upperx; xdata(idx+1:end)];
                            ydata = [uppery; ydata(idx+1:end)];
                        end

                    end
                else
                    xdata = [upperx; xdata(idx+1:end)];
                    ydata = [uppery; ydata(idx+1:end)];
                end

            end

%             %Use these plot features to verify that the pelvic/pectoral
%             %overlap has been handled properly.  Turn off for final running
%             %of code.
%             hold on
%             if exist('pelvicFin1','var'),
%                 plot(pelvicFin1(:,1),pelvicFin1(:,2),'r.')
%             end
%             if exist('pelvicFin2','var'),
%                 plot(pelvicFin2(:,1),pelvicFin2(:,2),'r.')
%             end
% %             % plot(xdata1,ydata1,'mo')
% %             % plot(xdata1,ydata1,'g--')
%             plot(xdata,ydata,'mo')
%             plot(xdata,ydata,'g--')
% %             plot(p1(1),p1(2),'g*')
% %             plot(upperx,uppery,'mo')
% %             plot(p2(1),p2(2),'g*')
%             return




            figure(1)
            imshow(im,'XData',[xmin xmax],'YData',[ymin ymax])    % display current movie frame, scaled to coordinates in pressure file
            %axis equal
            axis on
            set(gca,'YDir','normal')                                                    %set axes to correct direction (because they are flipped in images)
            hold on                                                                      %should see correct axes here
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
            %pause(0.5)  


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

