%Apr 18, 2019


%Reads in pressure fields, body depth information, fish outlines, and movie
%for each fish sequence, makes nice plots/animations of pressure fields,
%calculates forces acting on the fish, and saves out a bunch of different
%types of force, outline, etc information for further analysis or plotting.

%Will not calculate forces on pelvic fins

clear all

close all


%Set metadata: filepath information, individual fish to analyze, scale factors,
%location of axis origins in DaVis, numberformat in the filenames,
%delimiters, and increment between frames used for pressure code.

individuals = {'klbg7','klbg3'};
mm2pixScale = 0.226674; %0.254432;
DaVis_x0 = 128; %110;
DaVis_y0 = 126; %881;
xmin = -DaVis_x0*mm2pixScale;
xmax = xmin+(1024*mm2pixScale);
ymin = (-1024+DaVis_y0)*mm2pixScale;
ymax = ymin+(1024*mm2pixScale);

outlineBase = 'H:/fish_outlines_v2/';
pressureBase = 'F:/fish pressure_v2/bluegill/';
depthBase = 'H:/fish body depth/';


dt = 0.001;                                                  %time step used in the PIV

numformat = '%05d';                                         % format of numbers in file name. '%05d' is five digit format with leading zeros

delimiter = '\t';                                             % delimiter between columns in PIV file
headerlines = 1;   

%Get ColorBrewer palette.  Note that cbrewer needs to be loaded first.
mycolorsStepCp=cbrewer('div','RdBu',29);
%Flips palette so that red = high values
fci = 0;
for row = 0:1:28
    fci = fci+1;
    flippedcolorsStepCp(fci,:) = mycolorsStepCp(29-row,:);
end
flippedcolorsStepCp = [flippedcolorsStepCp(1:13,:); ones(4,3)*0.9686; flippedcolorsStepCp(17:29,:)];
% return



%loop through each individual which will be analyzed.
for i = 1:1:length(individuals),
    
    %Report current individual - progress marker for while code is running
    individuals{i}
    
    %Find all sequences for that individual
    d = dir([outlineBase individuals{1,i}]);;
    isub = [d(:).isdir];
    folders = {d(isub).name};
    folders(ismember(folders,{'.','..'}))=[];
    
    %Get corresponding flow speed
    if strcmp(individuals{i},'klbg2'),
        flowBase = 0.1;
        latArea = 0.00229634;
    elseif strcmp(individuals{i},'klbg4'),
        flowBase = 0.097;
        latArea = 0.002140139;
    elseif strcmp(individuals{i},'klbg5'),
        flowBase = 0.096;
        latArea = 0.002082938;
    elseif strcmp(individuals{i},'klbg7'),
        flowBase = 0.093;
        latArea = 0.001834378;
    elseif strcmp(individuals{i},'klbg3'),
        flowBase = 0.115;
        latArea = 0.002814104;
    end

    
    %Get scale
    scale = mm2pixScale;
    
    
    %Make body outline path
    bodyoutlinepath = [depthBase individuals{i} '_body.xlsx'];
    %Find body length/depth and pectoral fin length/depth info.  Note
    %height_dict.m must be loaded first.
    [bodylength, bodydepth] = height_dict_body_only(bodyoutlinepath);
    
    %For each sequence (saved to unique folders),
    for f = 1:1:size(folders,2),
        
        %Report which sequence we're on - for tracking code progress
        folders{f}
        
        %Extract the first and last frame number from the filename
        frameStr = strsplit(folders{f},'frm');
        frameStr = strsplit(frameStr{2},'-');
        first = str2num(frameStr{1});
        last = str2num(frameStr{2});
        
        %Extract flow speed from filename and convert to meters/sec
        flowBL = strsplit(folders{f},'BLs');
        flowBL = strsplit(flowBL{1},'_');
        flow = str2num(flowBL{2})*flowBase;
        
        
        %set filepath and name stub for pressure fields and outlines
        pressureStub = [pressureBase individuals{i} '/' folders{f} '/queen2_dT10ms_'];
        outlineStub = [outlineBase individuals{i} '/' folders{f} '/' folders{f} '_'];
        
        vidfr=1;
        
        %Find all outline files for this trial, and use number of files to
        %determine how many frames we have to analyze
        d2 = dir([outlineBase individuals{i} '/' folders{f}]);
        isub2 = [d2(:).isdir];
        isub2 = abs(1-isub2);
        numFrames = sum(isub2)-1;
        
        %Set a count index
        count = 1;
        %Initialize storage for total x & y forces, positive x & y forces, and
        %negative x & y forces
        totx = zeros(numFrames,1);
        toty = zeros(numFrames,1);
        totposx = zeros(numFrames,1);
        totposy = zeros(numFrames,1);
        totnegx = zeros(numFrames,1);
        totnegy = zeros(numFrames,1);
        
        
        %Initialize storage for midlines
        midlineSaveX = [];
        midlineSaveY = [];
        
        %Clear variables from previous trial
        clear pressures
        clear boundariesx
        clear boundariesy
        clear forcexall
        clear forceyall
        clear normsx
        clear normsy
        clear locs
        clear dists

        
        %for every 10th frame of the sequence,
        for filenum = first:10:last-11,
            %clear variables from last frame
            close all

            clear forcex
            clear forcey
            

            clear pelvics
            clear pelvicFin1
            clear pelvicFin2
            clear p1
            clear p2

            
            %Display current filenum so we can track code progress
            filenum
            
            %Make sure folder where will save pressure & force results exists, and if not,
            %create it
            savefolder = [pressureBase individuals{i} '/' folders{f}];
            if ~exist(savefolder,'dir'),
                mkdir(savefolder);
            end
            
            
            %Read in pressure files and figure out x and y dimensions
            pressure = dlmread([pressureStub num2str(filenum,numformat) '.dat'],',',0,0);
            pressXDim = length(find(pressure(:,2) == pressure(1,2)));
            pressYDim = length(find(pressure(:,1) == pressure(1,1)));

            %Read in fish outline
            interface = xlsread([outlineStub num2str(filenum,numformat) '.xlsx'],'Outline');
            %Move origin in outline files to match origin in pressure files
            xoutline = interface(:,2) + xmin;
            youtline = interface(:,4) + ymin;
            outline = [xoutline youtline];
            

% return
            %Check for pelvic fins
            [status, sheets] = xlsfinfo([outlineStub num2str(filenum,numformat) '.xlsx']);
            check = find(strcmp(sheets, 'Pelvic Fins'));
            if not(isempty(check)),
                pelvics = xlsread([outlineStub num2str(filenum,numformat) '.xlsx'],'Pelvic Fins');
            end

            if exist('pelvics','var'),
                %Find the coordinates for upper(1) and lower(2) pelvic fins
                if pelvics(1) ~= 0,
                    pelvicFin1 = outline(pelvics(1):pelvics(2),:);

                    %Check to see if the pelvic fins overlap with pectoral fins
                    %(by seeing if the arc of the pelvic fin crosses the rest
                    %of the body/fins outline).  Repeat for upper and lower
                    %fins.
                    p1 = InterX(pelvicFin1', [outline(1:pelvics(1),1) outline(1:pelvics(1),2)]');
                    %Make sure the anteriormost point of the pelvic fin isn't marked
                    %as an overlap
                    rep = find(abs(p1(1,:) - pelvicFin1(1,1))<0.0001);
                    if not(isempty(rep)),
                        p1(:,rep) = [];
                    end
                %If pelvic fin isn't out, make pelvicFin1 an empty matrix
                else
                    pelvicFin1 = [];
                end
                if pelvics(3) ~= 0,
                    pelvicFin2 = outline(pelvics(3):pelvics(4),:);

                    p2 = InterX(pelvicFin2', [outline(pelvics(4):end,1) outline(pelvics(4):end,2)]');
                    %Make sure the anteriormost point of the pelvic fin isn't marked
                    %as an overlap
                    rep = find(abs(p2(1,:) - pelvicFin2(end,1))<0.0001);
                    if not(isempty(rep)),
                        p2(:,rep) = [];
                    end
                %If pelvic fin isn't out, make pelvicFin2 an empty matrix    
                else
                    pelvicFin2 = [];
                end
%                 return
                %find index of tail tip (used to split fish outline into
                %upper and lower halves later)
                [~,idx] = max(outline(:,1));

                %If we don't have an upper pelvic fin, keep whole upper
                %outline
                if not(exist('p1', 'var')),
                    outlineNoPelvics = outline(1:idx,:);
                %Or, if no intersections between upper pelvic fin and body were 
                %found, or exactly two intersections are found
                elseif isempty(p1) || size(p1,2) == 2,
                    %Drop pelvic fin and keep rest of outline
                    if outline(pelvics(1),1) > outline(pelvics(1)-1,1),
                        outlineNoPelvics = [outline(1:pelvics(1),:); outline(pelvics(2):idx,:)];
                    elseif outline(pelvics(1),1) < outline(pelvics(1)-1,1) && outline(pelvics(2),1) < outline(pelvics(1)-1,1),
                        outlineNoPelvics = [outline(1:pelvics(1)-1,:); outline(pelvics(2)+1:idx,:)];
                    else
                        outlineNoPelvics = [outline(1:pelvics(1)-1,:); outline(pelvics(2):idx,:)];
                    end
                %Or, if exactly one intersection is found (means posterior end
                %of pectoral fin overlaps with anterior end of pelvic),
                elseif size(p1,2) == 1,
                    %Determine if the last point on the pelvic fin and the
                    %first point on the body have nearly the same y value
                    dif = abs(outline(pelvics(1)-1,2) - outline(pelvics(1),2));
                    %If the y-values are similar,
                    if dif < 1,
                        %Delete the whole pelvic fin including and the previous point
                        %(previous will be an internal point)
                        outlineNoPelvics = [outline(1:pelvics(1)-1,:); outline(pelvics(2):idx,:)];
                    %Delete the pelvic fin
                    else
                        outlineNoPelvics = [outline(1:pelvics(1),:); outline(pelvics(2):idx,:)];
                    end
                elseif size(p1,2) > 2,
                    %Drop pelvic fin and keep rest of outline
                    outlineNoPelvics = [outline(1:pelvics(1)-1,:); outline(pelvics(2):idx,:)];
                %Keep the whole upper outline 
                else
                    outlineNoPelvics = outline(1:idx,:);
                end

                %If we don't have a lower pelvic fin, keep whole lower
                %outline
                if not(exist('p2','var')),
                    outlineNoPelvics = [outlineNoPelvics; outline(idx+1:end,:)];
                %Or, if no intersections between lower pelvic fin and body were 
                %found, or exactly two intersections are found
                elseif isempty(p2) || size(p2,2) == 2,
                    %Drop pelvic fin and keep rest of outline, if pectoral
                    %is larger
                    if outline(pelvics(4),1) > outline(pelvics(4)+1,1),
                        outlineNoPelvics = [outlineNoPelvics; outline(idx+1:pelvics(3),:); outline(pelvics(4):end,:)];
                    elseif outline(pelvics(4),1) < outline(pelvics(4)+1,1) && outline(pelvics(3),1) < outline(pelvics(4)+1,1),
                        outlineNoPelvics = [outlineNoPelvics; outline(idx+1:pelvics(3)-1,:); outline(pelvics(4)+1:end,:)];
                    else
                        %If pelvic fin ends more anteriorly than pectoral,
                        %drop the last point of the pelvic too
                        outlineNoPelvics = [outlineNoPelvics; outline(idx+1:pelvics(3),:); outline(pelvics(4)+1:end,:)];
                    end
                %Or, if exactly one intersection is found (means posterior end
                %of pectoral fin overlaps with anterior end of pelvic),
                elseif size(p2,2) == 1,
                    %Determine if the last point on the pelvic fin and the
                    %first point on the body have nearly the same y value
                    dif = abs(outline(pelvics(4)+1,2) - outline(pelvics(4),2));
                    %If the y-values are similar,
                    if dif < 1,
                        %Delete the whole pelvic fin including and the previous point
                        %(previous will be an internal point)
                        outlineNoPelvics = [outlineNoPelvics; outline(idx+1:pelvics(3),:); outline(pelvics(4)+1:end,:)];
                    %Delete the pelvic fin
                    else
                        outlineNoPelvics = [outlineNoPelvics; outline(idx+1:pelvics(3),:); outline(pelvics(4):end,:)];
                    end
                %Or, if more than two intersections are found (pectoral passes through
                %a LARGE pelvic fin),    
                elseif size(p2,2) > 2,
                    %Drop pelvic fin and keep rest of outline
                    outlineNoPelvics = [outlineNoPelvics; outline(idx+1:pelvics(3),:); outline(pelvics(4)+1:end,:)];
                %Keep the whole lower outline
                else
                    outlineNoPelvics = [outlineNoPelvics; outline(idx+1:end,:)];
                end
            %If we don't have pelvic fins, just keep outline the way it is
            else
                outlineNoPelvics = outline;
            end

            %Find the index of the tail tip - for splitting the outline
            %into upper and lower halves
            [~, oNPidx] = max(outlineNoPelvics(:,1));
            
            
            outline_forMid = [outlineNoPelvics(:,1)-xmin outlineNoPelvics(:,2)-ymin];
% return
            %Use the pelvic-fin-free outline to get the fish's midline.  Note: fish_midline_fctn.m must
            %be loaded
            midline = fish_midline_fctn(filenum, outline_forMid, scale, numformat, delimiter, headerlines, DaVis_x0, DaVis_y0);

% return
            %Interpolate more points onto the midline using arc-length
            %interpolation, and yielding equidistant points
            midlineXtra = interparc(linspace(0,1,100),midline(:,1),midline(:,2),'pchip');
            midlineXtra = midlineXtra/1000;
            %Split midline data into x and y coordinate lists
            midlineXtraX = midlineXtra(:,1);
            midlineXtraY = midlineXtra(:,2);
            
            %Add current midline to the master array of midlines
            midlineSaveX = [midlineSaveX midlineXtraX];
            midlineSaveY = [midlineSaveY midlineXtraY];
            

            %For each point, get a distance traveled from the nose measure
            distOnMid = zeros(100,1);
            for pt = 2:1:100,
                distOnMid(pt) = distOnMid(pt-1)+sqrt((midlineXtraX(pt-1)-midlineXtraX(pt))^2+(midlineXtraY(pt-1)-midlineXtraY(pt))^2);
            end
            
%             plot(midline(:,1),midline(:,2))
%             hold on
%             plot(midlineXtraX*1000,midlineXtraY*1000,'k.')
%             return
            



            %Interpolate on the arc described by the upper body outline with
            %evenly spaced points
            upperXtra = interparc(linspace(0,1,100),outlineNoPelvics(1:oNPidx,1),outlineNoPelvics(1:oNPidx,2),'pchip');


%     return

            
            lowerhalf = flipud([outlineNoPelvics(oNPidx:end,:); outlineNoPelvics(1,:)]);
            lowerXtra = interparc(linspace(0,1,100),lowerhalf(:,1),lowerhalf(:,2),'pchip');


            %Assemble upper and lower outlines into one complete outline
            outlineXtra = [upperXtra; flipud(lowerXtra(2:end-1,:))];
            
            %Convert outline into pixels
            outlineXtraPix = (outlineXtra(:,1)-xmin)/scale;
            outlineXtraPix = [outlineXtraPix (outlineXtra(:,2)-ymin)/scale];

            %Create a mask from the interpolated outline
            imOutline = poly2mask(outlineXtraPix(:,1),outlineXtraPix(:,2),1024,1024);
%             imshow(imOutline)
%             axis on
%             return

            
            

            
%             %Plot original and interpolated outlines to verify the above
%             %worked
%             plot(outline(:,1),outline(:,2))
%             hold on
%             plot(outlineNoPelvics(:,1),outlineNoPelvics(:,2),'go')
%             plot(outlineXtra(:,1),outlineXtra(:,2),'b.')
%             plot(midline(:,1),midline(:,2))
%             return
% % pause(1.5)
% % continue
            



            
            %Inflate body outline via snake
            
            



            %for checking kinematics - turn on plotter and video save features

      
            %Make binary image of body
            xpix = (outlineXtra(:,1)-xmin)/scale;
            ypix = (outlineXtra(:,2)-ymin)/scale;
            imBody = poly2mask(xpix,ypix,1024,1024);
%             figure(2)
%             imshow(imBody)
%             axis on
%             hold on
%             set(gca,'YDir','normal')            %make sure y-axis is oriented correctly (positive up)
%             return

            %get distances away from boundary
            Daway = bwdist(imBody);

            imBody = imcomplement(imBody);


            %get distance from boundary
            D = -bwdist(imBody);
            D = D + Daway;
            [rows cols] = size(D);
            [xdim ydim] = meshgrid(1:1:cols,1:1:rows);

            iterations = 5;

            %alpha = elasticity weight, beta = flexibility weight, gamma = step size for moving x and
            %y, and kappa is weight of external gradient
            numpts = size(outlineXtra,1);
            alpha = 0.1;
            beta = 0.2;
            gamma = 0.4;
            kappa = 0.5;

            xs = xpix;
            ys = ypix;

            % return

            for iter = 1:1:iterations,
                %set up elastic force on each point - force keeping points 
                %evenly spread out.  Basically, saying here that we want the
                %average distance between points to be about the same.
                F_elasticx = alpha*2*((circshift(xs,1)-xs)+(circshift(xs,-1)-xs));
                F_elasticy = alpha*2*((circshift(ys,1)-ys)+(circshift(xs,-1)-xs));
  

                %set up stiffness force on each point - force keeping points 
                %from making too bendy a line.  Only works >2 pts from end.
                %will do this by preventing curvature (1/radcurv) from being
                %too big.  Like above, we'll do this by saying we want the
                %average curvature between adjacent points to not be very 
                %different.  This means we need the curvature at j-1, j, and
                %j+1.
                
                %curvature 
                m1 = (ys-circshift(ys,1))./(xs-circshift(xs,1));
                m2 = (circshift(ys,-1)-ys)./(circshift(xs,-1)-xs);
                xc = (m1.*m2.*(circshift(ys,1)-circshift(ys,-1))+m2.*(circshift(xs,1)+xs)-m1.*(xs+circshift(xs,-1)))./(2*(m2-m1));
                yc = -(1./m1).*(xc-(circshift(xs,1)+xs)/2)+(circshift(ys,1)+ys)/2;
                rad = sqrt((xs-xc).^2+(ys-yc).^2);

                F_stiff = beta*2*((1./circshift(rad,1)-1./rad)+(1./circshift(rad,-1)-1./rad));
                    

                %Last, set up force to maximize distance from the boundary.
                %If distance at x+1 is bigger, move that way, etc.

                %get distance value at x+1 and x-1
                Dlf = interp2(xdim,ydim,D,xs-2,ys);
                Drt = interp2(xdim,ydim,D,xs+2,ys);
                %get distance value at y+1 and y-1
                Dup = interp2(xdim,ydim,D,xs,ys-2);
                Ddn = interp2(xdim,ydim,D,xs,ys+2);

                %get forces
                F_distx = (kappa/2)*(Drt-Dlf);
                F_disty = (kappa/2)*(Ddn-Dup);



                xs = xs + gamma*(nansum([F_distx,F_stiff,F_elasticx],2));
                ys = ys + gamma*(nansum([F_disty,F_stiff,F_elasticy],2));


                
%                 if iter == 1,
%                     figure(3)
%                 end
%                 plot(xpix,ypix,'c-','LineWidth',2)
%                 hold on
%                 plot(xs,ys,'g*')
% %                 axis([0 1024 0 1024])
%                 axis square
%                 hold off
%                 pause(0.01)
% %                 return
                
            end
            
            boundx = xs*scale+xmin;
            boundy = ys*scale+ymin;
            boundx = double(boundx);
            boundy = double(boundy);
            
%             figure(1)
%             hold on
%             plot(boundx,boundy)
%             plot(boundx,boundy,'.')
% %             axis([xmin xmax ymin ymax])
%             axis equal
% %             return




            
            
            
            %Convert boundary coordinates into meters
            boundx = boundx/1000;
            boundy = boundy/1000;

            surfdx = boundx - circshift(boundx,1);                                  % compute x-spacing between boundary coordinates
            surfdy = boundy - circshift(boundy,1);                                  % compute y-spacing between boundary coordinates

            surfposx = (boundx + circshift(boundx,1))./2;                           % compute x-midpoint of each object surface facet
            surfposy = (boundy + circshift(boundy,1))./2;                           % compute y-midpoint of each object surface facet
            
            surfnormx = surfdy;                                                     % compute x-component of surface normal
            surfnormy = -surfdx;                                                    % compute y-component of surface normal

            surfunitnormx = surfnormx./(sqrt(surfnormx.^2 + surfnormy.^2));         % compute x-component of unit surface normal
            surfunitnormy = surfnormy./(sqrt(surfnormx.^2 + surfnormy.^2));         % compute y-component of unit surface normal
% return
            
            
            
%             
%             figure(5)
%             plot(boundx,boundy)
%             hold on

            
            
            
            
            %Assemble the final outline from the body & fin outlines (because
            %some points will be internal)
            %Start with body outline
            assembled_bound = [boundx boundy];
            %Create a dictionary in a cell structure, labeling each outline
            %point as belonging to the body
            dict = {assembled_bound};
            for p = 1:1:size(assembled_bound,1),
                dict{1,2}{p,1} = 'body';
            end
            
            
            
            %If pelvic fins are extended,
            if exist('pelvics','var'),
                %If we have an upper pelvic fin,
                if ~isempty(pelvicFin1),
                    %Find the points in the outline that are internal to
                    %the pelvic fin outline
                    dropOut = inpolygon(assembled_bound(:,1),assembled_bound(:,2),pelvicFin1(:,1)/1000,pelvicFin1(:,2)/1000);

                    %Get the indices of the first and last point to delete
                    %from the body outline
                    startCutOut = find(dropOut==1,1);
                    stopCutOut = find(dropOut==1,1,'last');
                    
                    %if we have start and stop indices,
                    if ~isempty(startCutOut) && ~isempty(stopCutOut),
                        %drop the points internal to the pelvic fin
                        assembled_bound = [assembled_bound(1:startCutOut-1,:); assembled_bound(stopCutOut+1:end,:)];
                        %Drop the associated points in boundary spacing and
                        %segment midpoints
                        surfdx = [surfdx(1:startCutOut-1); surfdx(stopCutOut+1:end)];
                        surfdy = [surfdy(1:startCutOut-1); surfdy(stopCutOut+1:end)];
                        surfposx = [surfposx(1:startCutOut-1); surfposx(stopCutOut+1:end)];
                        surfposy = [surfposy(1:startCutOut-1); surfposy(stopCutOut+1:end)];
                        surfnormx = [surfnormx(1:startCutOut-1); surfnormx(stopCutOut+1:end)];
                        surfnormy = [surfnormy(1:startCutOut-1); surfnormy(stopCutOut+1:end)];
                        surfunitnormx = [surfunitnormx(1:startCutOut-1); surfunitnormx(stopCutOut+1:end)];
                        surfunitnormy = [surfunitnormy(1:startCutOut-1); surfunitnormy(stopCutOut+1:end)];

                    end
                end

                %If we have a lower pelvic fin,
                if ~isempty(pelvicFin2),
                    %Find the points in the outline that are internal to
                    %the pelvic fin outline
                    dropOut = inpolygon(assembled_bound(:,1),assembled_bound(:,2),pelvicFin2(:,1)/1000,pelvicFin2(:,2)/1000);

                    %Get the indices of the first and last point to delete
                    %from the body outline
                    startCutOut = find(dropOut==1,1);
                    stopCutOut = find(dropOut==1,1,'last');
                    
                    %if we have start and stop indices,
                    if ~isempty(startCutOut) && ~isempty(stopCutOut),
                        %drop the points internal to the pelvic fin
                        assembled_bound = [assembled_bound(1:startCutOut-1,:); assembled_bound(stopCutOut+1:end,:)];
                        %Drop the associated points in boundary spacing and
                        %segment midpoints
                        surfdx = [surfdx(1:startCutOut-1); surfdx(stopCutOut+1:end)];
                        surfdy = [surfdy(1:startCutOut-1); surfdy(stopCutOut+1:end)];
                        surfposx = [surfposx(1:startCutOut-1); surfposx(stopCutOut+1:end)];
                        surfposy = [surfposy(1:startCutOut-1); surfposy(stopCutOut+1:end)];
                        surfnormx = [surfnormx(1:startCutOut-1); surfnormx(stopCutOut+1:end)];
                        surfnormy = [surfnormy(1:startCutOut-1); surfnormy(stopCutOut+1:end)];
                        surfunitnormx = [surfunitnormx(1:startCutOut-1); surfunitnormx(stopCutOut+1:end)];
                        surfunitnormy = [surfunitnormy(1:startCutOut-1); surfunitnormy(stopCutOut+1:end)];
                    end
                end
            end
            
            
            

%             
%             %plot the outline and boundary coordinates to make sure fins
%             %were dropped correctly
%             plot(outline(:,1),outline(:,2))
%             hold on
%             plot(boundx*1000,boundy*1000,'.')
%             plot(assembled_bound(:,1)*1000,assembled_bound(:,2)*1000)
%             return


            %GET HEIGHT ASSIGNMENTS%
            %Initialize a new entry area in the dictionary
            dict{1,3} = zeros(size(dict{1,1},1),1);
            dict{1,4} = zeros(size(dict{1,1},1),1);
            %For each point in the dictionary,
            for pt = 1:1:size(dict{1,1},1),
                %If labeled "body",
                if strcmp(dict{1,2}{pt,1},'body'),
                    %Get the distance along midline by finding the distance
                    %associated with the nearest midline point
                    distances = sqrt((midlineXtra(:,1)-dict{1,1}(pt,1)).^2+(midlineXtra(:,2)-dict{1,1}(pt,2)).^2);
                    [~, closest] = min(distances);
                    %Normalize by body length
                    distVal = distOnMid(closest)/max(distOnMid);
                    %Body length can't be less than 0 or more than 1, so if
                    %normalization does (happens because of dilation), manually
                    %set distances to 0 or 1.
                    if distVal > 1,
                        distVal = 1;
                    elseif distVal < 0,
                        distVal = 0;
                    end
                    %Interpolate the body length/depth data to get the body
                    %depth at the current point, and store in dictionary
                    dict{1,3}(pt,1) = interp1(bodylength,bodydepth,distVal);
                    dict{1,4}(pt,1) = distVal;
                end
            end
            %Initialize storage for the list of "heights" (body depths) for
            %each point in the assembled boundary
            heights = zeros(size(assembled_bound,1),1);
            location = {assembled_bound};
            distValues = zeros(size(assembled_bound,1),1);
            %Look up correct values in the dictionary
            for pt = 1:1:size(assembled_bound,1),
                [~, distIdx] = ismember(assembled_bound(pt,:), dict{1,1}, 'rows');
                heights(pt,1) = dict{1,3}(distIdx,1);
                location{1,2}{pt,1} = dict{1,2}(distIdx,1);
                distValues(pt,1) = dict{1,4}(distIdx,1);
            end


            surfpress = griddata(pressure(:,1),pressure(:,2),pressure(:,7),surfposx,surfposy,'nearest');       % compute pressure at midpoint of each surface facet
            time(count) = (count-1)*0.001;                                   %Time at each step


            %Calculate forces at each boundary point (area of rectangle with
            %dimensions = spacing*bodydepth)*(direction of normal
            %vector)*(pressure)/1000 to get units in Pa
            forcex(1:size(assembled_bound,1),1) = sqrt(surfdx.^2+surfdy.^2).*surfunitnormx.*surfpress.*heights/1000;
            forcey(1:size(assembled_bound,1),1) = sqrt(surfdx.^2+surfdy.^2).*surfunitnormy.*surfpress.*heights/1000;

            %Get rid of any really outrageously wrong vectors - 10 standard
            %deviations above mean
            xfix = find(abs(forcex) > 10*std(abs(forcex(~isnan(forcex))))+mean(abs(forcex(~isnan(forcex)))));
            yfix = find(abs(forcey) > 10*std(abs(forcey(~isnan(forcey))))+mean(abs(forcey(~isnan(forcey)))));
            if ~isempty(xfix),
                forcex(xfix) = NaN;
                forcey(xfix) = NaN;
            end
            if ~isempty(yfix),
                forcex(yfix) = NaN;
                forcey(yfix) = NaN;
            end

            %Sum up total positive and total negative x & y forces.
            totposx(count,1) = sum(forcex(forcex>0));
            totposy(count,1) = sum(forcey(forcey>0));
            totnegx(count,1) = sum(forcex(forcex<0));
            totnegy(count,1) = sum(forcey(forcey<0));

            %sum up total x & y forces
            totx(count,1) = nansum(forcex);
            toty(count,1) = nansum(forcey);
            
% return


            %Set up arrays of data for saving
            if filenum == first,
                pressures(:,count) = surfpress;
                boundariesx(:,count) = assembled_bound(:,1);
                boundariesy(:,count) = assembled_bound(:,2);
                forcexall(:,count) = forcex;
                forceyall(:,count) = forcey;
                normsx(:,count) = surfunitnormx;
                normsy(:,count) = surfunitnormy;
                locs(:,count) = location{1,2};
                dists(:,count) = distValues;
            elseif size(pressures,1) == size(surfpress,1),
                pressures(:,count) = surfpress;
                boundariesx(:,count) = assembled_bound(:,1);
                boundariesy(:,count) = assembled_bound(:,2);
                forcexall(:,count) = forcex;
                forceyall(:,count) = forcey;
                normsx(:,count) = surfunitnormx;
                normsy(:,count) = surfunitnormy;
                locs(:,count) = location{1,2};
                dists(:,count) = distValues;
            %Below, add NaN to bottom of list of points if the number of
            %points in the list is less than the number of rows in the data
            %array, or, add NaN to the bottom of all columns in the data
            %array if the number of points in the list is bigger than the
            %number of rows in the data array
            elseif size(pressures,1) < size(surfpress,1),
                pressures(size(pressures,1)+1:size(surfpress,1),:) = NaN;
                boundariesx(size(boundariesx,1)+1:size(assembled_bound,1),:) = NaN;
                boundariesy(size(boundariesy,1)+1:size(assembled_bound,1),:) = NaN;
                forcexall(size(forcexall,1)+1:size(forcex,1),:) = NaN;
                forceyall(size(forceyall,1)+1:size(forcey,1),:) = NaN;
                normsx(size(normsx,1)+1:size(surfunitnormx,1),:) = NaN;
                normsy(size(normsy,1)+1:size(surfunitnormy,1),:) = NaN;
                locs(size(locs,1)+1:size(location{1,2},1),:) = {'0'};
                dists(size(dists,1)+1:size(distValues,1),:) = NaN;
                pressures(:,count) = surfpress;
                boundariesx(:,count) = assembled_bound(:,1);
                boundariesy(:,count) = assembled_bound(:,2);
                forcexall(:,count) = forcex;
                forceyall(:,count) = forcey;
                normsx(:,count) = surfunitnormx;
                normsy(:,count) = surfunitnormy;
                locs(:,count) = location{1,2};
                dists(:,count) = distValues;
            elseif size(pressures,1) > size(surfpress,1),
                surfpress(size(surfpress,1)+1:size(pressures,1),:) = NaN;
                assembled_bound(size(assembled_bound,1)+1:size(boundariesx,1),:) = NaN;
                
                forcex(size(forcex,1)+1:size(forcexall,1),:) = NaN;
                forcey(size(forcey,1)+1:size(forceyall,1),:) = NaN;
                surfunitnormx(size(surfunitnormx,1)+1:size(normsx,1),:) = NaN;
                surfunitnormy(size(surfunitnormy,1)+1:size(normsy,1),:) = NaN;
                location{1,2}(size(location{1,2},1)+1:size(locs,1),:) = {'0'};
                distValues(size(distValues,1)+1:size(dists,1),:) = NaN;
                pressures(:,count) = surfpress;
                boundariesx(:,count) = assembled_bound(:,1);
                boundariesy(:,count) = assembled_bound(:,2);
                forcexall(:,count) = forcex;
                forceyall(:,count) = forcey;
                normsx(:,count) = surfunitnormx;
                normsy(:,count) = surfunitnormy;
                locs(:,count) = location{1,2};
                dists(:,count) = distValues;
            end

        figure('Color',[1 1 1]);
        axis on
        hold on      
        plot([outlineNoPelvics(:,1); outlineNoPelvics(1,1)],[outlineNoPelvics(:,2); outlineNoPelvics(1,2)],'LineWidth',1)
        plot(midline(:,1),midline(:,2),'r-','LineWidth',1)
        if filenum == first,
            pltxmin = midline(1,1)-5;
            pltxmax = midline(size(midline,1),1)+20;
            pltymin = midline(1,2)-30;
            pltymax = midline(1,2)+30;
        end
        axis([pltxmin pltxmax pltymin pltymax])
        daspect([1 1 1])
        axis off
        hold off

% return

        print([pressureBase individuals{i} '/' folders{f} '/kinematics_'  num2str(filenum,numformat)], '-dtiffn','-r300');

        

        


        



            %make fish silhouette
            xdata = xoutline;
            ydata = youtline;
            xdata(end+1) = xdata(1,1);
            ydata(end+1) = ydata(1,1);
            if exist('pelvics','var'),


                %Find the coordinates for upper(1) and lower(2) pelvic fins
                if pelvics(1) >0,
                    pelvicFin1 = interface(pelvics(1):pelvics(2),:);
                    pelvicFin1 = [pelvicFin1(:,2)+xmin pelvicFin1(:,4)+ymin];
                end
                if pelvics(3) > 0,
                    pelvicFin2 = interface(pelvics(3):pelvics(4),:);
                    pelvicFin2 = [pelvicFin2(:,2)+xmin pelvicFin2(:,4)+ymin];
                end
            end
                

                

            silhouette = [xdata ydata];

            %Make sure folder where will save silhouettes for figure-making exists, and if not,
            %create it
            silhouetteFolder = [pressureBase individuals{i} '/' folders{f} '/silhouettes for figures'];
            if ~exist(silhouetteFolder,'dir'),
                mkdir(silhouetteFolder);
            end
%             dlmwrite([silhouetteFolder '/silhouette_' num2str(filenum,numformat) '.dat'], silhouette, '\t');


            
            %Plot the pressure contour in the ColorBrewer palette.
            figure(4)

            contourf(reshape(pressure(:,1),pressYDim,pressXDim),reshape(pressure(:,2),pressYDim,pressXDim),reshape(pressure(:,7)/(0.5*1000*flow*flow),pressYDim,pressXDim),50,'LineColor','none')
            colormap(flippedcolorsStepCp)
            hold on
            fill(silhouette(:,1)/1000,silhouette(:,2)/1000,[0 0 0]);      %Plot fish silhouette from the fish outline data
            c = colorbar;

            caxis([-0.5 0.5])

            xlabel('position [m]')
            ylabel('position [m]')
            c.Label.String = 'C_{p}';

            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15 0.15 0.7 0.7]);

            axis equal

            daspect([1 1 1])


           
%             return

            %Save out figure
            print([pressureBase individuals{i} '/' folders{f} '/field_'  num2str(filenum,numformat)], '-dtiffn','-r300');
                       



            
            [~,assembled_max]=max(assembled_bound(:,1));
            %Plot the force vectors
            figure(5)
            %plot(assembled_bound(:,1),assembled_bound(:,2),'.','color',[0.651 0.3, 0.102])         %Plot force calculation boundary
            hold on
            fill(silhouette(:,1)/1000,silhouette(:,2)/1000,[0 0 0]);      %Plot fish silhouette from the fish outline data
%             return

            %Plot force vectors, with vector base on the boundary
            quiver(assembled_bound(1:assembled_max,1),assembled_bound(1:assembled_max,2),1000*0.1*forcex(1:assembled_max,1),1000*0.1*forcey(1:assembled_max,1),'AutoScale','off','color',[0.004 0.522 0.443]);
            quiver(assembled_bound(assembled_max+1:end,1),assembled_bound(assembled_max+1:end,2),1000*0.1*forcex(assembled_max+1:end,1),1000*0.1*forcey(assembled_max+1:end,1),'AutoScale','off','color',[0.651 0.3, 0.102]);
            
            


            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15 0.15 0.7 0.7]);
            axis([min(pressure(:,1))-0.02 max(pressure(:,1))+0.02  min(pressure(:,2))-0.02 max(pressure(:,2))+0.02])
%             axis equal

            %Make a scale bar (mN)
            line([0.05*range(pressure(:,1))+min(pressure(:,1)-0.02); 0.05*range(pressure(:,1))+min(pressure(:,1))-0.02+0.01],[0.9*range(pressure(:,2))+min(pressure(:,2)+0.02); 0.9*range(pressure(:,2))+min(pressure(:,2))+0.02],'color',[0 0 0],'LineWidth',1.5);
            str = '0.01 mN';
            text(0.05*range(pressure(:,1))+min(pressure(:,1))-0.02+0.013, 0.9*range(pressure(:,2))+min(pressure(:,2))+0.02,str,'FontSize',16);

            
            xlabel('position [m]')
            ylabel('position [m]')
            daspect([1 1 1])
            
            
            
% return
% % %     
            %Save out figure
            print([pressureBase individuals{i} '/' folders{f} '/forces_'  num2str(filenum,numformat)], '-dtiffn','-r300');


            count = count+1;

% return


        end
% return



        %Set an index equal to current frame number
        currFr = first;
        %Create a video file from the pressure field image stack
        v = VideoWriter([pressureBase individuals{i} '/' folders{f} '/field_incr10.avi'],'Uncompressed AVI');
        open(v);
        if rem((last-first),10)~=0,
            for j = 1:1:numFrames-1,
                im = imread([pressureBase individuals{i} '/' folders{f} '/field_' num2str(currFr,numformat) '.tif']);
                writeVideo(v,im);
                currFr = currFr + 10;
            end
        else
            for j = 1:1:numFrames-2,
                im = imread([pressureBase individuals{i} '/' folders{f} '/field_' num2str(currFr,numformat) '.tif']);
                writeVideo(v,im);
                currFr = currFr + 10;
            end
        end
        close(v);
    
        currFr = first;
        %Create a video file from the force figure image stack
        v = VideoWriter([pressureBase individuals{i} '/' folders{f} '/forces_incr10.avi'],'Uncompressed AVI');
        open(v);
        if rem((last-first),10)~=0,
            for j = 1:1:numFrames-1,
                im = imread([pressureBase individuals{i} '/' folders{f} '/forces_' num2str(currFr,numformat) '.tif']);
                writeVideo(v,im);
                currFr = currFr + 10;
            end
        else
            for j = 1:1:numFrames-2,
                im = imread([pressureBase individuals{i} '/' folders{f} '/forces_' num2str(currFr,numformat) '.tif']);
                writeVideo(v,im);
                currFr = currFr + 10;
            end
        end

        close(v);
    

        %Write kinematics proof video
        %Set an index equal to current frame number
        currFr = first;
        %Create a video file from the pressure field image stack
        v = VideoWriter([pressureBase individuals{i} '/' folders{f} '/kinematics_incr10.avi'],'Uncompressed AVI');
        open(v);
        if rem((last-first),10)~=0,
            for j = 1:1:numFrames-1,
                im = imread([pressureBase individuals{i} '/' folders{f} '/kinematics_' num2str(currFr,numformat) '.tif']);
                writeVideo(v,im);
                currFr = currFr + 10;
            end
        else
            for j = 1:1:numFrames-2,
                im = imread([pressureBase individuals{i} '/' folders{f} '/kinematics_' num2str(currFr,numformat) '.tif']);
                writeVideo(v,im);
                currFr = currFr + 10;
            end
        end
        close(v);

        %Save out force data (total, positive, negative x & y)
        xlswrite([pressureBase individuals{i} '/' folders{f} '/netForceX_incr10.xlsx'],totx)
        xlswrite([pressureBase individuals{i} '/' folders{f} '/netForceY_incr10.xlsx'],toty)
        xlswrite([pressureBase individuals{i} '/' folders{f} '/posForceX_incr10.xlsx'],totposx)
        xlswrite([pressureBase individuals{i} '/' folders{f} '/posForceY_incr10.xlsx'],totposy)
        xlswrite([pressureBase individuals{i} '/' folders{f} '/negForceX_incr10.xlsx'],totnegx)
        xlswrite([pressureBase individuals{i} '/' folders{f} '/negForceY_incr10.xlsx'],totnegy)


        %Save out midlines
        xlswrite([pressureBase individuals{i} '/' folders{f} '/midlines.xlsx'],midlineSaveX,'xvals');
        xlswrite([pressureBase individuals{i} '/' folders{f} '/midlines.xlsx'],midlineSaveY,'yvals');


        %Save out pressure, calculation boundaries, forces at each point in
        %boundary, and normal unit vectors for each force vector
        dlmwrite([pressureBase individuals{i} '/' folders{f} '/pressures.dat'], pressures, '\t');
        dlmwrite([pressureBase individuals{i} '/' folders{f} '/xvalsBound.dat'], boundariesx, '\t');
        dlmwrite([pressureBase individuals{i} '/' folders{f} '/yvalsBound.dat'], boundariesy, '\t');
        dlmwrite([pressureBase individuals{i} '/' folders{f} '/xForcesAll.dat'], forcexall, '\t');
        dlmwrite([pressureBase individuals{i} '/' folders{f} '/yForcesAll.dat'], forceyall, '\t');
        dlmwrite([pressureBase individuals{i} '/' folders{f} '/xUnitNormals.dat'], normsx, '\t');
        dlmwrite([pressureBase individuals{i} '/' folders{f} '/yUnitNormals.dat'], normsy, '\t');
        dlmwrite([pressureBase individuals{i} '/' folders{f} '/distOnMidline.dat'], dists, '\t');

        fileID = fopen([pressureBase individuals{i} '/' folders{f} '/locationofpoints.dat'], 'w');
        [nrows,ncols]=size(locs);
        for row = 1:1:nrows,
            for col = 1:1:ncols,
                if class(locs{row,col}) == 'cell',
                    locs{row,col} = locs{row,col}{1,1};
                end
            end
        end
        for row=1:1:nrows,
            fprintf(fileID,'%s\n',strjoin(locs(row,:),'\t'));
        end
        fclose(fileID);
        
        
%         return
    end

end