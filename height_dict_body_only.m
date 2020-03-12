%Make height dictionary

%Assumes that the fish is facing left, and the body is nearly horizontal.

%Assumes that the manually-digitized outlines are saved in Excel files of a
%particular structure (see below).

%Assumes InterX is loaded.


function [bodylength, bodydepth] = height_dict_body_only(bodyoutlinepath)



%read body pic metadata
%scale is in pix/mm
scale = xlsread(bodyoutlinepath,'Image data','A2');
xDimPix = xlsread(bodyoutlinepath,'Image data','A6');
yDimPix = xlsread(bodyoutlinepath,'Image data','B6');
xDimReal = xlsread(bodyoutlinepath,'Image data','A7');
yDimReal = xlsread(bodyoutlinepath,'Image data','B7');


%Read outlines of the body
outline = xlsread(bodyoutlinepath,'Outline');

%Indicate where the nose point is in the body outline
xnose_real = outline(1,1);
ynose_real = outline(1,2);

%Find the tail fork in the body outline
fork = -1;
i=0;
%Search along upper part of fish traveling down the list of x values.  
%When hit the point where x swaps from increasing to decreasing,
%indicate we're on the upper lobe.  Here, start search for fork.
while fork<0,
    
    i=i+1;
    if outline(i,1) > outline(i+1,1),
        fork = 0;
        break
    end
    
end
%Search the caudal fin for the point where x stops decreasing and starts
%increasing again.  This point is the fork.
while fork<1,
    i=i+1;
    if outline(i,1) < outline(i+1,1),
        fork = i;
        break
    end
end

%Indicate "tail" position as the fork.
xtail_real = outline(fork,1);
ytail_real = outline(fork,2);

%We need to cut the lobes of the caudal fin off of the body outline in order to
%draw a shape that is easy for the midline extractor to find a midline for.
%Cut based on the shape of the tail fork.  Here, get the directions the
%upper and lower lobes diverge from the fork.
slope1 = (outline(fork+1,2)-outline(fork,2))/(outline(fork+1,1)-outline(fork,1));
slope2 = (outline(fork,2)-outline(fork-1,2))/(outline(fork,1)-outline(fork-1,1));

%Divide the outline into upper and lower halves of the body.
upperhalf = flipud(outline(1:fork,1:2));
lowerhalf = outline(fork+1:end,1:2);

%Make the cut by finding the point on the body outline that matches where we
%desire to cut the best.  This means that the slope between the point on
%the outline and the fork point is similar to the slope for cutting.
upperTestSlopes = (upperhalf(:,2) - ytail_real)./(upperhalf(:,1) - xtail_real) - slope1;
lowerTestSlopes = (ytail_real - lowerhalf(:,2))./(xtail_real - lowerhalf(:,1)) - slope2;
[~,ui] = min(abs(upperTestSlopes));
[~,li] = min(abs(lowerTestSlopes));
upperEnd = fork+1-ui;
lowerEnd = li+fork;

% for i = 2:1:length(upperhalf)-1,
%    
%     testSlope = (upperhalf(i,2)-ytail_real)/(upperhalf(i,1)-xtail_real);
%     
%     if abs(testSlope - slope1) < 0.5,
%         upperEnd = fork+1-i;
%         break
%     end
% end
% 
% 
% for i = 1:1:length(lowerhalf)-1,
%    
%     testSlope = (ytail_real-lowerhalf(i,2))/(xtail_real-lowerhalf(i,1));
%     
%     if abs(testSlope - slope2) < 0.5,
%         lowerEnd = i+fork;
%         break
%     end
% end

%Make "cut" body outline x and y points lists.
xdata = [outline(1:upperEnd,1); xtail_real; outline(lowerEnd:end,1)];
ydata = [outline(1:upperEnd,2); ytail_real; outline(lowerEnd:end,2)];


%from the snake_midline script (not calling as a function due to its
%different assumptions)



%define nose and "tail" (really fork) in pixels
nose_x = xnose_real/scale;
nose_y = (yDimReal-ynose_real)/scale;
tail_x = xtail_real/scale;
tail_y = (yDimReal-ytail_real)/scale;



%Make binary image of fish - fish is white, background is black
xpix = xdata/scale;
ypix = (yDimReal-ydata)/scale;
im = poly2mask(xpix,ypix,yDimPix,xDimPix);
% imshow(im)
% axis on
% set(gca,'YDir','normal')            %make sure y-axis is oriented correctly (positive up)
% return


%get distances away from black/white boundary.  bwdist measures only black
%pixels
Daway = -bwdist(im);

%Invert image colors so that we can get distances from boundary inside fish
im = imcomplement(im);
%get distance from boundary
D = bwdist(im);

%Put all distance information together
D = D + Daway;
[rows, cols] = size(D);
[xdim, ydim] = meshgrid(1:1:cols,1:1:rows);


%set number of iterations will use to find midline
iterations = 10;

%Set up snake parameters
%alpha = elasticity weight, beta = flexibility weight, gamma = step size for moving x and
%y, and kappa is weight of external gradient
numpts = 20;
alpha = zeros(numpts,1);
for i = 1:1:numpts,
    if i < 5,
        alpha(i) = 0.5;
    elseif i < numpts*0.75,
        alpha(i) = 0.1;
    else
        alpha(i) = 0.05;
    end
end
beta = 0.8;
gamma = 0.9;
kappa = 0.5;



%initial line that we'll move to get the midline
xs = linspace(nose_x,tail_x,numpts);

init_index = round(xs);
ys = [];
for i = 1:1:length(init_index),
    [~, idx] = max(D(:,init_index(i)));
    ys = [ys idx];
end
%make sure line starts at nose and ends at tail
ys(1) = nose_y;
ys(end) = tail_y;

xs = xs';
ys = ys';


%return

%for each iteration,
for i = 1:1:iterations,
    %for each point the midline-to-be,
    for j = 1:1:length(xs),
        if j == 1,
            %don't move the nose point
            continue
        elseif j==length(xs),
            %don't move the tail point
            continue
        else
            %set up elastic force on each point - force keeping points 
            %evenly spread out.  Basically, saying here that we want the
            %average distance between points to be about the same.
            F_elasticx = alpha(j)*2*((xs(j-1)-xs(j))+(xs(j+1)-xs(j)));
            F_elasticy = alpha(j)*2*((ys(j-1)-ys(j))+(ys(j+1)-ys(j)));

            %set up stiffness force on each point - force keeping points 
            %from making too bendy a line.  Only works >2 pts from end.
            %will do this by preventing curvature (1/radcurv) from being
            %too big.  Like above, we'll do this by saying we want the
            %average curvature between adjacent points to not be very 
            %different.  This means we need the curvature at j-1, j, and
            %j+1.

            if j < 3,
                F_stiff = 0.1;
            elseif j> length(xs)-3,
                F_stiff = 0.1;
            else
                %Curvature at j-1
                m1 = (ys(j-1)-ys(j-2))/(xs(j-1)-xs(j-2));
                m2 = (ys(j)-ys(j-1))/(xs(j)-xs(j-1));
                xc = (m1*m2*(ys(j-2)-ys(j))+m2*(xs(j-2)+xs(j-1))-m1*(xs(j-1)+xs(j)))/(2*(m2-m1));
                yc = -(1/m1)*(xc-(xs(j-2)+xs(j-1))/2)+(ys(j-2)+ys(j-1))/2;
                rad1 = sqrt((xs(j-1)-xc)^2+(ys(j-1)-yc)^2);

                %Curvature at j
                m1 = (ys(j)-ys(j-1))/(xs(j)-xs(j-1));
                m2 = (ys(j+1)-ys(j))/(xs(j+1)-xs(j));
                xc = (m1*m2*(ys(j-1)-ys(j+1))+m2*(xs(j-1)+xs(j))-m1*(xs(j)+xs(j+1)))/(2*(m2-m1));
                yc = -(1/m1)*(xc-(xs(j-1)+xs(j))/2)+(ys(j-1)+ys(j))/2;
                rad2 = sqrt((xs(j)-xc)^2+(ys(j)-yc)^2);

                %Curvature at j+1
                m1 = (ys(j+1)-ys(j))/(xs(j+1)-xs(j));
                m2 = (ys(j+2)-ys(j+1))/(xs(j+2)-xs(j+1));
                xc = (m1*m2*(ys(j)-ys(j+2))+m2*(xs(j)+xs(j+1))-m1*(xs(j+1)+xs(j+2)))/(2*(m2-m1));
                yc = -(1/m1)*(xc-(xs(j)+xs(j+1))/2)+(ys(j)+ys(j+1))/2;
                rad3 = sqrt((xs(j+1)-xc)^2+(ys(j+1)-yc)^2);

                %Get stiffness force
                F_stiff = beta*2*((1/rad1-1/rad2)+(1/rad3-1/rad2));
            end




            %Last, set up force to maximize distance from the boundary.
            %If distance at x+1 is bigger, move that way, etc.

            %get distance value at x+1 and x-1
            Dlf = interp2(xdim,ydim,D,xs(j)-2,ys(j));
            Drt = interp2(xdim,ydim,D,xs(j)+2,ys(j));
            %get distance value at y+1 and y-1
            Dup = interp2(xdim,ydim,D,xs(j),ys(j)-2);
            Ddn = interp2(xdim,ydim,D,xs(j),ys(j)+2);

            %get distance forces
            F_distx = (kappa/2)*(Drt-Dlf);
            F_disty = (kappa/2)*(Ddn-Dup);


            %apply force to each point and take a proportional step
            xs(j) = xs(j) + gamma*(nansum([F_distx,F_stiff,F_elasticx]));
            ys(j) = ys(j) + gamma*(nansum([F_disty,F_stiff,F_elasticy]));



        end

    end


%         figure(1)
%         plot(xpix,ypix,'c-','LineWidth',2)
%         hold on
%         plot(xs,ys,'g*')
%         %axis([0 xDimPix 0 yDimPix])
%         axis equal
%         hold off
%         pause(0.01)

    %return
    %M(i) = getframe(gcf);                                             % store frame for movie


end

%movie2avi(M, 'midlinesnake.avi')

%return



%Convert xs and ys (midline) from pixels to mm
xs = (xs*scale);
ys = yDimReal - (ys*scale);


% plot(outline(:,1),outline(:,2))
% hold on
% % plot(xs,ys,'r.')
% axis([0 xDimReal 0 yDimReal])

%To get body depths, we need more internal points than the 20 on the
%midline.  Interpolate to find more.
xsLong = (linspace(xnose_real,xtail_real,50))';
ysLong = smooth(spline(xs,ys,xsLong));

% plot(xsLong,ysLong,'k.')
% hold off


%Need to make an accurate cut to get rid of the tail lobes - because there
%are two of them, we need to deal with them independently of the rest of
%the body/caudal fin.  Make the cut at the fork, on a line perpendicular to
%the midline.
%Get perpendicular slope
tailLobesCutSlope = -1/((ysLong(end) - ysLong(end-1))/(xsLong(end) - xsLong(end-1)));

%Get the slope between the fork point and upper/lower outline points.
%We'll need to travel from nose to tail along the outline, and find the
%slope that matches where we want to cut the best.
upperSlope = (upperhalf(:,2)-ytail_real)./(upperhalf(:,1)-xtail_real);
lowerSlope = (ytail_real-lowerhalf(:,2))./(xtail_real-lowerhalf(:,1));

[~, upperClosest] = min(abs(upperSlope-tailLobesCutSlope));
[~, lowerClosest] = min(abs(lowerSlope-tailLobesCutSlope));

if outline(fork+1-upperClosest,1) > outline(fork,1),
    upperClosest = upperClosest + 1;
    trigUp = 1;
end
if outline(lowerClosest+fork,1) > outline(fork,1),
    lowerClosest - lowerClosest - 1;
    trigLo = 1;
end

%Get indices of the outline points where the upper lobe starts in the
%outline, and where the lower lobe ends
upperClosest = fork+1-upperClosest;
lowerClosest = lowerClosest+fork;

%Get the lobes' x and y values as a list
xTailLobe = outline(upperClosest:lowerClosest,1);
yTailLobe = outline(upperClosest:lowerClosest,2);

% plot(xTailLobe,yTailLobe,'g-')


%Find the body depth at the nose
heightNose = sqrt((outline(1,1)-outline(end,1))^2+(outline(1,2)-outline(end,2))^2);

%Get the upper and lower body outline (minus tail lobes)
if exist('trigUp','var'),
    upper = outline(1:upperClosest+1,:)';
else
    upper = outline(1:upperClosest,:)';
end
if exist('trigLo','var'),
    lower = outline(lowerClosest-1:end,:)';
else
    lower = outline(lowerClosest:end,:)';
end

%Initialize storage for body length and depth
bodylength = zeros(length(xsLong-1),1);
bodydepth = zeros(length(xsLong-1),1);
%Assign body depth at the nose and tail
bodydepth(1) = heightNose;

%Initialize storage for x and y values on the lower body outline (for plotting
%purposes only)
xforPlotBody = zeros(length(xsLong),1);
xforPlotBody(1) = outline(end,1);
yforPlotBody = zeros(length(xsLong),1);
yforPlotBody(1) = outline(end,2);

%Initialize storage for the angle (from horizontal) of the perpendicular
%lines we used to measure body depth (for plotting purposes)
angleBody = zeros(length(xsLong),1);
%Set angle for the nose and tail
if outline(1,1) > outline(end,1),
    theta = -(90-acosd(abs(outline(1,1)-outline(end,1))/heightNose));
else
    theta = 90 - acosd(abs(outline(1,1)-outline(end,1))/heightNose);
end
angleBody(1) = theta;



%for each internal point on the midline,
for i = 2:1:length(xsLong)-1,
    
    %Find the slope of the perpendicular line
    slope = polyfit(xsLong(i-1:i),ysLong(i-1:i),1);
    slope = -1/(slope(1));
    
    %Make a giant perpendicular line
    xnew = [xsLong(i)-100; xsLong(i)+100];
    ynew = [(slope*(xnew(1)-xsLong(i))+ysLong(i)); (slope*(xnew(2)-xsLong(i))+ysLong(i))];
    
    new = [xnew ynew]';
    
    %Find where giant perpendicular line intersects the upper and lower
    %body outlines
    P1 = InterX(upper,new);
    P2 = InterX(lower,new);

    %Mainly for the last internal point.  Sometimes, the perpendicular line
    %just misses the outline (tail subtracted).  If that's the case, search
    %again with one more point on the outline, and see if we can then make
    %a measurement.
    if isempty(P1),
        upperAdd = [upper outline(size(upper,2)+1,:)'];
        P1 = InterX(upperAdd,new);
    end
    if isempty(P2),
        addIdx = size(outline,1) - size(lower,2);
        lowerAdd = [outline(addIdx,:)' lower];
        P2 = InterX(lowerAdd,new);
    end
    
    %Save x and y values of the lower body intersection points, and the
    %angle from horizontal (for plotting purposes)
    xforPlotBody(i) = P2(1);
    yforPlotBody(i) = P2(2);
    
    %Calculate depth of body
    d = sqrt((P1(2)-P2(2))^2+(P1(1)-P2(1))^2);
    
    if P1(1) > P2(1),
        theta = -(90-acosd(abs(P1(1)-P2(1))/d));
    else
        theta = 90-acosd(abs(P1(1)-P2(1))/d);
    end
    angleBody(i) = theta;    
    
    bodydepth(i) = d;
    
    bodylength(i) = bodylength(i-1) + abs(sqrt((ysLong(i)-ysLong(i-1))^2+(xsLong(i)-xsLong(i-1))^2));
    

end



%Now will look at the caudal lobes.  We need to extend the midline at the
%tail outward so that we can keep finding perpendicular lines on which to
%measure body depth.

%Get the spacing between points on the midline
spacing = xsLong(2)-xsLong(1);

%Extend the midline 10 steps onward, but since the lobes are small and
%depth changes rapidly, use 30 points instead of 10
xextend = linspace(xtail_real,10*spacing+xtail_real,30);
%The midline will be extended as a straight line continuation of the last
%line segment in the midline.  Get the slope of the last segment.
slopeTail = polyfit(xsLong(49:50),ysLong(49:50),1);
slopeTail = slopeTail(1);

%Use slope, fork point, and x-values in line formula to get y values.
yextend = slopeTail*(xextend-xtail_real)+ytail_real;

%We listed all the points in the lobes in one list before.  Divide into two
%lists at the fork point (so upper and lower lobes are listed separately).
idx = find(xTailLobe == xtail_real);
tailTop = [xTailLobe(1:idx) yTailLobe(1:idx)]';
tailBot = [xTailLobe(idx:end) yTailLobe(idx:end)]';

%Get the perpendicular slope to the extended midline
slope = -1/slopeTail;
%Get the distance we travel on the extended midline every time we take a
%step
step = sqrt((xextend(2)-xextend(1))^2+(yextend(2)-yextend(1))^2);

% %Initialize storage for the x and y values from the ventral side of each
% %lobe (for plotting purposes)
% xforPlotTailUpper = [];
% yforPlotTailUpper = [];
% angleTailUpper = [];
% 
% xforPlotTailLower = [];
% yforPlotTailLower = [];
% angleTailLower = [];
% 
% %Store the depths specific to each lobe (for plotting purposes)
% htLower = [];
% htUpper = [];
% 
% plot(outline(:,1),outline(:,2))
% hold on
% axis([0 xDimReal 0 yDimReal])





%Travel down the extended midline.
for i = 1:1:length(xextend),
    
    %Make a giant perpendicular line
    xnew = [xextend(i)-100; xextend(i)+100];
    ynew = [(slope*(xnew(1)-xextend(i))+yextend(i)); (slope*(xnew(2)-xextend(i))+yextend(i))];
    
    new = [xnew ynew]';
    
    %Find points where the giant midline intersect the two lobes
    P1 = InterX(tailTop,new);
    P2 = InterX(tailBot,new);
    
    
    if size(P1,2) < 2 && size(P2,2) < 2,
        if i == 1,
            %If we're at the tail fork, measure depth the same way as on
            %body.  Add to body list.
            xforPlotBody(end) = P2(1);
            yforPlotBody(end) = P2(2);

            d = sqrt((P1(2)-P2(2))^2+(P1(1)-P2(1))^2);

%             if P1(1) > P2(1),
%                 theta = -(90-acosd(abs(P1(1)-P2(1))/d));
%             else
%                 theta = 90-acosd(abs(P1(1)-P2(1))/d);
%             end
%             angleBody(end) = theta; 
%             plot(xnew,ynew,'m-')
        else
            %If we have no intersection points on either lobe, we've gotten to the
            %end of the tail.  Stop measuring.
            break
        end
        
    elseif size(P1,2) < 2 && size(P2,2) > 1,
        %If we didn't find multiple intersection points on the upper lobe
        %(either tangent to or off the lobe), but we did intersect the
        %lower lobe, focus on lower lobe only.
        
        %Get distance between intersection points
        d = sqrt((P2(2,2)-P2(2,1))^2+(P2(1,2)-P2(1,1))^2);
        
%         %Figure out which point in the intersection points list is closest
%         %to the fork point.  Option 1: the first point is.  Option 2: the
%         %second point is.
%         op1 = sqrt((P2(1,1)-xextend(i))^2+(P2(2,1)-yextend(i))^2);
%         op2 = sqrt((P2(1,2)-xextend(i))^2+(P2(2,2)-yextend(i))^2);
%         
%         [~, idx] = min([op1 op2]);
%         
%         %Mark the point closest to the fork appropriately.
%         if idx == 1,
%             xclose = P2(1,1);
%             yclose = P2(2,1);
%             xfar = P2(1,2);
%             yfar = P2(2,2);
%         else
%             xclose = P2(1,2);
%             yclose = P2(2,2);
%             xfar = P2(1,1);
%             yfar = P2(2,1);
%         end
%         
%         %Find the angle of the perpendicular line, for plotting purposes
%         %only.
%         if xclose > xfar
%             theta = 360-(90-acosd(abs(xclose-xfar)/d));
%         else
%             theta = 90-acosd(abs(xclose-xfar)/d);
%         end
%         
%         %Save list of x and y values of the ventral intersection point, and
%         %the depth of the specific lobe
%         angleTailLower = [angleTailLower; theta];
%         xforPlotTailLower = [xforPlotTailLower; xfar];
%         yforPlotTailLower = [yforPlotTailLower; yfar];
%         htLower = [htLower; d];
%         
% %         %For debugging...plot the perpendicular line
% %         plot(xnew,ynew,'r-')


    elseif size(P2,2) < 2 && size(P1,2) > 1,
        %If we didn't find multiple intersection points on the lower lobe
        %(either tangent to or off the lobe), but we did intersect the
        %upper lobe, focus on upper lobe only.
        
        %Get distance between intersection points
        d = sqrt((P1(2,2)-P1(2,1))^2+(P1(1,2)-P1(1,1))^2);
        
%         %Figure out which point in the intersection points list is closest
%         %to the fork point.  Option 1: the first point is.  Option 2: the
%         %second point is.
%         op1 = sqrt((P1(1,1)-xextend(i))^2+(P1(2,1)-yextend(i))^2);
%         op2 = sqrt((P1(1,2)-xextend(i))^2+(P1(2,2)-yextend(i))^2);
%         
%         [~, idx] = min([op1 op2]);
%         
%         %Mark the point closest to the fork appropriately.
%         if idx == 1,
%             xclose = P1(1,1);
%             yclose = P1(2,1);
%             xfar = P1(1,2);
%             yfar = P1(2,2);
%         else
%             xclose = P1(1,2);
%             yclose = P1(2,2);
%             xfar = P1(1,1);
%             yfar = P1(2,1);
%         end
% 
%         %Find the angle of the perpendicular line, for plotting purposes
%         %only.
%         if xfar > xclose
%             theta = 360-(90-acosd(abs(xfar-xclose)/d));
%         else
%             theta = 90-acosd(abs(xfar-xclose)/d);
%         end
%         
%         %Save list of x and y values of the ventral intersection point, and
%         %the depth of the specific lobe
%         angleTailUpper = [angleTailUpper; theta];
%         xforPlotTailUpper = [xforPlotTailUpper; xclose];
%         yforPlotTailUpper = [yforPlotTailUpper; yclose];
%         htUpper = [htUpper; d];
%         
% %         %For debugging...plot the perpendicular line
% %         plot(xnew,ynew,'r-')

        
        
    else
        %If we intersect both lobes,
        
%         %Get depth on the two lobes (for plotting purposes)
%         dLower = sqrt((P2(2,2)-P2(2,1))^2+(P2(1,2)-P2(1,1))^2);
%         dUpper = sqrt((P1(2,2)-P1(2,1))^2+(P1(1,2)-P1(1,1))^2);
        
        %Get total depth
        d = sqrt((P2(2,2)-P2(2,1))^2+(P2(1,2)-P2(1,1))^2) + sqrt((P1(2,2)-P1(2,1))^2+(P1(1,2)-P1(1,1))^2);
        
%         %Figure out which points are closest to the fork (for plotting
%         %purposes)
%         op1P2 = sqrt((P2(1,1)-xextend(i))^2+(P2(2,1)-yextend(i))^2);
%         op2P2 = sqrt((P2(1,2)-xextend(i))^2+(P2(2,2)-yextend(i))^2);
%         
%         [~, idxP2] = min([op1P2 op2P2]);
%         
%         op1P1 = sqrt((P1(1,1)-xextend(i))^2+(P1(2,1)-yextend(i))^2);
%         op2P1 = sqrt((P1(1,2)-xextend(i))^2+(P1(2,2)-yextend(i))^2);
%         
%         [~, idxP1] = min([op1P1 op2P1]);
%         
%         %Angle will be the same for top and bottom lobe, so only assign
%         %"close" and "far" on the lower lobe.
%         if idxP2 == 1,
%             xclose = P2(1,1);
%             yclose = P2(2,1);
%             xfar = P2(1,2);
%             yfar = P2(2,2);
%         else
%             xclose = P2(1,2);
%             yclose = P2(2,2);
%             xfar = P2(1,1);
%             yfar = P2(2,1);
%         end
%         
%         %Find the angle of the line.
%         if xclose > xfar
%             theta = 360-(90-acosd(abs(xclose-xfar)/dLower));
%         else
%             theta = 90-acosd(abs(xclose-xfar)/dLower);
%         end
%         
%         %Save the x and y values of the ventral intersection on each lobe,
%         %and the angle of the perpendicular measurement line (for plotting
%         %purposes)
%         angleTailLower = [angleTailLower; theta];
%         xforPlotTailLower = [xforPlotTailLower; xfar];
%         yforPlotTailLower = [yforPlotTailLower; yfar];
%         
%         angleTailUpper = [angleTailUpper; theta];
%         xforPlotTailUpper = [xforPlotTailUpper; P1(1,idxP1)];
%         yforPlotTailUpper = [yforPlotTailUpper; P1(2,idxP1)];
%         
%         %Store depth on each lobe (for plotting purposes)
%         htUpper = [htUpper; dUpper];
%         htLower = [htLower; dLower];
%         
%         if i==1,
%             xforPlotBody(end) = [];
%             yforPlotBody(end) = [];
%         end
% %         %For debugging...plot the perpendicular line
% %         plot(xnew,ynew,'g-')



    end

    %Add the length along the body and the total lobe depth to the lists
    if i == 1,
        bodydepth(end) = d;
    else
        bodydepth = [bodydepth; d];
    end

    if i == 1,
        next = bodylength(end-1) + abs(sqrt((ytail_real-ysLong(end-1))^2+(xtail_real-xsLong(end-1))^2));
        bodylength(end) = next;
    else
        next = bodylength(end)+step;
        bodylength = [bodylength; next];
    end
    

end

%Normalize the body length to a decimal percentage
bodylength = bodylength/max(bodylength);


% %Plot the body and tail depths on top of the original outline as rectangles
% plot(outline(:,1), outline(:,2))
% hold on
% width = 0.25;
% %For each point along the fish's length,
% for i = 1:1:length(bodylength),
%     %If we're on the body (not lobes),
%     if i <= length(xforPlotBody),
%         %Find the x and y value of the ventral intersection point, and the
%         %body depth there
%         xpos = xforPlotBody(i);
%         ypos = yforPlotBody(i);
%         ht = bodydepth(i);
%         %Make a list of points for the rectangle (rectangle is vertical.
%         %Will rotate it to the angle of the perpendicular line)
%         xlist = [xpos; xpos - width/2; xpos - width/2; xpos + width/2; xpos + width/2; xpos];
%         ylist = [ypos; ypos; ypos + ht; ypos + ht; ypos; ypos];
%         %Because rotation occurs around the origin, translate the rectangle
%         %onto the origin.
%         xlist = xlist - xpos;
%         ylist = ylist - ypos;
%         %Rotate the rectangle
%         xrot = xlist*cos(deg2rad(angleBody(i))) - ylist*sin(deg2rad(angleBody(i)));
%         yrot = xlist*sin(deg2rad(angleBody(i))) + ylist*cos(deg2rad(angleBody(i)));
%         %Translate the rectangle back to its starting point
%         xrot = xrot + xpos;
%         yrot = yrot + ypos;
%         %Plot the rectangle
%         plot(xrot,yrot,'k-')
%     else
%         %Look at the lobes.
%         
%         idx = i - length(xforPlotBody);
%         
%         %Figure out how many measurements we have on each lobe (may not be
%         %equal!)
%         upperLen = length(xforPlotTailUpper);
%         lowerLen = length(xforPlotTailLower);
%         
%         %Start with upper lobe.  For each measurement, do the same thing we
%         %did to body points.
%         if idx <= upperLen,
%             xpos = xforPlotTailUpper(idx);
%             ypos = yforPlotTailUpper(idx);
%             ht = htUpper(idx);
%             xlist = [xpos; xpos - width/2; xpos - width/2; xpos + width/2; xpos + width/2; xpos];
%             ylist = [ypos; ypos; ypos + ht; ypos + ht; ypos; ypos];
%             xlist = xlist - xpos;
%             ylist = ylist - ypos;
%             xrot = xlist*cos(deg2rad(angleTailUpper(idx))) - ylist*sin(deg2rad(angleTailUpper(idx)));
%             yrot = xlist*sin(deg2rad(angleTailUpper(idx))) + ylist*cos(deg2rad(angleTailUpper(idx)));
%             xrot = xrot + xpos;
%             yrot = yrot + ypos;
%             plot(xrot,yrot,'k-')
%         end
%         %Now do lower lobe.  For each measurement, do the same thing we
%         %did to body points.
%         if idx <= lowerLen,
%             xpos = xforPlotTailLower(idx);
%             ypos = yforPlotTailLower(idx);
%             ht = htLower(idx);
%             xlist = [xpos; xpos - width/2; xpos - width/2; xpos + width/2; xpos + width/2; xpos];
%             ylist = [ypos; ypos; ypos + ht; ypos + ht; ypos; ypos];
%             xlist = xlist - xpos;
%             ylist = ylist - ypos;
%             xrot = xlist*cos(deg2rad(angleTailLower(idx))) - ylist*sin(deg2rad(angleTailLower(idx)));
%             yrot = xlist*sin(deg2rad(angleTailLower(idx))) + ylist*cos(deg2rad(angleTailLower(idx)));
%             xrot = xrot + xpos;
%             yrot = yrot + ypos;
%             plot(xrot,yrot,'k-')
%         end
%         
%     end
%     
% 
% end








end





