%Get fish midline & pect fins for current frame

% %parameters for testing/debugging
% frame = 43;
% % outlineXtra = [interface(:,2) interface(:,4)];
% outlineXtra = outline_forMid;
% mm2pixScale = 0.226674;
% numformat = '%05';
% delimiter = '/t';
% headerlines = 1;
% DaVis_x0 = 128;
% DaVis_y0 = 126;
% mov = VideoReader('H:/fish vids for pressure fields/trout/clips for steady analysis_long/trout7/trout7_2.5BLs_1_ventral_frm43-240.avi');          % read in movie
% vidfr = 1;

function [midline] = fish_midline_fctn(frame, outlineXtra, mm2pixScale, numformat, delimiter, headerlines, DaVis_x0, DaVis_y0)

fr = frame;
outline = outlineXtra;
%pix scale factor (mm/pix)
scale = mm2pixScale;
xmin = -DaVis_x0*scale;
xmax = xmin+(1024*scale);
ymin = (-1024+DaVis_y0)*scale;
ymax = ymin+(1024*scale);




xdata = outline(:,1);
ydata = outline(:,2);

% plot(xdata,ydata)
% return    



xoutline = xdata;
xoutline = [xoutline; xoutline(1)];
youtline = ydata;
youtline = [youtline; youtline(1)];



%     plot(xdata,ydata)
%     hold on
%     plot(xoutline,youtline,'b.')
%     
%     return






%Get midline


%Find nose point when it was digitized
nose_x_real = xdata(1);
nose_y_real = ydata(1);
xoutline = [xdata; xdata(1)];
youtline = [ydata; ydata(1)];

%Repeat above for tail
[maxxnum, maxidx] = max(xdata);
maxynum = ydata(maxidx);

if maxidx == 1 || maxidx == length(xdata),
    xtailseq = [xdata(end-2:end); xdata(1:3)];
    ytailseq = [ydata(end-2:end); ydata(1:3)];
else
    xtailseq = [xdata(maxidx:end); xdata(1:maxidx-1)];
    xtailseq = [xtailseq(end-2:end); xtailseq(1:4)];
    ytailseq = [ydata(maxidx:end); ydata(1:maxidx-1)];
    ytailseq = [ytailseq(end-2:end); ytailseq(1:4)];
end


%get curvature
curvtail = zeros(length(xtailseq),1);
for j = 2:1:length(xtailseq)-1,
    m1 = (ytailseq(j)-ytailseq(j-1))/(xtailseq(j)-xtailseq(j-1));
    m2 = (ytailseq(j+1)-ytailseq(j))/(xtailseq(j+1)-xtailseq(j));
    xc = (m1*m2*(ytailseq(j-1)-ytailseq(j+1))+m2*(xtailseq(j-1)+xtailseq(j))-m1*(xtailseq(j)+xtailseq(j+1)))/(2*(m2-m1));

    if abs(m1) == Inf,
        yc = mean([ytailseq(j) ytailseq(j-1)]);
        xc = (xtailseq(j)+xtailseq(j+1))/2-m2*(yc-((ytailseq(j)+ytailseq(j+1))/2));
        rad = sqrt((xtailseq(j)-xc)^2+(ytailseq(j)-yc)^2);
    elseif abs(m2) == Inf,
        yc = mean([ytailseq(j) ytailseq(j+1)]);
        xc = (xtailseq(j-1)+xtailseq(j))/2-m1*(yc-((ytailseq(j-1)+ytailseq(j))/2));
        rad = sqrt((xtailseq(j)-xc)^2+(ytailseq(j)-yc)^2);
    elseif m1 ~= 0,
        yc = -(1/m1)*(xc-(xtailseq(j-1)+xtailseq(j))/2)+(ytailseq(j-1)+ytailseq(j))/2;
        rad = sqrt((xtailseq(j)-xc)^2+(ytailseq(j)-yc)^2);
    elseif m1 == 0 && m2 ==0,
        rad = 100000000000;
    else
        yc = -(1/m2)*(xc-(xtailseq(j)+xtailseq(j+1))/2)+(ytailseq(j)+ytailseq(j+1))/2;
        rad = sqrt((xtailseq(j)-xc)^2+(ytailseq(j)-yc)^2);
    end
    curvtail(j) = 1/rad;        

end



%Assume true tail tip is halfway between max and next highest curv
%(adjacent)
[~, curvmaxidx] = max(curvtail);
maxx = xtailseq(curvmaxidx);
maxy = ytailseq(curvmaxidx);
[~, idx] = max([curvtail(curvmaxidx-1),curvtail(curvmaxidx+1)]);
if idx == 1,
    xnum = xtailseq(curvmaxidx-1);
    ynum = ytailseq(curvmaxidx-1);
else
    xnum = xtailseq(curvmaxidx+1);
    ynum = ytailseq(curvmaxidx+1);
end

tail_x_real = mean([maxx,xnum]);
tail_y_real = mean([maxy,ynum]);

%add new tail tip point to list of x and y data
if idx == 1,
    addidx = find(xdata==xnum);
    xdata = [xdata(1:addidx); tail_x_real; xdata(addidx+1:end)];
    ydata = [ydata(1:addidx); tail_y_real; ydata(addidx+1:end)];
else
    addidx = find(xdata==maxx);
    xdata = [xdata(1:addidx); tail_x_real; xdata(addidx+1:end)];
    ydata = [ydata(1:addidx); tail_y_real; ydata(addidx+1:end)];
end



%define nose and tail in pix
nose_x = nose_x_real/scale;%+x0;
nose_y = nose_y_real/scale;%+1024-y0;
tail_x = tail_x_real/scale;%+x0;
tail_y = tail_y_real/scale;%+1024-y0;



%Make binary image of fish
xpix = xdata/scale;%+x0;
ypix = ydata/scale;%+1024-y0;
im = poly2mask(xpix,ypix,1024,1024);
% imshow(im)
% axis on
% hold on
% set(gca,'YDir','normal')            %make sure y-axis is oriented correctly (positive up)
% plot(nose_x,nose_y,'r*')
% plot(tail_x,tail_y,'g*')
% % return

%get distances away from boundary
Daway = -bwdist(im);

im = imcomplement(im);




%get distance from boundary
D = bwdist(im);
D = D + Daway;
[rows cols] = size(D);
[xdim ydim] = meshgrid(1:1:cols,1:1:rows);




iterations = 10;

%alpha = elasticity weight, beta = flexibility weight, gamma = step size for moving x and
%y, and kappa is weight of external gradient
numpts = 20;
alpha = zeros(numpts);
for i = 1:1:numpts,
    if i < 5,
        alpha(i) = 0.5;
    elseif i < numpts*0.75,
        alpha(i) = 0.1;
    else
        alpha(i) = 0.05;
    end
end
beta = 0.3;
gamma = 0.7;
kappa = 0.5;



%initial line
xs = linspace(nose_x,tail_x,numpts);

init_index = round(xs);
ys = [];
for i = 1:1:length(init_index),
    [num idx] = max(D(:,init_index(i)));
    ys = [ys idx];
end
ys(1) = nose_y;
ys(end) = tail_y;

xs = xs';
ys = ys';

% plot(xs,ys,'r.')
% return

for i = 1:1:iterations,
%     i
    for j = 1:1:length(xs),
        if j == 1,
            continue
        elseif j==length(xs),
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

            %get forces
            F_distx = (kappa/2)*(Drt-Dlf);
            F_disty = (kappa/2)*(Ddn-Dup);



            xs(j) = xs(j) + gamma*(nansum([F_distx,F_stiff,F_elasticx]));
            ys(j) = ys(j) + gamma*(nansum([F_disty,F_stiff,F_elasticy]));



        end

    end


%         figure(1)
%         plot(xpix,ypix,'c-','LineWidth',2)
%         hold on
%         plot(xs,ys,'g*')
%         %axis([0 1024 0 1024])
%         hold off
%         pause(0.01)

%     return


end


% return



%Convert xs and ys from pix to mm
xs = (xs/size(im,2))*(xmax-xmin)+xmin;      % convert boundary coordinates from pixels to physical units
ys = (ys/size(im,1))*(ymax-ymin)+ymin;      % convert boundary coordinates from pixels to physical units

midline = [xs,ys];


%     xoutline = xoutline+xmin;
%     youtline = youtline+ymin;
%     figure(2)
%     plot(xoutline,youtline)
%     hold on
%     plot(xs,ys,'r-')
%     %plot(nose_x_real,nose_y_real,'bo')
%     axis([xmin xmax ymin ymax])
%     axis square
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15 0.15 0.7 0.7]);%[0 0 1 1]);
%     hold off
%     
% 
%     figure(3)
%     frame = flipud(read(mov,vidfr));
%     imshow(frame,'XData',[xmin xmax],'YData',[ymin ymax])
%     set(gca,'YDir','normal')
%     axis on
%     hold on
%     plot(xoutline,youtline,'LineWidth',2)
%     plot(xs,ys,'r-','LineWidth',2)
%     axis square
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15 0.15 0.7 0.7]);%[0 0 1 1]);
%     hold off
% 
% 
%     pause(0.1);
% return

end
