invX = vertcat(bFits(1).vox.Xdeg); uprX = vertcat(bFits(2).vox.Xdeg);
invY = vertcat(bFits(1).vox.Ydeg); uprY = vertcat(bFits(2).vox.Ydeg);
invSize = vertcat(bFits(1).vox.size); uprSize = vertcat(bFits(2).vox.size);

clear stepX stepY stepSize;
steps = 200;

for n = 1:length(invX)
    stepX(n,:) = linspace(uprX(n),invX(n),steps);
    stepY(n,:) = linspace(uprY(n),invY(n),steps);
    stepSize(n,:) = linspace(uprSize(n),invSize(n),steps);
end

figure;

%%% linear color steps between two colors (upright & inverted)
%col1 = condColors(2); col2 = condColors(1);
%stepCol = [linspace(col1(1),col2(1),steps);linspace(col1(2),col2(2),steps);linspace(col1(3),col2(3),steps)]';

%%% color space based on angle of shift
XYshift = vertcat(bFits(1).vox.XYdeg)-vertcat(bFits(2).vox.XYdeg);
    angles = rad2deg(atan2(XYshift(:,2),XYshift(:,1)));
    cmap = cmapang; % to account for our inverted axis
        vCol = cmap(ceil(angles./360*length(cmap)),:);
        

figure; clear M;
for n = 1:steps
    %subplot(4,5,n);
    %clf; 
    scratch_axis;    hold on;
    scatter(stepX(:,n),stepY(:,n),1,vCol.*(.75+.25*(n/steps)),'Filled'); hold on;
    %scatter(stepX(:,n),stepY(:,n),15,stepCol(n,:),'Filled');
    drawnow;
    M(n) = getframe;
end
hold on;scratch_axis;
scatter(stepX(:,n),stepY(:,n),15,vCol,'Filled');M(n) = getframe;

  % create the video writer with 1 fps
  writerObj = VideoWriter('/Users/sonia/Desktop/TalkPrep/myVideo','MPEG-4');
  writerObj.FrameRate = 20
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
writeVideo(writerObj,M);
% close the writer object
close(writerObj);


%%% upright coverage
figure; scratch_axis;
col1 = [0 135 203]/255;
    for v = 1:length(stepX)
    plotCircle(stepX(v,1),stepY(v,1),stepSize(v,1)/2,.33*col1,bFits(1).vox(v).r2/100*.5,'edge'); % center
    end
scatter(stepX(:,1),stepY(:,1),15,col1,'Filled'); hold on;
drawnow;

%%% inverted coverage
col2 = [71 47 146]/255;
figure;
scratch_axis;
scatter(stepX(:,n),stepY(:,n),15,col2,'Filled'); %hold on;
    for v = 1:length(stepX)
    plotCircle(stepX(v,n),stepY(v,n),stepSize(v,1)/2,.33*col2,bFits(1).vox(v).r2/100*.5,'edge'); % center
    end