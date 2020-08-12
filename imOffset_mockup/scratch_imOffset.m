clear all; close all;

which = {'outline' 'full' 'internal2' 'eyes'};


%%% for some basic params
load internalAvg.mat
faceSize = size(avgFace,1); ppd = faceSize/3.2;
imCenter = faceSize/2;

for n = 1:length(which)
    eval(['load ' which{n} 'Avg.mat;']);
    if max(avgFace(:)) > 1 avgFace = avgFace./255; end
    subplot(2,4,n);

imshow(avgFace); hold on;
[centX,centY] = weightedCent(avgFace);
h = hline(centY,'b-');set(h,'LineWidth',2)

offset(n) = (centY-imCenter)/ppd;
title(sprintf('%s: %.3fdeg from horiz.',which{n},offset(n)));

subplot(2,4,n+4); imshow(flipud(avgFace));hold on;
[centX,centY] = weightedCent(flipud(avgFace)); 
h = hline(centY,'b-');set(h,'LineWidth',2)
end
