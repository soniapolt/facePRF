
%load('/Volumes/projects/invPRF/fixPRF/prfSets/fixPRF_tempCSSn_photo_bilat_r2-20.mat');
clear all; close all;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),'fixPRF',20,'photo','tempCSSn',{'lh' 'rh'},[]));
r = cellNum('mFus_faces',info.ROIs); c = 2; %upright



tw = [];sn = []; maxW = [];
for s = 1:6
    tempW = vertcat(subj(s).roi(r).fits(c).vox.tempW);
    sn = [sn s*ones(1,length(subj(s).roi(7).fits(2).vox))];
    
    [a,b] = max(tempW');
    maxW = [maxW b];
    tw = [tw; tempW];
    
end

[Y,loss] = tsne(tw);
figure;
gscatter(Y(:,1),Y(:,2),sn'); title('tSNE by subject number');
figure; 
gscatter(Y(:,1),Y(:,2),maxW'); title('tSNE by max weight');

