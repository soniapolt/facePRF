% creates and saves a bootstrap estimates for relevant params on a pRF-sets
% struct. saves it back in the PRF-sets .mat file. i've spun this out of
% any plotting functions because it will be relatively time intensive
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

expt = 'fixPRF';
% assumes all subjects and all ROIs...
% to add later - bootstrapping within subjects

minR2 = ['r2-20'];          % cutoff for vox selection

whichStim = 'outline';
whichModel = 'kayCSS';
hems = {'lh' 'rh'};
bootPars = {'gain' 'r2' 'Y' 'X' 'size' 'eccen'};
plotIt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prfSet = pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);
info.boot = 'condition-wise bootstrapping';
info.numBoot = 1000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bootstrap o'clock                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(roi)
for b = 1:length(bootPars)
    for c = 1:length(roi(1).fits)
       pars = getPar(bootPars{b},roi(r).fits(c),1);
       
       roi(r).fits(c).boot(b).parName = bootPars{b};
       [roi(r).fits(c).boot(b).median ...
        roi(r).fits(c).boot(b).CI ...
        roi(r).fits(c).boot(b).dist] = bootstrapCI(pars,[],info.numBoot,[]);
    end
    
    if plotIt
        if b == 1 && c == 1 niceFig([.1 .1 .8 .6]); pl = 1; end
        subplot(2, length(bootPars),pl);
        niceHist(roi(r).fits(c).boot(b).dist,condColors(c),1);
        title(sprintf('Bootstrap Results: Median = %.2f, CI = [%.2f %.2f]',roi(r).fits(c).boot(b).median,...
        roi(r).fits(c).boot(b).CI(1),roi(r).fits(c).boot(b).CI(2)))
    end
end
end

save(pRFfile,'info','subj','roi','boot');

    %Save the coverage info 
    % saveFile = fullfile(savePath, [hems{h} '_coverage_data_ve', num2str(ve_cutoff*100)]);
   