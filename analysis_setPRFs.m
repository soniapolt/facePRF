% for the plotCode update, we generate and save a struct of pRF fits for
% experiments that can be efficiently used and manipulated by other

clear all; close all;

info.subjs = {'MG' 'JG' 'TH' 'EM' 'DF' 'SP'};%
info.task = '';
info.expt = 'fixPRF';
info.setNotes = '';

info.minR2 = 20;          % cutoff for vox selection
info.ROIs= standardROIs;

info.whichStim = 'photo';
info.whichModel ='cssExpN';%  'kayCSS';%%'cssShift';%
info.hems = {'rh' 'lh'};
info.fitSuffix = ''; %'_orig';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pRF loading - collapse across subjects    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear fits; clear roi;

for r = 1:length(info.ROIs)
    init = 1;
for s = 1:length(info.subjs)
    sinit = 1;
        [session, numRuns] = vpnlSessions(info.expt,info.subjs{s}); % OPTIONAL: SESSNUM, TASK

        for h = 1:length(info.hems)
            
            thisROI = [info.hems{h} '_' info.ROIs{r}];
            [~, fitsName] = fitsDirs(dirOf(pwd),info.expt,session,info.whichStim,vpnlROI(thisROI,session(1:2),info.expt),info.whichModel,info.fitSuffix);
            if exist(fitsName)>0
                load(fitsName);
            if sinit   
                subj(s).roi(r).fits = fits; sinit = 0;% keep subjects in distinct fields
            else
                for c = 1:length(fits)
                    subj(s).roi(r).fits(c).vox = [subj(s).roi(r).fits(c).vox fits(c).vox];
                end
            end
            % add subject/stim info to each voxel representation
            for f = 1:length(fits)
                for v = 1:length(fits(f).vox)
                    fits(f).vox(v).stim = fits(f).stim; end
            end
            
            if init 
                roi(r).fits = fits; init = 0; % across-subject struct
            else
                for c = 1:length(fits)
                    roi(r).fits(c).vox = [roi(r).fits(c).vox fits(c).vox];
                end
            end
             
            else fprintf('Missing fits in: %s...\n',fitsName);
            end
        end % hems
    end %rois
end %subjs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate fit params into more meaningful terms for invPRF expt              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(info.ROIs)
    % across subjects
        for c = 1:length(roi(r).fits)
            if isfield(fits,'expN')
            roi(r).fits(c).vox = readPRFs(roi(r).fits(c).vox,fits(1).ppd,fits(1).res, fits(1).expN);else
            roi(r).fits(c).vox = readPRFs(roi(r).fits(c).vox,fits(1).ppd,fits(1).res,[]);end
    % each subject    
        for s = 1:length(info.subjs)
            if isfield(fits,'expN')
            subj(s).roi(r).fits(c).vox = readPRFs(subj(s).roi(r).fits(c).vox,fits(1).ppd,fits(1).res, fits(1).expN);else
            subj(s).roi(r).fits(c).vox = readPRFs(subj(s).roi(r).fits(c).vox,fits(1).ppd,fits(1).res,[]);end
        end
        end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % which voxels are we plotting? trim edge values, R2 cutoff                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % across subjects
    plotVox = 1:length(roi(r).fits(1).vox);
    for v = 1:length(roi(r).fits(1).vox)
        for c = 1:length(roi(r).fits)
            if  roi(r).fits(c).vox(v).r2 < info.minR2 ...
                    || trimPRFs(roi(r).fits(c).vox(v).params,roi(r).fits(c).vox(v).betas,fits(1).ppd,fits(1).res);
                plotVox(find(plotVox==v)) = []; end
        end
    end
 
    for c = 1:length(roi(r).fits)
        roi(r).fits(c).vox = roi(r).fits(c).vox(plotVox);
    end
    
    % each subject
    for s = 1:length(info.subjs)
    plotVox = 1:length(subj(s).roi(r).fits(1).vox);
    for v = 1:length(subj(s).roi(r).fits(1).vox)
        for c = 1:length(subj(s).roi(r).fits)
            if  subj(s).roi(r).fits(c).vox(v).r2 < info.minR2 ...
                    || trimPRFs(subj(s).roi(r).fits(c).vox(v).params,subj(s).roi(r).fits(c).vox(v).betas,fits(1).ppd,fits(1).res);
                plotVox(find(plotVox==v)) = []; end
        end
    end
 
    for c = 1:length(subj(s).roi(r).fits)
        subj(s).roi(r).fits(c).vox = subj(s).roi(r).fits(c).vox(plotVox);
    end    
    end
end

saveFile = pRFfile(dirOf(pwd),info.expt,info.minR2,info.whichStim,info.whichModel,info.hems,info.fitSuffix,info.task);
fprintf('Saving %s...\n',saveFile);
save(saveFile,'roi','subj','info');
