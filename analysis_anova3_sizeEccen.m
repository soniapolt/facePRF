% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 3/11/20: THIS IS THE SIZE ANOVA CODE
% 5/6/20: changed size to == 1 sigma/sqrt(n)

clear all; close all;


expt = 'fixPRF';

minR2 = 20;          % cutoff for vox selection
whichANOVA = 'EVC';
ROIs= [standardROIs(whichANOVA)];%('face');% 'V3' 'hV4' 'IOG_faces' 'pFus_faces','mFus_faces'

whichStim = 'outline';
whichModel = 'kayCSS';%'cssExpN';%
fitSuffix = '';

hems = {'lh' 'rh'};
tests = {'slope' 'intercept'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
checkDir([dirOf(pwd) '/stats/' expt '/anova3/']);
fid = fopen([dirOf(pwd) '/stats/' expt '/anova3/ANOVA3_sizeEccen_' whichANOVA '-' num2str(minR2) '_' whichModel '_' whichStim '.txt'],'w+');
fprintf(fid,'\n**************\n%s\n**************\n',strTogether(ROIs));

fontSize = 11; titleSize = 14;
for t = 1:length(tests)
    anovaData= []; rmSubjs = [];
    
    
    for h = 1:length(hems)
        pf = pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems{h},fitSuffix);
        load(pf);fprintf(fid,'pRF file: %s\n',pf);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure 1: vox scatter + fit lines
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for r = 1:length(ROIs)
            niceFig([.1 .1 .8 .8],fontSize,1);
            
            ROInum = cellNum(ROIs{r},info.ROIs);
            sFit = struct;
            for s = 1:length(subj)
                fits = subj(s).roi(ROInum).fits;
                
                for c = 1:length(fits)
                    if c == 1 mult = .25; else mult = 1; end
                    
                    subplot(2,length(info.subjs)/2,s);
                    
                    x = [fits(c).vox.eccen]';
                    y = [fits(c).vox.size]'/2;
                    X = [x ones(length(fits(c).vox),1)];
                    
                    if ~isempty(X)
                        [sFit(c).h(s,:),R2(s)] = fitl1line(X,y);
                        h1 = plot(x,X*sFit(c).h(s,:)','Color',condColors(s,1)*mult); hold on;
                        hold on; scatter(x,y,5,condColors(s,1)*mult,'filled'); hold on; title([hems{h} ' ' ROIs{r} ' subj ' info.subjs{s}]);
                        set(h1,'LineWidth',1); %set(extLine,'LineStyle',':');
                        alpha(.2);
                        
                        anovaData = [anovaData; sFit(c).h(s,t) h r c s];
                    else rmSubjs(end+1) = s;
                    end
                end
                xl = xlim; xlim([0.25 6]);  ylim([0.25 6]);
                set(gca,'TickDir','out');
                xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (1*SD/sqrt(N)) (dva)'],'FontSize',fontSize);
                axis square;
                
            end
            
        end
    end
    
    fprintf('\n**************\n%s\n**************\n', strTogether(ROIs));
    
    % check for missing data and remove those subjects from the comparison
    [anovaData] = anova_rmSubjs(anovaData,rmSubjs);
    % run anova
    factNames = {'hem' 'ROI' 'condition'};
    result = rmAnova3(anovaData,factNames,0);
    %print output
    anova3_text(fid,result,tests{t});
end
%save('results/sizeEccen_anovaData.mat','anovaData');

%%% do anovaData style


playSound;