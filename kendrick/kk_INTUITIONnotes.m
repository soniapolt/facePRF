% constants
res = 128;          % resolution of simulations
sd = res/(3*4);     % sd for the gaussian (+/- 2 sd will take up 1/3 of side)
sizes = [res/16 res/8 res/4 res/2];  % size of square stimuli
goout = 16;          % how many on each side
spacing = res/32;   % gap in-between
transzone = 2;      % transition zone size in pixels

% calc
c = (1+res)/2;      % what is the center
grid = c-goout*spacing : spacing : c+goout*spacing;  % grid of where the square centers live

% construct square stimuli
xx = []; yy = [];
stimulus = zeros(res,res,length(grid),length(grid),length(sizes));
for p=1:length(sizes), p
  for q=1:length(grid), q
    for r=1:length(grid)
      [im,xx,yy] = makesquareimage(res,grid(q),grid(r),sizes(p),transzone,xx,yy);
      stimulus(:,:,q,r,p) = im;
    end
  end
end

% % inspect stimuli
% viewmovie(reshape(permute(squish(permute(stimulus,[3 4 5 1 2]),3),[2 3 1]),res,res,1,[]),[],[],[0 1]);

% make the PRF model
model0 = makegaussian2d(res,[],[],sd,sd);

% the response
rfun = @(expt,ft) reshape((squish(permute(squish(permute(ft*stimulus,[3 4 5 1 2]),3),[2 3 1]),2)' * vflatten(l1unitlength(makegaussian2d(res,[],[],sd*sqrt(expt),sd*sqrt(expt))))) .^ expt,[length(grid) length(grid) length(sizes)]);

% visualize model
figureprep([100 100 150 150]); hold on;
set(gca,'YDir','reverse');
axis([.5 res+.5 .5 res+.5]); axis equal tight; axis([.5 res+.5 .5 res+.5]);
imagesc(model0,[0 1]);
%%drawellipse(c,c,0,2*sd,2*sd,[],[],'w--');
set(gca,'XTick',[.5 (1+res)/2 res+.5],'XTickLabel',mat2cellstr([-.5 0 .5]));
set(gca,'YTick',[.5 (1+res)/2 res+.5],'YTickLabel',mat2cellstr([.5 0 -.5]));
figurewrite('prfmodel',[],-1,'~/inout',0);

% visualize stimuli
for p=1:length(sizes)
  figureprep([100 100 150 150]); hold on;
  set(gca,'YDir','reverse');
  axis([.5 res+.5 .5 res+.5]); axis equal tight; axis([.5 res+.5 .5 res+.5]);
  imagesc(stimulus(:,:,(end+1)/2,(end+1)/2,p),[0 1]);
  [xxx,yyy] = meshgrid(grid,grid);
  set(scatter(xxx(:),yyy(:),25,'b.'),'CData',[.5 .5 .5]);
  set(gca,'XTick',[.5 (1+res)/2 res+.5],'XTickLabel',mat2cellstr([-.5 0 .5]));
  set(gca,'YTick',[.5 (1+res)/2 res+.5],'YTickLabel',mat2cellstr([.5 0 -.5]));
  figurewrite('stimulus%d',p,-1,'~/inout');
end
    % for p=1:length(sizes)
    %   set(drawrectangle(c,c,sizes(p),[],'-'),'Color',rand(1,3));
    % end

% do it
ees = [1 .2 .01];   %[0.125];%[1 0.5 0.25 0.125];
for p=1:length(ees)
  temp = feval(rfun,ees(p),1);
  figureprep([100 100 700 150]);

  subplot(1,3,1);
  imagesc(makeimagestack(temp,[0 1],[],-1));
  axis equal tight; colormap(gray); axis off;

  subplot(1,3,2); hold on;
  h = [];
  for q=1:length(sizes)
    h(q) = plot(grid,temp((end+1)/2,:,q),'-','Color',subscript(gray(length(sizes)+1),{q ':'}));
  end
%  legend(h,'Location','NorthEastOutside');
%%    straightline([c-2*sd c+2*sd],'v','k--');
  axis([.5 res+.5 0 1.2]);
  set(gca,'XTick',[.5 (1+res)/2 res+.5],'XTickLabel',mat2cellstr([-.5 0 .5]));
  set(gca,'YTick',[0 .5 1]);
  xlabel('x-position');
  ylabel('Response');

  subplot(1,3,3); hold on;
  bar(flatten(temp((end+1)/2,(end+1)/2,:)));
  axis([0 length(sizes)+1 0 1.2]);
  set(gca,'XTick',1:length(sizes));
  set(gca,'YTick',[0 .5 1]);
  xlabel('Size');
  ylabel('Response');

  figurewrite('big%d',p,-1,'~/inout');

end

    % final = [];
    % for p=1:length(sizes)
    %   wt = bsxfun(@times,stimulus(:,:,:,:,p),reshape(resp(:,:,p),1,1,length(grid),length(grid)));  % res x res x grid x grid
    %   im = sum(sum(wt,3),4);
    %   den = sum(sum(stimulus(:,:,:,:,p),3),4);
    %   final(:,:,p) = im./den;
    % end
    % 
    % figure; imagesc(makeimagestack(final));

%%%%%

notes:
- delete horizontal bar on bar chart
- 200% for big image
- 90% for top row
- avoid crowding the 0 (five clicks up)
- regular hyphen
- bring x-axis to front
- italicize n
- use , not parentheses
* convert all those images to grayscale!
