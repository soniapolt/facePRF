function scratch_axis
% plot an empty polar axis for figure making

set(gcf,'color',[1 1 1]);%,'Units', 'Normalized', 'OuterPosition', [.2 .2 .7 1]);
axis square;axis([-7 7 -7 7]);
    for m = floor([0:1:5.5])
        a = plotCircle(0,0,m,[0 0 0],.75,'edge'); % center
        set(a,'Linewidth',1); hold on;
    end
    hline(0,'k:'); vline(0,'k:');
    set(gca,'visible','off'); 
end

