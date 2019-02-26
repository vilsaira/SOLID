function SOLID_savePNG(oname, modZ, b, thresholdLower, thresholdUpper)    
    uniqueb = unique(b);
    for i = 1 : length(uniqueb)
        shell = (b == uniqueb(i));
        modZt = modZ(:,shell);
        modZt(isnan(modZt)) = 0;
        fig = figure('Name', 'SOLID', 'visible','off');
        ax = axes;
        im = imagesc(modZt);
        xlabel('Image volume #');
        ylabel('Slice #');
        cb = colorbar;
        ylabel(cb, 'Modified Z-score');
        title(['SOLID, b-value' num2str(uniqueb(i))]);
        ax.FontSize = 14;
        ax.CLim = [thresholdLower, thresholdUpper];
        ax.YDir = 'normal';
        axis equal tight;
        print(fig, strcat(oname, '_modZ2D_b', num2str(uniqueb(i)), '.png'), '-dpng', '-r300');
        close(fig);
    end
end