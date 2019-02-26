function SOLID_GUI_saveManualLabels(oname, modZ, modZorig)
    modZdiffs = modZ ~= modZorig;            
    if any(modZdiffs(:))
        fid = fopen(strcat(oname, '_modZ2D_manual_edits.txt'), 'w');
        for i = 1:size(modZdiffs, 1)
            fprintf(fid, repmat('%f, ', [1, size(modZdiffs,2)-1]), modZdiffs(i, 1:end-1));
            fprintf(fid, '%f\n', modZdiffs(i, end));
        end
        fclose(fid);
        fid = fopen(strcat(oname, '_modZ2D_automatic.txt'), 'w');
        for i = 1:size(modZorig, 1)
            fprintf(fid, repmat('%f, ', [1, size(modZorig,2)-1]), modZorig(i, 1:end-1));
            fprintf(fid, '%f\n', modZorig(i, end));
        end
        fclose(fid); 
    end
end