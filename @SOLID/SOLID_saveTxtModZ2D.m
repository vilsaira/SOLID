function SOLID_saveTxtModZ2D(oname, modZ)
    fid = fopen(strcat(oname, '_modZ2D.txt'), 'w');
    for i = 1:size(modZ, 1)
        fprintf(fid, repmat('%f, ', [1, size(modZ,2)-1]), modZ(i, 1:end-1));
        fprintf(fid, '%f\n', modZ(i, end));
    end
    fclose(fid);
end