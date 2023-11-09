TOPLEVEL = pwd;

[~,rdirs] = QcmP.post.get_list('R');

for r = 1:length(rdirs)

    cd(rdirs(r))

    [~,udirs] = QcmP.post.get_list('U');

    for u = 1:length(udirs)

        fprintf("> %s/%s\n",rdirs(r),udirs(u))
        cd(udirs(u))

        try
            load('Z2marker.txt')
        catch
            Z2marker = NaN;
        end
        Z2bulk = mode(round(Z2marker));
        fprintf(">> Z2bulk = %d\n",Z2bulk)
        writematrix(Z2bulk);

        cd('..')

    end

    cd('..')

end