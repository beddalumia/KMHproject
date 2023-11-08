TOPLEVEL = pwd;

[~,rdirs] = QcmP.post.get_list('R');

for r = 1:length(rdirs)

    fprintf("    > %s\n",rdirs(r))
    cd(rdirs(r))

    [~,udirs] = QcmP.post.get_list('U');

    for u = 1:length(udirs)

        fprintf("        >> %s\n",udirs(u))
        cd(udirs(u))

        build_spin_coordinates(TOPLEVEL)

        cd('..')

    end

    cd('..')

end

function build_spin_coordinates(topdir)

    try
        flake = load('flake.txt');
    catch
        warning("'flake.dat' not found in %s",pwd)
        cd(topdir)
        error("Stopping here, in the top directory: %s",pwd)
    end

    spin_flake = zeros(size([flake;flake]));

    spin_flake(1:2:end,:) = flake;
    spin_flake(2:2:end,:) = flake;

    writematrix(spin_flake,'spin_flake.dat');
    writematrix(spin_flake(:,1),'spin_xvals.dat');
    writematrix(spin_flake(:,2),'spin_yvals.dat');
    disp("            spin coordinates have been built!")

end