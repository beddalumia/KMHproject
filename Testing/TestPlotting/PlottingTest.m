clear all
clc

        %%% Structure (for ed_kane_mele) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obs{1}    Density (at half filling = 1)                   %
        % obs{2}    Docc (Double Occupancy \in [0,1])               %
        % obs{3}    Nup (intensive)                                 %
        % obs{4}    Ndown (intensive)                               %
        % obs{5}    Magn(Nup-Ndown)                                 %
        % obs{6}    S2 (Impurity magnetic dipole: <S^2>)            %
        % obs{7}    Egs (Ground-State Energy)                       %
        % obs{8}    Sz2_11 (Impurity magnetic dipole:<Sz^2>)        %
        % obs{9}    N2_11 (??????)                                  %
        % obs{10}   Z: Quasiparticle Weight for spin-up             %
        % obs{11}   Z: Quasiparticle Weight for spin-down           %
        % obs{12}   Sig1_s1(\Im\Sigma(iw=0) for spin-up)            %
        % obs{13}   Sig1_s1(\Im\Sigma(iw=0) for spin-down)          %
        % obs{14}   Nph (??????)                                    %
        % obs{15}   Wph (??????)                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Select observable (iOBS==0 means everything)

        iOBS = 11;

% We don't have a SOI-values list, but we can obtain that by just
% inspecting the subdirectories...
[SOI_list, SOI_names] = get_list('SOI');
Nlines = length(SOI_list);
for iSOI = 1:Nlines
    lineID = SOI_names(iSOI);
    cd(lineID); fprintf(lineID);
    clear('ids','obs','U_list');
    load('observables_line.mat','ids','obs','U_list');
    if iOBS~=0
        figName = strcat(lineID,' | ',ids{iOBS})';
        figure("Name",figName);
        plot(U_list,obs{iOBS},'LineWidth',2); drawnow
    else
        for iOBS = 1:length(ids)
            figName = strcat(lineID,' | ',ids{iOBS})';
            figure("Name",figName);
            plot(U_list,obs{iOBS},'LineWidth',2); drawnow
        end
    end
    cd('..'); fprintf('..DONE!\n\nPress any key to continue\n\n'); pause
end

