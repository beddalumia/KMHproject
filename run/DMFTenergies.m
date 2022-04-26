%% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

global ignUlist
ignUlist = false; %  true  |  false  

[SOI_list, SOI_names] = postDMFT.get_list('SOI');
Nlines = length(SOI_list);
totalEnergies = cell(Nlines,1);
for iSOI = 1:Nlines
    SOI = SOI_list(iSOI);
    lineID = SOI_names(iSOI);
    fprintf('************\n');
    fprintf(lineID);
    fprintf('\n************\n');
    cd(lineID);
    if isfile('U_list.txt')
        U_list = load('U_list.txt');
    else
        U_list = [];
    end
    fprintf('> E_kin.');
    [kins,U_list] = postDMFT.kinetic_line(U_list); fprintf('..DONE\n');
    fprintf('> E_pot.');
    [ids,ens,U_list] = postDMFT.energy_line(U_list); fprintf('..DONE\n');
    pots = ens{1}; % The total potential energy is always the 1st value!
    totalEnergies{iSOI} = kins + pots;
    save('energy_line.mat','ids','ens','kins','U_list');
    % FIGURE ..............................................................
    if iSOI==1
        figure("Name",'DMFT Energy');
        ax = axes;
    end
    % Get the map data
    U = U_list;
    SOI = SOI_list(iSOI)*ones(length(U),1);
    z = totalEnergies{iSOI};
    % Plot the map
    %Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
    Plt = plot3(ax,SOI,U,z,'g','LineWidth',2); hold on
    % .....................................................................
    cd('..');
end

% Title, legend, all of that
title(ax,'DMFT Energy');
xlabel(ax,'\lambda_{SO} / t');
ylabel(ax,'U / t');
ax.Box = 'on';
%colormap(ax,'copper');
%cb = colorbar(ax);
view(ax,-70,52);
fig = gcf;
fig.Renderer='Painters';
%set(gca, 'ZScale', 'log');

