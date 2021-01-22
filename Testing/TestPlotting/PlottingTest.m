clear all
clc

%% Single Phase-Lines

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

%% Full Phase Diagram | Just one channel

varID = 10;  % \in [1,15]

[SOI_list, SOI_names] = get_list('SOI');
Nlines = length(SOI_list);
phaseVAR = cell(Nlines,1);
transLine = zeros(2,Nlines);
for iSOI = 1:Nlines
    lineID = SOI_names(iSOI);
    cd(lineID);
    clear('ids','obs','U_list');
    load('observables_line.mat','ids','obs','U_list');
    if iSOI==1
        figure("Name",ids{varID});
        ax = axes;
    end
    % Get the line data
    U = U_list;
    SOI = SOI_list(iSOI)*ones(length(U),1);
    phaseVAR{iSOI} = obs{varID};
    z = phaseVAR{iSOI};
    ztrim = z(z<1e-2);
    ztrans = max(ztrim);
    transID = find(z==ztrans);
    transLine(2,iSOI) = U(transID);
    transLine(1,iSOI) = SOI(transID);
    Sct = scatter(ax,SOI,U,30,z,'filled','MarkerFaceAlpha',1); hold on
    %Plt = plot3(ax,SOI,U,z,'g','LineWidth',2); hold on
    cd('..');
end
% Plot transition line
plot(transLine(1,:),transLine(2,:),'r','LineWidth',2.5);
% Title, legend, all of that
title(ax,ids{varID});
xlabel(ax,'\lambda_{SO} / t');
ylabel(ax,'U / t');
ax.Box = 'on';
colormap(ax,'copper');
cb = colorbar(ax);
%view(ax,-70,52);
fig = gcf;
fig.Renderer='Painters';
clc

%% Full Phase Diagram | Spin-Resolved

varID = 'ImSigma';  % 'Z' | 'ImSigma' | 'spinDensity'

[SOI_list, SOI_names] = get_list('SOI');
Nlines = length(SOI_list);
phaseVARup = cell(Nlines,1); phaseVARdown = phaseVARup;
figure("Name",varID);
% Axes for spinUP
axUP = axes;
for iSOI = 1:Nlines
    lineID = SOI_names(iSOI);
    cd(lineID);
    clear('ids','obs','U_list');
    load('observables_line.mat','ids','obs','U_list');
    if strcmp(varID,'Z')==1
        iOBS = 10;
    elseif strcmp(varID,'ImSigma')==1
        iOBS = 12;
    elseif strcmp(varID,'spinDensity')==1
        iOBS = 3;
    else
        error('Invalid VAR!');
    end
    % Get the line data
    U = U_list;
    SOI = SOI_list(iSOI)*ones(length(U),1);
    phaseVARup{iSOI} = obs{iOBS};
    colorUP = abs(phaseVARup{iSOI});
    Sup = scatter(axUP,SOI,U,100,colorUP,'filled','MarkerFaceAlpha',1);
    hold on
    cd('..');
end
% Axes for spinDW
axDW = axes;
for iSOI = 1:Nlines
    lineID = SOI_names(iSOI);
    cd(lineID);
    clear('ids','obs','U_list');
    load('observables_line.mat','ids','obs','U_list');
    if strcmp(varID,'Z')==1
        iOBS = 10;
    elseif strcmp(varID,'ImSigma')==1
        iOBS = 12;
    elseif strcmp(varID,'spinDensity')==1
        iOBS = 3;
    else
        error('Invalid VAR!');
    end
    % Get the line data
    U = U_list;
    SOI = SOI_list(iSOI)*ones(length(U),1);
    phaseVARdown{iSOI} = obs{iOBS+1};
    colorDW = abs(phaseVARdown{iSOI});
    Sdw = scatter(axDW,SOI,U,100,colorDW,'filled','MarkerFaceAlpha',.37);
    hold on
    cd('..');
end
% Link axes together
linkaxes([axUP,axDW])
% Hide the top axes
axDW.Visible = 'off';
axDW.XTick = [];
axDW.YTick = [];
% Give each one its own colormap
colormap(axUP,'autumn')
colormap(axDW,'winter')
% Then add colorbars and get everything lined up
set([axUP,axDW],'Position',[.17 .11 .685 .815]);
cbUP = colorbar(axUP,'Position',[.04 .11 .0675 .815]);
cbDW = colorbar(axDW,'Position',[.88 .11 .0675 .815]);
% Title, legend, all of that
title(axUP,varID);
xlabel(axUP,'\lambda_{SO} / t');
ylabel(axUP,'U / t');
axUP.Box = 'on';
clc

