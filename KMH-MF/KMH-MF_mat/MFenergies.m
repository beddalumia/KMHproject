clear all
clc

whichMF = 'AFMxy';  % 'AFMz' | 'AFMxy' | 'AFMxyz'
doBands = true;    %  true  |  false
doVector = false;   %  true  |  false
doRaster = false;   %  true  |  false

DATA = '../../../Data/KMH-MF_Data/';
cd([DATA,whichMF]);

[SOI_list, SOI_names] = postDMFT.get_list('SOI');
Nlines = length(SOI_list);
phaseVAR = cell(Nlines,1);
transLine = zeros(2,Nlines);
for iSOI = 1:1%Nlines
    SOI = SOI_list(iSOI);
    lineID = SOI_names(iSOI);
    fprintf('************\n');
    fprintf(lineID);
    fprintf('\n************\n');
    cd(lineID);
    clear('ids','ordpms','U_list');
    load('order_parameter_line.mat','ids','ordpms','U_list');
    Npar = length(ordpms);
    Npoints = length(U_list);
    for iHubb = 1:Npoints
        U = U_list(iHubb);
        UDIR= sprintf('U=%f',U);
        cd(UDIR);
        %
        Emb(iSOI,iHubb) = 0;
        %
        Sx = 0;
        Sy = 0;
        Sz = 0;
        %
        Rx = 0;
        Ry = 0;
        Rz = 0;
        %
        Pz = 0;
        Nel = 0; % --> Hartree-Fock convention on chemical potential!
        %
        if      Npar == 2
            Sz = ordpms{1}(iHubb)/2;
            Rz = ordpms{2}(iHubb)/2;
        elseif  Npar == 4
            Sx = ordpms{1}(iHubb)/2;
            Sy = ordpms{2}(iHubb)/2;
            Rx = ordpms{3}(iHubb)/2;
            Ry = ordpms{4}(iHubb)/2;
        elseif  Npar == 6
            Sx = ordpms{1}(iHubb)/2;
            Sy = ordpms{2}(iHubb)/2;
            Sz = ordpms{3}(iHubb)/2;
            Rx = ordpms{4}(iHubb)/2;
            Ry = ordpms{5}(iHubb)/2;
            Rz = ordpms{6}(iHubb)/2;
        end
        %
        Emb(iSOI,iHubb) = Emb(iSOI,iHubb) + (Sz+Rz)^2 - (Pz+Nel)^2 + (Sz-Rz)^2 - (Pz-Nel)^2;
        Emb(iSOI,iHubb) = Emb(iSOI,iHubb) + (Sx+Rx)^2 + (Sy+Ry)^2 + (Sx-Rx)^2 + (Sy-Ry)^2;
        Emb(iSOI,iHubb) = Emb(iSOI,iHubb) * U/8;
        %
        Eigenbands = load('Eigenbands.nint');
        cd('..');
        Ncell = length(Eigenbands);
        Ev = Eigenbands(1:floor(Ncell/2),:);
        fprintf('Computing GS energy for U=%f..',U);
        Emf(iSOI,iHubb) = 2*sum(Ev(:,2))/Ncell; % spin-degeneracy x normalization
        fprintf('.DONE!\n');
        %% Plotting Bands
        %------------------------------------------------------------------
        if doBands
            id = sprintf('Mean-Field Bands [SOI=%f U=%f]',SOI,U);
            figure("Name",id);
            Ec = Eigenbands(floor(Ncell/2)+1:end,:);
            scatter(Ev(:,1),Ev(:,2),'r'); hold on
            scatter(Ec(:,1),Ec(:,2),'b');
            % Title, legend, all of that
            title(id);
            xticks([0 2.418399 4.836798 7.245524])
            xlim([0,7.245524]);
            xticklabels({'\Gamma','K','K`','\Gamma'})
            ylabel('\epsilon / t');
            ax = gca; ax.Box = 'on';
        
           % [These two lines ensure filling of the fig]
            InSet = get(ax, 'TightInset');
            set(ax, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3),...
                1-InSet(2)-InSet(4)]);
            if doRaster
               filename = sprintf('MF_bands_SOI=%f_U=%f.png',SOI,U);
               fprintf('Printing %s..',filename);
               print(gcf,filename,'-dpng','-r600');
            elseif doVector
               filename = sprintf('MF_bands_SOI=%f_U=%f.pdf',SOI,U);
               fprintf('Printing %s..',filename);
               print(gcf,filename,'-dpdf','-fillpage');
            else
               fprintf('SKIP printing!\n');
            end
            fprintf('.DONE!\n');
        end
        %------------------------------------------------------------------
    end
    cd('..');
end

CODE = '../../../KMproj[git]/KMH-MF/KMH-MF_mat/';
cd(CODE);

% Get the map data
[X,Y] = meshgrid(SOI_list,U_list);
Z = Emf' + Emb';
% Plot the map data
surf(X,Y,Z);

