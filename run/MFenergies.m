clear all
clc

doBands  = false;    %  true  |  false
doVector = false;    %  true  |  false
doRaster = false;    %  true  |  false

% Select mode ('line' | 'map')
mode = 'map';

% Dirty path selector
DATA = '../../Data/KMH-MF_Data/';
cd(DATA)

whichMF = {'AFMz','AFMxy','AFMxyz'};

for iMF = 1:2

    cd(whichMF{iMF});

    [SOI_list, SOI_names] = postDMFT.get_list('SOI');
    Nlines = length(SOI_list);

    for iSOI = 1:Nlines
        SOI = SOI_list(iSOI);
        lineID = SOI_names(iSOI);
        fprintf('************\n');
        fprintf(lineID);
        fprintf('\n************\n');
        cd(lineID);
        clear('ids','ordpms','U_list');
        load('order_parameter_line.mat','ids','ordpms','U_list');
        Emfb(iSOI,:) = load('mf_energy.txt');      % mean-field bandstructure
        try
        Ekin(iSOI,:) = load('kinetic_energy.txt'); % dmft-like kinetic energy
        catch
        Ekin(iSOI,:) = NaN;
        end
        Npar = length(ordpms);
        Npoints = length(U_list);
        for iHubb = 1:Npoints
            U = U_list(iHubb);
            UDIR= sprintf('U=%f',U);
            cd(UDIR);
            %
            Emag(iSOI,iHubb) = 0;
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
            Ne = 2; % --> In the *unit-cell*
            %
            if      Npar == 2
                Sz = ordpms{1}(iHubb);
                Rz = ordpms{2}(iHubb);
            elseif  Npar == 4
                Sx = ordpms{1}(iHubb);
                Sy = ordpms{2}(iHubb);
                Rx = ordpms{3}(iHubb);
                Ry = ordpms{4}(iHubb);
            elseif  Npar == 6
                Sx = ordpms{1}(iHubb);
                Sy = ordpms{2}(iHubb);
                Sz = ordpms{3}(iHubb);
                Rx = ordpms{4}(iHubb);
                Ry = ordpms{5}(iHubb);
                Rz = ordpms{6}(iHubb);
            end
            %
            fprintf('Computing "< >< >" term for U=%.2f ..',U);
            Emag(iSOI,iHubb) = Emag(iSOI,iHubb) + (Sz+Rz)^2 - (Pz+Ne)^2 + (Sz-Rz)^2 - (Pz-Ne)^2;
            Emag(iSOI,iHubb) = Emag(iSOI,iHubb) + (Sx+Rx)^2 + (Sy+Ry)^2 + (Sx-Rx)^2 + (Sy-Ry)^2;
            Emag(iSOI,iHubb) = Emag(iSOI,iHubb) * U/32;
            Emag(iSOI,iHubb) = Emag(iSOI,iHubb) + U/04; % Hartree-Fock correction (<-> Ne=0)
            %
            fprintf('.DONE!\n');
            %% Plotting Bands
            %------------------------------------------------------------------
            if doBands
                Eigenbands = load('EigenbandsKM.mf');
                id = sprintf('Mean-Field Bands [SOI=%f U=%f]',SOI,U);
                figure("Name",id);
                Ncell = length(Eigenbands);
                Ev = Eigenbands(1:floor(Ncell/2),:);
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
                % PRINTING (if requested)
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
            cd('..');
        end
        cd('..');

        Etot(iSOI,:) = Emfb(iSOI,:) + Emag(iSOI,:); % mean-field bands + order parameters
        Epot(iSOI,:) = Etot(iSOI,:) - Ekin(iSOI,:); % should correspond to <Uloc*nup*ndw>

        if strcmp(mode,'line')
            figure('Name',whichMF{iMF})
            plot(U_list,Etot(iSOI,:))
            hold on
            plot(U_list,Ekin(iSOI,:))
            plot(U_list,Epot(iSOI,:))
        end

    end

    cd('..')

    if strcmp(mode,'map')
        % Get the map data
        [X,Y] = meshgrid(SOI_list,U_list); % Assaad's convention
        Z = Etot'; 
        % Plot the map data
        figure('Name',whichMF{iMF})
        switch whichMF{iMF}
            case 'AFMxy'
                AFMxy = Z;
                surf(X,Y,AFMxy,'FaceAlpha',0.5);
                colormap(plotDMFT.colorlab.matplotlib.plasma)
            case 'AFMz'
                AFMz  = Z;
                surf(X,Y,AFMz,'FaceAlpha',0.5);
                colormap(plotDMFT.colorlab.matplotlib.viridis)
            otherwise
                AFM   = Z;
                surf(X,Y,AFM,'FaceAlpha',0.5);
                colormap(plotDMFT.colorlab.matplotlib.cividis)
        end
        view(136,30)
    end

end

if size(AFMxy) == size(AFMz)
    figure('Name','GSenergy difference @ Hartree-Fock')
    title('Groundstate Energy Difference [AFMxy - AFMz]','Interpreter','latex')
    xlabel('$U/t$','Interpreter','latex')
    ylabel('$\lambda_{SO}/t$','Interpreter','latex')
    meshz(X,Y,AFMxy-AFMz,'EdgeColor','k'); view(136,30)
end

% Dirty path reset
CODE = '../../KMproj[git]/run/';
cd(CODE);


