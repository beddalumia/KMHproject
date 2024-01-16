clear,clc

%% INPUT

doSmooth = 0;    %  true  |  false  -->  Smoothing on energy differences
doBands  = 0;    %  true  |  false  -->  plot mean-field bandstructure
doCheck  = 1;    %  true  |  false  -->  check formula for < >< > terms

mode = 'map';        % 'line' | 'map'   -->  many lines or phase-diagram

whichMF = {'AFMz','AFMxy','AFMxyz'}; % which mean-field decoupling to look

%% MAIN

% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/MF/';
cd(DATA)

for iMF = 1:2

    cd(whichMF{iMF});
    disp(whichMF{iMF});

    [SOI_list, SOI_names] = QcmP.post.get_list('SOI');
    Nlines = length(SOI_list);

    for iSOI = 1:Nlines
        SOI = SOI_list(iSOI);
        lineID = SOI_names(iSOI);
        fprintf('••••••••••••\n');
        fprintf(lineID{:});
        fprintf('\n••••••••••••\n');
        cd(lineID{:});
        Emfb(iSOI,:) = load('mf_bands_energy.txt');  % mean-field bandstructure
        Ekin(iSOI,:) = load('kinetic_energy.txt');   % dmft-like kinetic energy
        Epot(iSOI,:) = load('potential_energy.txt'); % matsubara potental energy
        Emag(iSOI,:) = load('magnetic_energy.txt');  % magnetic ordering energy
        Etot(iSOI,:) = Emfb(iSOI,:) + Emag(iSOI,:);  % mean-field bands + order parameters
        Upot(iSOI,:) = Etot(iSOI,:) - Ekin(iSOI,:);  % should correspond to <Uloc*nup*ndw>

        if strcmp(mode,'line')
            figure('Name',whichMF{iMF})
            plot(U_list,Etot(iSOI,:))
            hold on
            plot(U_list,Ekin(iSOI,:))
            plot(U_list,Epot(iSOI,:))
        end

        if doBands
            plot_bands(SOI)
        end

        if doCheck
            check_magnetic_energy(Emag(iSOI,:));
        end

        U_list = QcmP.post.get_list('U');

        cd('..');

    end

    cd('..')

    if strcmp(mode,'map')
        % Get the map data
        [X,Y] = meshgrid(SOI_list,U_list); % Assaad's convention
        Z = Etot'; 
        if doSmooth, Z = smoothdata(Z,'SmoothingFactor',0.05); end
        % Plot the map data
        QcmP.plot.import_colorlab
        figure('Name',whichMF{iMF})
        switch whichMF{iMF}
            case 'AFMxy'
                AFMxy = Z;
                surf(X,Y,AFMxy,'FaceAlpha',0.5);
                colormap(palette.cmocean('-matter'))
            case 'AFMz'
                AFMz  = Z;
                surf(X,Y,AFMz,'FaceAlpha',0.5);
                colormap(palette.cmocean('-algae'))
            otherwise
                AFM   = Z;
                surf(X,Y,AFM,'FaceAlpha',0.5);
                colormap(palette.cmocean('tempo'))
        end
        view(136,30)
        close
    end

end

if strcmp(mode,'map') & size(AFMxy) == size(AFMz)
    figure('Name','GSenergy difference @ Hartree-Fock')
    title('Groundstate Energy Difference [AFMxy - AFMz]','Interpreter','latex')
    xlabel('$U/t$','Interpreter','latex')
    ylabel('$\lambda_{SO}/t$','Interpreter','latex')
    meshz(X,Y,AFMxy-AFMz,'EdgeColor','none','FaceColor','interp');%,str2rgb('powder blue'));
    palette.cmocean('matter');
    view(110,30)
    xlim([0,0.3]);
    xlabel('$\lambda_\mathrm{so}/t$','Interpreter','latex');
    ylim([0,10]);
    ylabel('$U/t$','Interpreter','latex');
    zlim([-0.0915,0]);
    zlabel('$E_\parallel - E_\perp$','Interpreter','latex');
    caxis([-0.0915,0]);
    grid off
%     figure('Name','GSenergy difference @ Hartree-Fock')
%     title('Groundstate Energy Difference [AFMxy - AFMz]','Interpreter','latex')
%     xlabel('$U/t$','Interpreter','latex')
%     ylabel('$\lambda_{SO}/t$','Interpreter','latex')
%     imagescn(SOI_list,U_list,AFMxy-AFMz);
%     palette.cmocean('matter');
%     colorbar
%     caxis([-0.1,0])
end

% Dirty path reset
cd(CODE);

% Export to LaTeX
addpath ../lib/m2tex/src
matlab2tikz('filename','../fig/MF_energy.tex','width','5cm','heigth','4cm');
rmpath ../lib/m2tex/src

%% contains

function plot_bands(SOI)
 doVector = false;    %  true  |  false  -->  save bands to vectorized file
 doRaster = false;    %  true  |  false  -->  save bands to rasterized file
    U_list = postDMFT.get_list('U');
    Npoints = length(U_list);
    for iHubb = 1:Npoints
        U = U_list(iHubb);
        UDIR= sprintf('U=%f',U);
        cd(UDIR);
        %% Plotting Bands
        Eigenbands = load('EigenbandsKMH.mf');
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
        cd('..');
    end
end

function Emag = check_magnetic_energy(Emag_ref)
    load('order_parameter_line.mat','ids','ordpms','U_list');
    Npar = length(ordpms);
    Npoints = length(U_list);
    for iHubb = 1:Npoints
        U = U_list(iHubb);
        UDIR= sprintf('U=%f',U);
        cd(UDIR);
        %
        Emag(iHubb) = 0;
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
        Emag(iHubb) = Emag(iHubb) + (Sz+Rz)^2 - (Pz+Ne)^2 + (Sz-Rz)^2 - (Pz-Ne)^2;
        Emag(iHubb) = Emag(iHubb) + (Sx+Rx)^2 + (Sy+Ry)^2 + (Sx-Rx)^2 + (Sy-Ry)^2;
        Emag(iHubb) = Emag(iHubb) * U/32;
        Emag(iHubb) = Emag(iHubb) + U/04; % Hartree-Fock correction (<-> Ne=0)
        %
        fprintf('.DONE!\n');
        %
        fprintf(' REF = %f\n DIFF = %f\n',Emag_ref(iHubb),Emag(iHubb)-Emag_ref(iHubb));
        %
        cd('..');
    end
end