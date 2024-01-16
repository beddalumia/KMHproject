clear,clc

%% INPUT

ignConv = true;    %  true  |  false  ->  Ignores converged-U lists

whichMF = {'AFMz','AFMxy','AFMxyz'}; % which mean-field decoupling

% Dirty path selector
CODE = fileparts(mfilename('fullpath'));
DATA = '../../Data/MF/';
cd(DATA);

for iMF = 1:2

    cd(whichMF{iMF});
    disp(whichMF{iMF});

    % We don't have a SOI-values list, but we can obtain that by just
    % inspecting the subdirectories...
    [SOI_list, SOI_names] = QcmP.post.get_list('SOI');
    Nlines = length(SOI_list);
    for iSOI = 1:Nlines
        SOIDIR = SOI_names(iSOI);
        fprintf(strcat(SOIDIR));
        cd(SOIDIR);
        if isfile('U_list.txt') && not(ignConv)
            U_list = load('U_list.txt');
        else
            U_list = QcmP.post.get_list('U');
        end
        [ids,ordpms] = QcmP.post.order_parameter_line('U',U_list); fprintf('...DONE\n');
        %save('order_parameter_line.mat','ids','ordpms','U_list');
        QcmP.post.custom_line('mf_bands_energy.dat','U',U_list);
        QcmP.post.custom_line('potential_energy.dat','U',U_list);
        QcmP.post.custom_line('magnetic_energy.dat','U',U_list);
        try
            QcmP.post.kinetic_line('U'); % Builds the 'kinetic_energy.txt' file
        end
        cd('..');
    end

    cd('..');

end

% Dirty path reset
cd(CODE);




