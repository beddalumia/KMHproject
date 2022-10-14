driver = 'ed_km_flakes';
doMPI = true;

% SETUP ARRAY OF FLAKE-SIZE VALUES
flake_size = 1:1:10;
Nsize = length(flake_size);

% GET SLURM ARRAY TASK ID (if undefined would be NaN)
aID  = str2double(getenv('SLURM_ARRAY_TASK_ID'));

if isnan(aID)
    % COMPUTE ALL LINES -> 1:Nsize
    iSTART = 1;
    iSTOP  = Nsize;
else
    % COMPUTE SIZE(aID) ONLY -> aID:aID
    iSTART = aID;
    iSTOP  = aID;
end

% Evil trick
if aID == 1
    jID = str2double(getenv('SLURM_ARRAY_JOB_ID'));
    sbatching = sprintf('sbatch --dependency=afternotok:%d rdmft.sh',jID)
    system(sbatching);
end

%% R-Line: outer loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSIZE = iSTART:iSTOP

    R  = flake_size(iSIZE);    % Input Spin-Orbit 

    rDIR= sprintf('R%d',R);    % Make a folder named 'R<r>', where '<r>' is
    mkdir(rDIR);               % the given integer "radius" for the flake.
    cd(rDIR);                  % Enter the R-folder

    copyfile ../input*         % Copy inside the **external** input file

    %% Phase-Line: inner loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Uold    = NaN;			% Restart point
    Umin    = 0.00;         % Phase-line iSTART
    Ustep   = 0.10;			% Phase-line step
    Umax    = 10.0;         % Phase-line end

    % Evil trick
    if isfile('U_conv.txt')    
        load U_conv.txt
        U_conv = sort(U_conv);
        Umin    = U_conv(end);  % Where we were...
        Uold    = Umin-Ustep;   % New restart point
    end

    if abs(Umin-Umax) > Ustep/2
        % Aim at completing the line
        runDMFT.dry_line(driver,doMPI,Uold,Umin,Ustep,Umax,'radius',R,'wmixing',0.25);
    else
        if isfile('U_list.txt')
            load U_list.txt
            if length(U_list) == length(U_conv)
                fprintf(2,"\n\nWARNING: We are already done with R=%d line!\n\n",R)
                continue
            end
        end
        % Converge missing points
        fprintf(2,"\n\nWARNING: Ready to refresh missing points in R=%d line:\n\n",R)
        U_todo = setdiff(U_list,U_conv);
        fprintf(2,formattedDisplayText(U_todo));
        for iHubbard = 1:length(U_todo)
            runDMFT.refresh_point(driver,doMPI,100,true,'wmixing',0.1)
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd ..                          % Exit the R-folder

end