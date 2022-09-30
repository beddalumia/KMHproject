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

%% SOC-Line: outer loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSIZE = iSTART:iSTOP

    R  = flake_size(iSIZE);    % Input Spin-Orbit 

    rDIR= sprintf('R=%f',R);   % Make a folder named 'R=...', where '...'
    mkdir(rDIR);               % is the given value for SpinOrb interaction
    cd(rDIR);                  % Enter the R-folder

    copyfile ../input*         % Copy inside the **external** input file

    %% Phase-Line: inner loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Uold    = NaN;			% Restart point
    Umin    = 0.00;         % Phase-line iSTART
    Ustep   = 0.10;			% Phase-line step
    Umax    = 10.0;         % Phase-line end

    runDMFT.dry_line(driver,doMPI,Uold,Umin,Ustep,Umax,'radius',R,'wmixing',0.25);

    %runDMFT.refresh_line(driver,doMPI,100,true,'bath_type','replica','wmixing',0.1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd ..                          % Exit the R-folder

end