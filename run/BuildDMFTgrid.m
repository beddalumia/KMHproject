%% INPUT
driver = 'ed_km_afmxy';
doMPI = true;

% SETUP ARRAY OF SOI VALUES
spin_orbit = linspace(0,0.3,10);
Nsoi  = length(spin_orbit);

% GET SLURM ARRAY TASK ID (if undefined would be NaN)
aID  = str2double(getenv('SLURM_ARRAY_TASK_ID'));

if isnan(aID)
    % COMPUTE ALL SOI LINES -> 1:Nsoi
    iSTART = 1;
    iSTOP  = Nsoi;
else
    % COMPUTE SOI(aID) ONLY -> aID:aID
    iSTART = aID;
    iSTOP  = aID;
end

%% SOC-Line: outer loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSOI = iSTART:iSTOP

    SOI  = spin_orbit(iSOI);       % Input Spin-Orbit 

    soiDIR= sprintf('SOI=%f',SOI); % Make a folder named 'SOI=...', where '...'
    mkdir(soiDIR);                 % is the given value for SpinOrb interaction
    cd(soiDIR);                    % Enter the SOI-folder

    copyfile ../input*             % Copy inside the **external** input file

    %% Phase-Line: inner loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Uold    = NaN;			% Restart point
    Umin    = 0.00;         % Phase-line iSTART
    Ustep   = 0.10;			% Phase-line step
    Umax    = 10.0;         % Phase-line end

    runDMFT.dry_line(driver,doMPI,Uold,Umin,Ustep,Umax,'t2',SOI,'wmixing',0.25);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd ..                          % Exit the SOI-folder

end
