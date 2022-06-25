%% INPUT
driver = 'mf_km_afm';
doMPI = false;  % >> MF-code is not MPI-safe <<

spin_orbit = [0,0.02,0.04,0.06,0.08,0.1,0.2,0.3];
Nsoi = length(spin_orbit);
Nopms   = 4;    % #{order parameters}

%% SOC-Line: outer loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSOI = 1:Nsoi

    SOI  = spin_orbit(iSOI);       % Input Spin-Orbit 

    soiDIR= sprintf('SOI=%f',SOI); % Make a folder named 'SOI=...', where '...'
    mkdir(soiDIR);                 % is the given value for SpinOrb interaction
    cd(soiDIR);                    % Enter the SOI-folder

    copyfile ../input*             % Copy inside the **external** input file

    %% Phase-Line: inner loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Uold    = NaN;			% Restart point
    Umin    = 0.00;         % Phase-line start
    Ustep   = 0.25;			% Phase-line step
    Umax    = 8.00;         % Phase-line end

    runDMFT.dry_line(driver,doMPI,Uold,Umin,Ustep,Umax,'t2',SOI,'nparams',Nopms);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd ..                          % Exit the SOI-folder

end
