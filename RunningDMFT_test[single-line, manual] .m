% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Phase Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Here = pwd;                    % Retrieve and store current directory path

SOI  = 0;                      % Input Spin-Orbit 
Umax = 1; Ustep   = 0.25;      % Input Hubbard 

notConverged = 0;		       % Convergence-fail *counter*
notConvThreshold = 5;          % Maximum #{times} we accept DMFT to fail

U = 0; Uold = -1;
doUpdate = true;
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

oldDIR=sprintf('../U=%f',Uold);% -----------------------------------------
if isfolder(oldDIR)            % If it exist a "previous" folder: 
copyfile(oldDIR);              % Copy everything from last dmft evaluation
end                            % -----------------------------------------

copyfile ../inputKANEMELE.conf % Copy inside the *external* input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
driver = 'ed_kane_mele';
HUBBARD =sprintf(' uloc=%f',U);		% OVERRIDE of
T2 =sprintf(' t2=%f',SOI);			% PARAMETERS
OutLOG = ' | tee LOG.out';
dmft_ed_call = [driver,HUBBARD,T2,OutLOG];
ed_status = system(dmft_ed_call)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE WE NEED TO CATCH SOMEHOW A FAILED (unconverged) DMFT LOOP
if isfile('ERROR.README') 
    notConverged = notConverged + 1;
    delete ERROR.README
    ed_status = -1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder


if notConverged > notConvThreshold
   error('DMFT not converged: phase-span stops now!')
elseif notConverged > 0 && notConverged < notConvThreshold
   % You either reduce the U grain...or even do nothing 
   % and manage *the mixing* manually, by modifying
   % the **external** input-file at runtime (MIT still far!)
   % Ustep = Ustep/2; 		   % EITHER this
   doUpdate = false;              % OR this
end

if doUpdate
   Uold = U;
   U = U + Ustep;              % Hubbard update  
else
   U = Uold;
   U = U + Ustep;		% old-Hubbard update (if nonconverged!)
end

end                            % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

