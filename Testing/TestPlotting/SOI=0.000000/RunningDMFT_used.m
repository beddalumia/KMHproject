% Let MATLAB see the goddamn PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -> works only if matlab has been started from a unix terminal! (0^0~~,)
path = getenv('PATH');
path = [path ':/usr/local/bin'];
setenv('PATH', path) 
clear path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Phase Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Here = pwd;                   % Retrieve and store current directory path

SOI  = 0;                     % Input Spin-Orbit 
Umax = 10; Ustep   = 0.01;    % Input Hubbard 

U = 9.51; Uold = 9.5;
while U <= Umax                % Hubbard loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~>

UDIR= sprintf('U=%f',U);       % Make a folder named 'U=...', where '...'
mkdir(UDIR);                   % is the given value for Hubbard interaction
cd(UDIR);                      % Enter the U-folder

if Uold >= 0                   % If it exist a "previous" folder:
oldDIR=sprintf('../U=%f',Uold);% -----------------------------------------
copyfile(oldDIR);              % Copy everything from last dmft evaluation
end                            % -----------------------------------------

copyfile ../inputKANEMELE.conf % Copy inside the external input file

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
dmft_ed_call = sprintf('ed_kane_mele uloc=%f t2=%f | tee LOG.out', U, SOI);
ed_status = system(dmft_ed_call)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..                          % Exit the U-folder

Uold = U;
U = U + Ustep;                 % Hubbard update  

end                           % <~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

