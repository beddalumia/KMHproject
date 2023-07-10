M = 0.5
t = 1
l = 0.01
f = pi/2

H = zeros(6,6);

for i = 1:6
   H(i,i)   = (-1)^(i+1) * M;
end

H = H - t * diag(ones(1,6-1),+1);
H = H - t * diag(ones(1,6-1),-1);
H = H + l * diag(ones(1,6-2),+2) * exp(+1i*f);
H = H + l * diag(ones(1,6-2),-2) * exp(-1i*f);

H(6,1) = -t;
H(1,6) = -t;

H(5,1) = +l*exp(+1i*f);
H(1,5) = +l*exp(-1i*f);

H(6,2) = +l*exp(+1i*f);
H(2,6) = +l*exp(-1i*f);

disp Re[H]:
disp(real(H))
disp Im[H]:
disp(imag(H))

[V,E] = eig(H);
E = diag(E);