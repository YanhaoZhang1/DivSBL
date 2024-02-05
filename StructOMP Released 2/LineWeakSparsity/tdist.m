function ft=tdist(x0, n, v);

% ft=(1+t.^2/v).^(-(v+1)/2);

% ft=pi*(1+t.^2);
% ft=1./ft;
% ft=ft.^(1/2);
t=1:n;
ft=(2+(t-x0).^2).^v;
ft=1./ft;

return