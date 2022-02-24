function [V] = V_cyl(Ca)
S = 0.05365; % Stroke
B = 0.06735; % Bore
l = 0.0847; % Connecting rod
cr = 8.5; % Compression ratio

Vd = (pi/4)*B^2*S;    % displacement volume
Vc = Vd/(cr-1);       % Clearance volume

r = S/2;                   % Radius of the crank
x = r*cos(Ca) + sqrt(l^2 - r^2*sin(Ca).^2); % Added/substracted length
d = (l+r-x);

V = pi*d*(B/2)^2 + Vc;
end