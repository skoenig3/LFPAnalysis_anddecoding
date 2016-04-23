function dY = odeFcn(t,y,m,k,c,Ft,time)
%differnetial equation for harmonic dampened oscillator
%t: time
%y: position/voltage
%w:harmonic frequency (omega) = sqrt(k/m)
%L: dampening coefficient/ratio (L < 1 for overdamppened, L = 1 for cirtically
%dampened, L > 1 for underdamppenend). L = c/[2*sqrt(m*k)]
%Ft: drive over time

FT = interp1(time,Ft,t); % Interpolate the data set (ft, ff) at times t

% if length(m) == 1
    eq = -1/m * (k * y(1) + c * y(2) - FT);
    dY = [y(2);eq];
% else %springs in series
%     eq1 = -1/m(1) * (k(1) *  y(1) + c(1) * y(2)  - k(2) * y(4) - FT);
%     eq2 = -1/m(2) * (-k(2) * y(1) + c(2) * y(4)  + k(2) *  y(4));
%     dY = [y(2); eq1; y(3); eq2];
% end
    
end