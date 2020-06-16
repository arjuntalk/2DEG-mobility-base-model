function mobility = getMobility(x, logNdd, logNbk, logDelta, logDrough)
% getMobility takes in density in units of 1E10/cm^2
% and return mobility in units of 1E6cm^2/Vs
% logNdd = log10(Ndd)
% logNbk = log10(Nbk)
% logDelta = log10(delta)
% logDrough = log10(drough)

    % Convert units to SI for density
    x = x*1E14;

    % Convert back to regular form
    Ndd = 10^logNdd;
    Nbk = 10^logNbk;
    delta = 10^logDelta;
    drough = 10^logDrough;

    depth = 110E-9;
    
    % This and the for loop is to 'vectorize' the code which is needed for
    % the fit function
    mobility = zeros(size(x));

    for ii = 1:length(x)
        density = x(ii);

        % Define constants to be used by a lot of the functions
        const.ee = 1.602E-19; % C
        const.mfree = 9.11E-31; % kg
        const.meff = const.mfree*0.067; % kg
        const.h = 6.626E-34; % J*s
        const.hbar = const.h/(2*pi); % J*s
        const.eps0 = 8.854E-12; % SI units
        const.epsr = 12.8;
        const.eps = const.eps0*const.epsr;

        ww = 20E-9; % Width of electron wave function
        kTF = 2*pi*const.meff*const.ee^2/(const.h^2*const.eps);

        preFactor = const.meff*const.ee^2/(2*const.eps*pi*const.hbar^2);

        kF = (2*pi*density).^(1/2);

        bb = 33*density*const.meff*const.ee^2/(8*const.eps*const.hbar^2);
        bb = bb^(1/3);

        % Calcualte the scattering from background impurities
        xVals = linspace(1E-8,pi,500);
        backAnon = @(xx) background(xx, 2*kF*sin(xx/2), bb, preFactor);
        one_by_tauB = trapz(xVals,arrayfun(backAnon, xVals));
        scaleFactor = const.meff/(const.hbar^3*pi)*...
            (const.ee^2/(2*const.eps))^2*...
            (0.5*Nbk/(8*kF^3));
        one_by_tauB = scaleFactor*one_by_tauB;

        % Calculate the scatter from delta doping
        xVals = linspace(0,pi,500);
        dopAnon = @(xx) doping(2*kF*sin(xx/2), kTF, bb, depth);
        one_by_tauD = trapz(xVals,arrayfun(dopAnon,xVals));
        scaleFactor = const.meff/(const.hbar^3*2*pi)*...
            (const.ee^2/(2*const.eps))^2*Ndd/kF^2;
        one_by_tauD = scaleFactor*one_by_tauD;

        % Calculate mobility from interface roughness
        xVals = linspace(1E-8,2*pi,500);
        interfaceAnon = @(xx) interface(const,2*kF*sin(xx/2),bb,drough);
        one_by_tauI = trapz(xVals,arrayfun(interfaceAnon,xVals));
        scaleFactor = const.meff*drough.^2*delta.^2*const.ee^4*...
            (Nbk*ww + density/2).^2/(const.hbar^3*const.eps^2*4*kF^2);
        one_by_tauI = scaleFactor*one_by_tauI;

    %     one_by_tauB = 0;
    %     one_by_tauD = 0;
    %     one_by_tauI = 0;
        one_by_tau = one_by_tauB + one_by_tauD + one_by_tauI;

        mobility(ii) = const.ee/(const.meff*one_by_tau);
    end
    
    % Convert mobility units
    mobility = mobility/1E2;
end

% Helper function for calculating background scattering
function res = background(xx, qq, bb, preFactor)
    piq = 1; % Holds only for q < 2kf
    Sq = 1 + preFactor*nu(bb,qq)*piq/qq;
    res = (1 - cos(xx))/(Sq^2*(sin(xx/2)^3));
end

% Helper function for calculating doping scattering
function res = doping(qq, kTF, bb, depth)
    res = exp(-2*qq*depth)/((qq + kTF*nu(bb,qq))^2)*(bb/(bb+qq))^6*qq^2;
end

% Helper function for calculating interface scattering
function res = interface(const, qq, bb, dd)
    temp = const.ee^2*const.meff/(2*const.eps*const.hbar^2*pi);
    res = qq^2*exp(-dd.^2*qq^2/4)/((1 + temp*nu(bb,qq)/qq)^2);
end

% Helper function for calculating nu(b,q)
function res = nu(bb, qq)
    res = 1/8*bb*(8*bb^2 + 9*bb*qq + 3*qq^2)*(bb+qq)^-3;
end

% Helper function for calculating dielectric function
function epsq = eps_of_q(const, qq)
    epsq = 1 + const.ee^2*2*const.meff/(4*pi*const.epsr*qq*const.hbar^2)*F(q,0);
end

% Helper function for calculating F(q,z)
function res = F(qq, zz)
    res = 1;
end