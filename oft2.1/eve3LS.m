function [result, lhs, rhs] = eve3LS(d, neff, lambda, fibreSpec, modeTask)
% Eigen-value equations for three-layer structure modes
%
% Known solutions: Monerie (LP modes), Tsao (TE and TM modes), Erdogan
% (hybrid cladding modes) plus some development solutions.
%
% See also: buildModes, traceMode, findModes, fibreMode, eve2LS

lhs = [];
rhs = [];

assert(nargin == 5,'%s: Invalid number of parameters', mfilename);
assert(~isnan(fibreSpec.coreCladdingRatio));

% Assert mode type for 3-layer structure
switch lower(modeTask.modetype)
    case {'monerie', 'tsaote', 'tsaotm', 'erdogan', 'zhang', 'tsaohybrid'}
        % ok
    otherwise
        error('Wrong mode type');
end

nu = modeTask.modeindex(1);

% convert wavelength to microns to make compatible with size in microns
lambda = lambda / 1000;

k0 = 2*pi ./ lambda;
beta = k0 .* neff;

n1 = refrIndex(fibreSpec.materials{1}, lambda);
n2 = refrIndex(fibreSpec.materials{2}, lambda);
n3 = refrIndex(fibreSpec.materials{3}, lambda);

[IXclad] = find(neff < n2);
[IXcore] = find(neff >= n2);

modetype = char(modeTask.modetype); % conversion from possible cell to string
modetype = modetype(1,:);

if strcmpi(modetype, 'tsaote') || strcmpi(modetype, 'tsaotm') || ...
        strcmpi(modetype, 'tsaohybrid') || strcmpi(modetype, 'zhang')
    r1 = d * fibreSpec.coreCladdingRatio / 2; % core radius;
    r2 = d / 2;
    a = r2 / r1;
    k = 2 * pi ./ lambda;
    
    U1 = sqrt(r1.^2 .* (k^2 * n1^2 - beta.^2) );
    W3 = sqrt( r2.^2 .* (beta.^2 - k^2 * n3^2) );
    V12 = sqrt( k^2 * r1.^2 * (n1^2 - n2^2) );
    V23 = sqrt( k^2 * r2.^2 * (n2^2 - n3^2) );
    
    s21 = n2^2 / n1^2;
    s23 = n2^2 / n3^2;
    
    J = besseljd(nu, U1) ./ (U1 .* besselj(nu, U1));
    K = besselkd(nu, W3) ./ (W3 .* besselk(nu, W3));
elseif strcmpi(modetype, 'monerie') || strcmpi(modetype, 'ari')
    a = d * fibreSpec.coreCladdingRatio / 2; % core radius;
    b = d / 2; % cladding radius
end

switch lower(modetype)
    %% monerie
    case 'monerie' % Monerie 1982
        c = fibreSpec.coreCladdingRatio;
        u = a .* sqrt(k0.^2 .* n1.^2 - beta.^2);
        uprim = d / 2 .* sqrt(k0.^2 .* n2.^2 - beta.^2);
        vprim = d / 2 .* sqrt(beta.^2 - k0.^2 .* n2.^2);
        v = d / 2 .* sqrt(beta.^2 - k0.^2 .* n3.^2);
        
        % Eq. (5)
        Jbar = @(m, x) besselj(m, x) ./ (x .* besselj(m+1, x));
        Ybar = @(m, x) bessely(m, x) ./ (x .* bessely(m+1, x));
        Kbar = @(m, x) besselk(m, x) ./ (x .* besselk(m+1, x));
        Ibar = @(m, x) besseli(m, x) ./ (x .* besseli(m+1, x));
        
        if ~isempty(IXcore)
            % Core mode
            resultCore = (Jbar(nu, u) - Kbar(nu, vprim .* c)) .* (Kbar(nu, v) + Ibar(nu, vprim)) ./ ...
                ((Jbar(nu, u) + Ibar(nu, vprim .* c)) .* (Kbar(nu, v) - Kbar(nu, vprim))) -...
                besseli(nu+1, vprim .* c) .* besselk(nu+1, vprim) ./ (besseli(nu+1, vprim) .* besselk(nu+1, vprim .*c));
            result(IXcore) = resultCore(IXcore);
        end;
        if ~isempty(IXclad)
            % Cladding mode
            resultCladding = (Jbar(nu, u) - Ybar(nu, uprim .* c)) .* (Kbar(nu, v) - Jbar(nu, uprim)) ./ ...
                ((Jbar(nu, u) - Jbar(nu, uprim .* c)) .* (Kbar(nu, v) - Ybar(nu, uprim))) - ...
                besselj(nu+1, uprim.*c) .* bessely(nu+1, uprim) ./ (besselj(nu+1, uprim) .* bessely(nu+1, uprim .* c));
            result(IXclad) = resultCladding(IXclad);
        end;
        
        %% ariane
    case 'ari' % LP modes in TLS fibre, according to Ariane Stiebeiner 2009
        h = sqrt(n1.^2 .* k0.^2 - beta.^2);
        q = sqrt(beta.^2 - n3^2 .* k0.^2);
        p = sqrt(n2.^2 .* k0.^2 - beta.^2);
        s = sqrt(beta.^2 - n2.^2 .* k0.^2);
        if ~isempty(IXcore)
            % Core mode
            resultCore = besseli(nu, s.*b) .* besselk(nu, s.*a) .*...
                (-s.*besseli(nu+1, s.*b)./n_clad.^2./besseli(nu, s.*b) - ...
                q .* besselk(nu+1, q.* b) ./ besselk(nu, q.* b)).*...
                (s .* besselk(nu+1, s.*a)./n_clad.^2./besselk(nu, s.*a)-...
                h.*besselj(nu+1, h.*a)./n_core^2./besselj(nu, h.*a)) - ...
                besseli(nu, s.*a) .* besselk(nu, s.*b) .*...
                (-s.*besseli(nu+1, s.*a)./n_clad.^2./besseli(nu, s.*a) - ...
                h .* besselj(nu+1, h.*a)./n_core.^2./besselj(nu, h.*a)) .*...
                (s.*besselk(nu+1, s.*b)./n_clad.^2./besselk(nu, s.*b)-...
                q.*besselk(nu+1, q.*b)./besselk(nu, q.*b));
            result(IXcore) = resultCore(IXcore);
        end;
        if ~isempty(IXclad)
            % Cladding mode
            resultCladding = besselj(nu, p.*b) .* bessely(nu, p.*a) .*...
                (p.*besselj(nu+1, p.*b)./n_clad.^2./besselj(nu, p.*b) - ...
                q .* besselk(nu+1, q.* b) ./ besselk(nu, q.* b)).*...
                (p .* bessely(nu+1, p.*a)./n_clad.^2./bessely(nu, p.*a)-...
                h.*besselj(nu+1, h.*a)./n_core^2./besselj(nu, h.*a)) - ...
                besselj(nu, p.*a) .* bessely(nu, p.*b) .*...
                (p.*besselj(nu+1, p.*a)./n_clad.^2./besselj(nu, p.*a) - ...
                h .* besselj(nu+1, h.*a)./n_core.^2./besselj(nu, h.*a)) .*...
                (p.*bessely(nu+1, p.*b)./n_clad.^2./bessely(nu, p.*b)-...
                q.*besselk(nu+1, q.*b)./besselk(nu, q.*b));
            result(IXclad) = resultCladding(IXclad);
        end;
        
        %% tsao
    case 'tsaote'
        if ~isempty(IXcore)
            % Core mode
            W2 = sqrt(r1.^2 .* (beta.^2 - n2^2 * k^2));
            pNu = -2/pi*( besseli(nu,a*W2).*besselk(nu,W2)-besseli(nu,W2).*besselk(nu,a*W2) );
            qNuOverU2 = 2./(pi*W2).*( besseli(nu,a*W2).*besselkd(nu,W2)-besselid(nu,W2).*besselk(nu,a*W2) );
            rNuOverAU2 = 2./(pi*a*W2) .* ( besselid(nu,a*W2).*besselk(nu,W2)-besseli(nu,W2).*besselkd(nu,a*W2) );
            sNuOverAU2Squared = -2./(pi*a*W2.^2) .* ( besselid(nu,a*W2).*besselkd(nu,W2)-besselid(nu,W2).*besselkd(nu,a*W2) );
            
            resultCore = J .* (rNuOverAU2 + K .* pNu) - (K .* qNuOverU2 + sNuOverAU2Squared);
            result(IXcore) = resultCore(IXcore);
        end;
        if ~isempty(IXclad)
            % Cladding mode
            U2 = sqrt(r1.^2 .* (k^2 * n2^2 - beta.^2));
            % x1 = sqrt( k^2 .* n1^2 * U1.^4 .* U2.^4 ./ (nu^2 * beta.^2 .* V12.^4) );
            % x2 = sqrt( k^2 * n3^2 * a^4 * U2.^4 .* W3.^4 ./ (nu^2 .* beta.^2 .* V23.^4) );
            pNu = besselj(nu, a .* U2) .* bessely(nu, U2) - besselj(nu, U2) .* bessely(nu, a .* U2);
            qNu = besselj(nu, a .* U2) .* besselyd(nu, U2) - besseljd(nu, U2) .* bessely(nu, a .* U2);
            rNu = besseljd(nu, a .* U2) .* bessely(nu, U2) - besselj(nu, U2) .* besselyd(nu, a .* U2);
            sNu = besseljd(nu, a .* U2) .* besselyd(nu, U2) - besseljd(nu, U2) .* besselyd(nu, a .* U2);
            
            resultCladding = J .* (rNu ./ (a .* U2) + K .* pNu) - (K .* qNu ./ U2 + sNu ./ (a .* U2.^2));
            result(IXclad) = resultCladding(IXclad);
        end;
        
    case 'tsaotm'
        if ~isempty(IXcore)
            % Core mode
            W2 = sqrt(r1.^2 .* (beta.^2 - n2^2 * k^2));
            pNu = -2/pi*( besseli(nu,a*W2).*besselk(nu,W2)-besseli(nu,W2).*besselk(nu,a*W2) );
            qNuOverU2 = 2./(pi*W2).*( besseli(nu,a*W2).*besselkd(nu,W2)-besselid(nu,W2).*besselk(nu,a*W2) );
            rNuOverAU2 = 2./(pi*a*W2) .* ( besselid(nu,a*W2).*besselk(nu,W2)-besseli(nu,W2).*besselkd(nu,a*W2) );
            sNuOverAU2Squared = -2./(pi*a*W2.^2) .* ( besselid(nu,a*W2).*besselkd(nu,W2)-besselid(nu,W2).*besselkd(nu,a*W2) );
            
            resultCore = J .* (s23 .* rNuOverAU2 + K .* pNu) - s21 .* (K .* qNuOverU2 + s23 .* sNuOverAU2Squared);
            result(IXcore) = resultCore(IXcore);
        end;
        if ~isempty(IXclad)
            % Cladding mode
            U2 = sqrt(r1.^2 .* (k^2 * n2^2 - beta.^2));
            % x1 = sqrt( k^2 .* n1^2 * U1.^4 .* U2.^4 ./ (nu^2 * beta.^2 .* V12.^4) );
            % x2 = sqrt( k^2 * n3^2 * a^4 * U2.^4 .* W3.^4 ./ (nu^2 .* beta.^2 .* V23.^4) );
            pNu = besselj(nu, a .* U2) .* bessely(nu, U2) - besselj(nu, U2) .* bessely(nu, a .* U2);
            qNu = besselj(nu, a .* U2) .* besselyd(nu, U2) - besseljd(nu, U2) .* bessely(nu, a .* U2);
            rNu = besseljd(nu, a .* U2) .* bessely(nu, U2) - besselj(nu, U2) .* besselyd(nu, a .* U2);
            sNu = besseljd(nu, a .* U2) .* besselyd(nu, U2) - besseljd(nu, U2) .* besselyd(nu, a .* U2);
            
            resultCladding = J .* (s23 .* rNu ./ (a .* U2) + K .* pNu) - s21 .* (K .* qNu ./ U2 + s23 .* sNu ./ (a .* U2.^2));
            result(IXclad) = resultCladding(IXclad);
        end;
        
    case 'tsaohybrid' % hybrid Tsao mode
        if ~isempty(IXcore)
            % Core mode, see Tsao 1989 p. 557, right column, case 2: ne > n2
            W2 = sqrt(r1.^2 .* (beta.^2 - n2^2 * k^2));
            pNu = -2/pi*( besseli(nu,a*W2).*besselk(nu,W2)-besseli(nu,W2).*besselk(nu,a*W2) );
            qNuOverU2 = 2./(pi*W2).*( besseli(nu,a*W2).*besselkd(nu,W2)-besselid(nu,W2).*besselk(nu,a*W2) );
            rNuOverAU2 = 2./(pi*a*W2) .* ( besselid(nu,a*W2).*besselk(nu,W2)-besseli(nu,W2).*besselkd(nu,a*W2) );
            sNuOverAU2Squared = -2./(pi*a*W2.^2) .* ( besselid(nu,a*W2).*besselkd(nu,W2)-besselid(nu,W2).*besselkd(nu,a*W2) );
            
            U2 = 1j * W2;
            x1 = sqrt( k^2 * n1^2 * U1.^4 .* U2.^4 ./ (nu^2 * beta.^2 .* V12.^4) );
            x2 = sqrt( k^2 * n3^2 * a^4 * U2.^4 .* W3.^4 ./ (nu^2 .* beta.^2 .* V23.^4) );
            
            %qNu = besselj(nu, a .* U2) .* besselyd(nu, U2) - besseljd(nu, U2) .* bessely(nu, a .* U2);
            %rNu = besseljd(nu, a .* U2) .* bessely(nu, U2) - besselj(nu, U2) .* besselyd(nu, a .* U2);
            %sNu = besseljd(nu, a .* U2) .* besselyd(nu, U2) - besseljd(nu, U2) .* besselyd(nu, a .* U2);
            resultCore = pNu.^2 + 2*x1.*x2 * n2^2/(n1 *n3).*(2./(pi*a.*U2.^2)).^2 + ...
                x1.^2.*x2.^2 .*( J.*(rNuOverAU2+K.*pNu) - (K.*qNuOverU2+sNuOverAU2Squared) ) .*...
                ( J.*(s23*rNuOverAU2+K.*pNu ) - s21 .* (K.*qNuOverU2 + s23.*sNuOverAU2Squared) ) - ...
                x1.^2.*( J.*pNu- s21*qNuOverU2 ) .* (J .*pNu - qNuOverU2) -...
                x2.^2.*( K.*pNu + s23 * rNuOverAU2 ) .* (K .* pNu + s23 .* rNuOverAU2);
            result(IXcore) = resultCore(IXcore);
        end;
        if ~isempty(IXclad)
            % Cladding mode
            U2 = sqrt(r1.^2 .* (k^2 * n2^2 - beta.^2));
            x1 = sqrt( k^2 .* n1^2 * U1.^4 .* U2.^4 ./ (nu^2 * beta.^2 .* V12.^4) );
            x2 = sqrt( k^2 * n3^2 * a^4 * U2.^4 .* W3.^4 ./ (nu^2 .* beta.^2 .* V23.^4) );
            pNu = besselj(nu, a .* U2) .* bessely(nu, U2) - besselj(nu, U2) .* bessely(nu, a .* U2);
            qNu = besselj(nu, a .* U2) .* besselyd(nu, U2) - besseljd(nu, U2) .* bessely(nu, a .* U2);
            rNu = besseljd(nu, a .* U2) .* bessely(nu, U2) - besselj(nu, U2) .* besselyd(nu, a .* U2);
            sNu = besseljd(nu, a .* U2) .* besselyd(nu, U2) - besseljd(nu, U2) .* besselyd(nu, a .* U2);
            resultCladding = pNu.^2 + 2*x1.*x2 * n2^2/(n1 *n3).*(2./(pi*a.*U2.^2)).^2 + ...
                x1.^2.*x2.^2 .*( J.*(rNu./(a*U2)+K.*pNu) - (K.*qNu./U2+sNu./(a*U2.^2)) ) .*...
                ( J.*(s23*rNu./(a*U2)+K.*pNu ) - s21 .* (K.*qNu./U2 + s23.*sNu ./ (a*U2.^2)) ) - ...
                x1.^2.*( J.*pNu- s21*qNu./U2 ) .* (J .*pNu - qNu ./ U2) -...
                x2.^2.*( K.*pNu + s23 * rNu ./ (a*U2) ) .* (K .* pNu + s23 .* rNu ./ (a * U2));
            result(IXclad) = resultCladding(IXclad);
        end;
        
    %% erdogan
    % Erdogan eigen-value equation for three-layer vector mode
    % According to Turan Erdogan (1997, 2000) 10.1364/JOSAA.14.001760,
    % 10.1364/JOSAA.17.002113. LHS corresponds to zeta, RHS to zeta prime.
    case 'erdogan' 
        if ~isempty(IXcore)
            result(IXcore) = NaN;
        end
        
        if ~isempty(IXclad)
            % result(IXclad) = erdogan(d, neff(IXclad), lambda, fibreSpec, nu);
            a1 = d / 2 * fibreSpec.coreCladdingRatio;
            a2 = d / 2;
            r = a2;
            
            Zo = 377; % sqrt(mu0 / eps0); % 377 Ohm ?
            sigma1 = 1i * nu * neff / Zo;
            sigma2 = 1i * nu * neff * Zo;
            
            u1 = k0 .* sqrt(n1.^2 - neff.^2);
            u2 = k0 .* sqrt(n2.^2 - neff.^2);
            w3 = k0 .* sqrt(neff.^2 - n3.^2);
            J = besseljd(nu, u1.*a1) ./ u1 ./ besselj(nu, u1.*a1);
            K = besselkd(nu, w3.*a2) ./ w3 ./ besselk(nu, w3.*a2);
            
            p = besselj(nu, u2.*r) .* bessely(nu, u2.*a1) - ...
                besselj(nu,u2.*a1) .* bessely(nu,u2.*r);
            q = besselj(nu, u2.*r) .* besselyd(nu, u2.*a1) - ...
                besseljd(nu,u2.*a1) .* bessely(nu,u2.*r);
            rho = besseljd(nu, u2.*r) .* bessely(nu, u2.*a1) - ...
                besselj(nu,u2.*a1) .* besselyd(nu,u2.*r);
            s = besseljd(nu, u2.*r) .* besselyd(nu, u2.*a1) - ...
                besseljd(nu,u2.*a1) .* besselyd(nu,u2.*r);
            
            u21 = 1./u2.^2 - 1./u1.^2;
            u32 = 1./w3.^2 + 1./u2.^2;
            
            
            nominatorLhs = u2 .* (J.*K + ...
                sigma1 .* sigma2 .* u21 .* u32 ./ ...
                (n2.^2 .* a1 .* a2)) .* p - K .* q + ...
                J .* rho - 1./u2 .* s;
            denominatorLhs = -u2 .* (u32 ./ n2.^2 ./ a2 .* J - ...
                u21 ./ n1.^2 ./ a1 .* K) .* p + ...
                u32 ./ n1.^2 ./ a2 .* q + u21 ./ n1.^2 ./ a1 .* rho;
            lhs = 1./sigma2 .* nominatorLhs ./ denominatorLhs;
            
            nominatorRhs = u2 .* (u32 ./ a2 .* J - ...
                n3.^2 .* u21 ./ n2.^2 ./ a1 .* K) .* p - ...
                u32 ./ a2 .* q - u21 ./ a1 .* rho;
            denominatorRhs = u2 .* (n3.^2 ./ n2.^2 .* J .* K + ...
                sigma1 .* sigma2 .* u21 .* u32 ./ ...
                (n1.^2 .* a1 .* a2)) .* p - ...
                n3.^2 ./ n1.^2 .* K .* q + J .* rho - ...
                n2.^2 ./ n1.^2 ./ u2 .* s;
            rhs = sigma1 .* nominatorRhs ./ denominatorRhs;
            
            resultCladding = 1i * (lhs - rhs);
            result(IXclad) = resultCladding(IXclad);
        end;
        
    %% zhang
    case 'zhang'
        if ~isempty(IXcore)
            sigma0 = (nu * beta(IXcore) / k);
            s21 = (n2/n1).^2;
            s23 = (n2/n3).^2;
            V12 = k * r1 * sqrt(n1^2-n2^2);
            V23 = k * r2 * sqrt(n2^2-n3^2);
            u1 = k * sqrt(n1^2 - neff(IXcore).^2);
            w2 = k * sqrt(neff(IXcore).^2 - n2^2);
            w3 = k * sqrt(neff(IXcore).^2 - n3^2);
            U1 = u1 .* r1;
            W2 = w2 .* r1;
            W3 = w3 .* r2;
            Jb = besseljd(nu, U1) ./ U1 ./ besselj(nu, U1);
            Kb = besselkd(nu, W3) ./ W3 ./ besselk(nu, W3);
            X1 = n1 * U1.^2 .* W2.^2 ./ sigma0 ./ V12.^2;
            X2 = n2 * a.^2 .* W2.^2 .* W3.^2 ./ sigma0 ./ V23.^2;
            Pnu = besseli(nu, w2.*r2) .* besselk(nu, w2.*r1) -...
                besseli(nu, w2.*r1) .* besselk(nu, w2.*r2);
            Qnu = besseli(nu, w2.*r2) .* besselkd(nu, w2.*r1) -...
                besselid(nu, w2.*r1) .* besselk(nu, w2.*r2);
            Rnu = besselid(nu, w2.*r2) .* besselk(nu, w2.*r1) -...
                besseli(nu, w2.*r1) .* besselkd(nu, w2.*r2);
            Snu = besselid(nu, w2.*r2) .* besselkd(nu, w2.*r1) -...
                besselid(nu, w2.*r1) .* besselkd(nu, w2.*r2);
            result(IXcore) = ...
                Pnu.^2 - 2*(1./a./W2.^2).^2*(n2^2/n1/n3).*X1.*X2 + (X1.*X2).^2.*(Jb.*(Kb.*Pnu - ...
                Rnu./a./W2) + (1./W2).*(Kb.*Qnu-Snu./a./W2)) .* ...
                (Jb.*(Kb.*Pnu - s23*Rnu./a./W2) + (s21./W2).*(Kb.*Qnu-s23*Snu./a./W2)) -...
                X1.^2.*(Jb.*Pnu+Qnu./W2).*(Jb.*Pnu+s21*Qnu./W2)-...
                X2.^2.*(Kb.*Pnu-Rnu./a./W2).*(Kb.*Pnu-s23*Rnu./a./W2);
        end
        result(IXclad) = NaN * ones(size(IXclad));
    otherwise
        error('Invalid mode type: %s', upper(modeSpec.modetype));
end;




