clc; clear;

%% Q2
%(a)
eta = 2.15; beta = 0.96; alpha = 0.26;
vu = 0.60; rho = 0.88; sigmaEps = 0.014; delta = 0.07;

% Compute kStar
temp = alpha/(1/beta-1+delta);
kStar = ((vu/((1-temp*delta)*eta))^vu * temp)^(1/(1-alpha));

kL = 0.5*kStar; kH = 1.4*kStar; Nk = 37;
kgrid = logspace(log(kL - kL + 1)/log(10), log(kH - kL + 1)/log(10.0), Nk)';
kgrid = kgrid + ones(size(kgrid))*(kL-1); kgrid = kgrid';

s = sprintf( ['k(2): %.4f, ' ...
    'k(36): %.4f '], ...
    kgrid(2), kgrid(36));

disp(s)

%%
%(b)
Nd = 1000; Nm = 7;
dgrid = linspace(kL, kH, Nd);
mgrid = linspace(kL, kH, Nm);

% Compute the weight chi and assinged indexes
condMet = dgrid > kStar;
i = find(condMet, 1, "first") - 1;
chi = (dgrid(i+1) - kStar) / (dgrid(i+1) - dgrid(i));

s = sprintf('(i, i+1) = (%i, %i)', ...
    i, i+1);

disp(s)

s = sprintf('chi = %.4f', ...
    chi);

disp(s)

%%
%(c)

% Tauchen method
Nz = 5;

sigma = sqrt(sigmaEps^2/(1-rho^2));
xgrid = linspace(-2.575*sigma, 2.575*sigma, Nz);
zgrid = exp(xgrid);

s = sprintf('zgrid:');
disp(s)
disp(zgrid)

w = xgrid(2) - xgrid(1); Pi = zeros(Nz, Nz);
pd = makedist('Normal', 'mu', 0, 'sigma', sigmaEps);

for iz = 1:Nz
    Pi(iz, 1) = cdf(pd, xgrid(1) - rho*xgrid(iz) + w/2);
    Pi(iz, Nz) = 1 - cdf(pd, xgrid(Nz-1) - rho*xgrid(iz) + w/2);

    for j = 2:Nz-1
        Pi(iz, j) = cdf(pd, xgrid(j) - rho*xgrid(iz) + w/2)...
            - cdf(pd, xgrid(j-1) - rho*xgrid(iz) + w/2);
    end
end

s = sprintf('Transition Matrix:');
disp(s)
disp(Pi)

%%
%(d)

rng(1234567)
T = 10000;
zSim = zeros(1, T); izSim = zSim;
zSim(1) = zgrid(3); izSim(1) = 3;


cumPi = cumsum(Pi, 2);
efSim = rand(1, T);

% Do the simulation

for t = 1:T-1
    cSumVec  = cumPi(izSim(t), 1:Nz);
    condMet = efSim(t+1) <= cSumVec;
    izSim(t+1) = find(condMet, 1, "first");
end

zSim = zgrid(izSim);

s = sprintf('No. of (izt = 1): %i', sum(izSim==1));
disp(s)
s = sprintf('No. of (izt = 3): %i', sum(izSim==3));
disp(s)
s = sprintf('No. of (izt = 5): %i', sum(izSim==5));
disp(s)

tempmat = corrcoef(zSim(1:T-1),zSim(2:T));
rhoSim = tempmat(1,2);

s = sprintf('Persistence of Process: %.2f', rhoSim);
disp(s)

%%
%(e)
tol = 1e-5; 
betaM_0 = [0.017 0.806; 0.034 0.810; 0.051 0.798; 0.068 0.795; 0.086 0.793];
betaP_0 = [0.924 -0.415; 0.900 -0.412; 0.875 -0.409; 0.849 -0.407; 0.824 -0.404];

v = zeros(Nk, Nz, Nm); Tv = v; G = v;
vm = zeros(Nk, Nz); H = zeros(Nk, Nz, Nm);
iter = 0; distance = 10*tol;



while(distance>tol)
    for im = 1:Nm
        for iz = 1:Nz            
            % Compute vm by interpolation
            EmNext = exp(betaM_0(iz, 1) + betaM_0(iz, 2)*log(mgrid(im)));  
            for ik = 1:Nk
                vm(ik, iz) = interpolation(v(ik, iz, :), EmNext, mgrid);
            end
    
        end
        % Compute H() matrix 
        for iz = 1:Nz
            for ik = 1:Nk      
                H(ik, iz, im) = Pi(iz, :)*vm(ik, :)';
            end
        end
    
    end
    
    E0 = zeros(Nz, Nm); 
    for im = 1:Nm
    
        for iz = 1:Nz
            % Set up prices by the forecasting rule
            p = exp(betaP_0(iz, 1) + betaP_0(iz, 2)*log(mgrid(im)));
            w = eta/p;
            % Implement GSS method
            obj = @(kNext)  p*kNext-beta*interpolation(H(:, iz, im), kNext, kgrid);
            [kopt, E0(iz, im)] = goldenSearch(kL, kH, obj, false);
            for ik = 1:Nk
                n = ((zgrid(iz)*vu*kgrid(ik)^alpha)/w)^(1/(1-vu));
                pi = zgrid(iz)*(kgrid(ik)^alpha)*(n^vu) -...
                                w*n + (1-delta)*kgrid(ik);
                % Compute Tv and store the decision rule
                Tv(ik, iz, im) = pi*p - E0(iz, im);
                G(ik, iz, im) = kopt;
            end
    
        end
    
    
    end
    
    distance = max(max(max(abs(Tv-v))));
    v = Tv;
    iter = iter + 1;
    if(mod(iter, 20) == 0)
        s1 = sprintf( 'Iteration %i: ||Tv-v|| = %.5f', ...
             iter, distance);
        disp(s1)

        s2 = sprintf( '(Tvmin, Tvmax) =  (%.2f, %.2f)', ...
           min(min(min(Tv))), max(max(max(Tv))));
        disp(s2) 

        s3 = sprintf( '(kmin, kmax) =  (%.4f, %.4f)', ...
           min(min(min(G))), max(max(max(G))));
        disp(s3) 
    end
end


%%
%(f)

s1 = sprintf('EV0(k19|z1, m4)= %.4f', H(19, 1, 4));
disp(s1)
s2 = sprintf('EV0(k19|z3, m4)= %.4f', H(19, 3, 4));
disp(s2)
s3 = sprintf('EV0(k19|z5, m4)= %.4f', H(19, 5, 4));
disp(s3)
s4 = sprintf('EV0(k19|z3, m1)= %.4f', H(19, 3, 1));
disp(s4)
s5 = sprintf('EV0(k19|z3, m7)= %.4f', H(19, 3, 7));
disp(s5)
s6 = sprintf('EV0(k1|z3, m4)= %.4f', H(1, 3, 4));
disp(s6)
s7 = sprintf('EV0(k37|z3, m4)= %.4f', H(37, 3, 4));
disp(s7)

%%
%%Piecewise Linear Interpolation Function
function [v] = interpolation(v0, kNext, kgrid)

    if(kNext <= kgrid(1))
        v = v0(1);
    elseif(kNext >= kgrid(end))
        v = v0(end);
    else
        condMet = kNext <= kgrid;
        j = find(condMet, 1, "first") - 1;
        w = (kgrid(j+1) - kNext) / (kgrid(j+1) - kgrid(j));
        v = w*v0(j) + (1 - w)*v0(j+1);
    end

end

%%
%%Golden Search Function
function [x, fx] = goldenSearch(a, b, f, display)



    iter = 0;
    tol = 10e-7; distance = 10*tol;
    r = (3- sqrt(5)) / 2;
    
    c = (1-r)*a + r*b;
    d = r*a + (1-r)*b;
    
    
    while(distance > tol)
        if(f(c)>=f(d))
            a = c;
            c = d;
            d = r*a + (1-r)*b;
    
        else
            b = d;
            d = c;
            c = (1-r)*a + r*b;
            
        end
        distance = b - a;

        if(display==true)
            s = sprintf( 'Iteration %i: [a, b, c, d] = [%.4f, %.4f, %.4f, %.4f], |b-a| = %.4f', ...
            iter, a, b, c, d, distance);
            disp(s)
        end

        iter = iter + 1;


    end
    
    x = (a+b) / 2;
    fx = f(x);
    if(display==true)
        s = sprintf( '(x^*, f(x^*)) = (%.4f, %.4f)', ...
        x, fx);
        disp(s)
    end
end

