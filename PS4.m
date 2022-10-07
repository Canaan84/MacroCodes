clc; clear;

%%

%%Q3
%(a)

eta = 2.15; beta = 0.96; alpha = 0.26; nu = 0.60;
rho = 0.88; sigmaEps = 0.014; delta = 0.07;

% Compute kStar
temp = alpha/(1/beta-1+delta);
kStar = ((nu/((1-temp*delta)*eta))^nu * temp)^(1/(1-alpha));

kL = 0.8*kStar; kH = 1.25*kStar; Nk = 125;
kgrid = logspace(log(kL - kL + 1)/log(10), log(kH - kL + 1)/log(10.0), Nk)';
kgrid = kgrid + ones(size(kgrid))*(kL-1); kgrid = kgrid';

Nd = 1000; Nm = 7;
dgrid = linspace(kL, kH, Nd);
mgrid = linspace(kL, kH, Nm);

Nz = 5;

sigma = sqrt(sigmaEps^2/(1-rho^2));
xgrid = linspace(-2.575*sigma, 2.575*sigma, Nz);
zgrid = exp(xgrid);

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


rng(1234567)
T = 10000;
zSim = zeros(1, T); izSim = zSim;
zSim(1) = zgrid(3); izSim(1) = 3;


cumPi = cumsum(Pi, 2);
efSim = rand(1, T);

%%
%(b)


tol = 1e-8; 
betaM_0 = [0.017 0.806; 0.034 0.810; 0.051 0.798; 0.068 0.795; 0.086 0.793];
betaP_0 = [0.924 -0.415; 0.900 -0.412; 0.875 -0.409; 0.849 -0.407; 0.824 -0.404];

v = zeros(Nk, Nz, Nm); Tv = v; G = v;
vm = zeros(Nk, Nz); H = zeros(Nk, Nz*Nm);
iter = 0; distance = 10*tol;



while(distance>tol)
    for im = 1:Nm
        for iz = 1:Nz
            
            
            EmNext = exp(betaM_0(iz, 1) + betaM_0(iz, 2)*log(mgrid(im)));
            
            condMet = mgrid > EmNext;
            i = find(condMet, 1, "first") - 1;
            wm = (mgrid(i+1) - EmNext) / (mgrid(i+1) - mgrid(i));
            
            for ik = 1:Nk
                for jz = 1:Nz
                    vm(ik, jz) = wm*v(ik, jz, i) + (1-wm)*v(ik, jz, i+1);
                end
                izm = Nm*(iz-1) + im;
        
                H(ik, izm) = Pi(iz, :)*vm(ik, :)';
    
            end
        end
    
    end
    
    E0 = zeros(Nz, Nm); 
    for im = 1:Nm
    
        for iz = 1:Nz
            izm = Nm*(iz-1) + im;
            % Set up prices by the forecasting rule
            p = exp(betaP_0(iz, 1) + betaP_0(iz, 2)*log(mgrid(im)));
            w = eta/p;
            
            % Implement GSS method
            obj = @(kNext)  p*kNext-beta*interpolation(H(:, izm), kNext, kgrid);
            [kopt, E0(iz, im)] = goldenSearch(kL, kH, obj, false);
            for ik = 1:Nk
                n = ((zgrid(iz)*nu*kgrid(ik)^alpha)/w)^(1/(1-nu));
                pi = zgrid(iz)*(kgrid(ik)^alpha)*(n^nu) -...
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


s = sprintf('H(k63|(z3, m4))= %.4f', H(63, Nm*(3-1)+4));
disp(s)


%%
%(c)



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
    tol = 10e-10; distance = 10*tol;
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