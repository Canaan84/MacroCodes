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

for t = 1:T-1
    cSumVec  = cumPi(izSim(t), 1:Nz);
    condMet = efSim(t+1) <= cSumVec;
    izSim(t+1) = find(condMet, 1, "first");
end

zSim = zgrid(izSim);

%%
%(b)


tol1 = 1e-8; 
betaM_0 = [0.017 0.806; 0.034 0.810; 0.051 0.798; 0.068 0.795; 0.086 0.793];
betaP_0 = [0.924 -0.415; 0.900 -0.412; 0.875 -0.409; 0.849 -0.407; 0.824 -0.404];

v = zeros(Nk, Nz, Nm); Tv = v; G = v;
vm = zeros(Nk, Nz); H = zeros(Nk, Nz*Nm);
iter = 0; distance = 10*tol1;



while(distance>tol1)
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

% Compute the weight chi and assinged indexes
condMet = dgrid > kStar;
i = find(condMet, 1, "first") - 1;
chi = (dgrid(i+1) - kStar) / (dgrid(i+1) - dgrid(i));


s = sprintf('(i, chi) = (%i, %.4f)', ...
    i, chi);

disp(s)

mu1 = zeros(1, Nd); 
mu1(i) = chi; mu1(i+1) = 1 - chi;

%%
%(d)
tol2 = 1e-10;
ktopt = zeros(1, T); mut = zeros(T, Nd);
mt = ktopt; pt = ktopt; wt = ktopt;
yt = ktopt; nt = ktopt; it = ktopt; ct = ktopt;
mut(1, :) = mu1;

for t = 1: T
    cL = 0.01; cH = 2.0; 
    mt(t) = mut(t, :)*dgrid';
    condMet = mgrid > mt(t);
    im = find(condMet, 1, "first") - 1;
    wL = (mgrid(im+1) - mt(t)) / (mgrid(im+1) - mgrid(im));
    izt = izSim(t);
    
    Ht = zeros(1, Nk);
    
    for ik = 1:Nk
        Ht(ik) = wL*H(ik, (izt-1)*Nm+im)+ (1-wL)*H(ik, (izt-1)*Nm+im+1);
    end
    
    
    
    iter = 0; distance = 10*tol2;
    while(distance> tol2)
        cM = (cL+cH)/2;
        fL = fC(cL, izt, Ht, mut(t, :), kgrid, dgrid, zgrid, kL, kH, Nd, ...
    beta, eta, alpha, nu, delta); 
        fM = fC(cM, izt, Ht, mut(t, :), kgrid, dgrid, zgrid, kL, kH, Nd, ...
    beta, eta, alpha, nu, delta);
    
        if fL*fM < 0
           cH = cM;
        else
           cL = cM;
        end
    
        distance = cH - cL;
        iter = iter + 1;
    end
    [fM, ktopt(t), yt(t), it(t), nt(t)] = fC(cM, izt, Ht, mut(t, :), kgrid, dgrid, zgrid, kL, kH, Nd, ...
    beta, eta, alpha, nu, delta);
    ct(t) = cM;
    pt(t) = 1/cM; wt(t) = eta/pt(t);
    
    
    condMet = dgrid > ktopt(t);
    i = find(condMet, 1, "first") - 1;
    chi = (dgrid(i+1) - ktopt(t)) / (dgrid(i+1) - dgrid(i));
    
    
    mut(t+1, i) = chi; mut(t+1, i+1) = 1 - chi;

    if(mod(t, 200) == 0)
        s1 = sprintf('Date %i: (iz, m) = (%i, %.4f)', ...
        t, izt, mt(t));

        disp(s1)
        s2 = sprintf('market clearing p = %.4f', ...
        pt(t));

        disp(s2)        

        s3 = sprintf('(y, i, c, n) = (%.3f, %.3f, %.3f, %.3f)', ...
        yt(t), it(t), ct(t), nt(t));

        disp(s3)    

        s4 = sprintf('X^D = %.4f', ...
        fM);

        disp(s4)   
    end

end

%%
%(e)
betaM_1 = zeros(Nz, 2); betaP_1 = zeros(Nz, 2); 
RAdj_M = zeros(Nz, 1); RAdj_P = zeros(Nz, 1);
for iz = 1:Nz
    tZ = find(izSim(1:end-1) == iz);
    X = log(mt(tZ)');
    Y = log(mt(tZ+1)');
    mdl = fitlm(X,Y);
    betaM_1(iz, :) = mdl.Coefficients.Estimate';
    RAdj_M(iz) = mdl.Rsquared.Adjusted;

    Y = log(pt(tZ+1)');
    mdl = fitlm(X,Y);
    betaP_1(iz, :) = mdl.Coefficients.Estimate';
    RAdj_P(iz) = mdl.Rsquared.Adjusted;
end



s1 = sprintf('betaM_1 =');
disp(s1)
disp(betaM_1)

s2 = sprintf('betaP_1 =');
disp(s2)
disp(betaP_1)

s3 = sprintf('RAdj_M =');
disp(s3)
disp(RAdj_M)

s4 = sprintf('RAdj_P =');
disp(s4)
disp(RAdj_P)

s1 = sprintf('|| betaM_1 - betaM_0 ||=');
disp(s1)
disp(abs(betaM_1 - betaM_0))

s2 = sprintf('|| betaP_1 - betaP_0 ||=');
disp(s2)  
disp(abs(betaP_1 - betaP_0))

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



%%
%%fC function
function [XD, ktopt, y, i, n] = fC(C, izt, Ht, mu, kgrid, dgrid, zgrid, kL, kH, Nd, ...
    beta, eta, alpha, nu, delta)
    
    obj = @(kNext) 1/C*kNext - beta*interpolation(Ht, kNext, kgrid);
    ktopt = goldenSearch(kL, kH, obj, false);
    
    yt = zeros(1, Nd);
    for jk = 1:Nd
        if(mu(jk)>0)
            w = eta*C;
            n = ((zgrid(izt)*nu*dgrid(jk)^alpha)/w)^(1/(1-nu));
            yt(jk) = zgrid(izt)*(dgrid(jk)^alpha)*(n^nu);
        end
    
    end
    y = mu*yt';
    i = mu*(ktopt - (1-delta)*dgrid)';
    Cs = y - i;
    XD = C - Cs;

end