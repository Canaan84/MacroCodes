%%Q1
%(d)

%(i)

delta = 0.08; kStar = 1.6;
alpha = [0.10, 0.25, 0.45, 0.70, 0.96, 1.0];
theta = zeros(1, length(alpha));
den = 1+ sum(cumprod(ones(1, length(alpha)) - alpha));
theta(1) = 1/den;

% Solve for the stationary distribution by formula
for j = 2:length(alpha)
    theta(j) = theta(j-1)*(1-alpha(j-1));

end

%(ii)
% Compute the number of adjustors and aggregate I
NumOfAdj = theta*alpha';

IVector = kStar* (ones(1, length(alpha)) - (1-delta).^(1:length(alpha)));
aggInv = sum(theta.*alpha.*IVector);
s = sprintf( 'No. of Adjustors: %.4f, Aggregate Investment: %.4f', ...
      NumOfAdj, aggInv);

disp(s)

%(iii)

kVector = kStar*(1-delta).^(0:length(alpha)-1);
plot(kVector, theta)
title('Stationary Distribution of Plants over Capital')
xlabel('k') 
ylabel('mu(k)')

%%
%%Q3

f = @(x) 20 + x.^2 - 8*x -log(x);

a = 0; b = 10;
% See the golden search section function below
[xStar, fxStar] = goldenSearch(a, b, f, true);

%%
%%Q4
%(a)

display = true;

q = 1;
% See the function Step 1 below
[v, G] = Step1(q, display);


%%
%(b)
% See the function Step 2 below
[mu, Gd] = Step2(G, display);

%%
%(c)
% See the function Step 3 below
B = Step3(mu, Gd, display);

%%
%(d)

display = false;
beta = 0.99322;
iter = 0;
qL = beta + 1e-4; qH = 1.02; tol3 = 1e-5;
distance3 = tol3*10; 
% Compute the excess demands for the initial bounds
[vL, GL] = Step1(qL, display);
[muL, GdL] = Step2(GL, display);
BL = Step3(muL, GdL, display);


[vH, GH] = Step1(qH, display);
[muH, GdH] = Step2(GH, display);
BH = Step3(muH, GdH, display);

BM = tol3*10; 

while(distance3>tol3 && abs(BM)>tol3)
    % Compute the excess demand for the median q
    % Perform bisection search
    qM = (qL+qH)/2;
    [vM, GM] = Step1(qM, display);
    [muM, GdM] = Step2(GM, display);
    BM = Step3(muM, GdM, display);

    if BL*BM < 0
       qH = qM;
    else
       qL = qM;
    end

    distance3 = qH - qL;

    iter = iter + 1;
    s = sprintf( 'Iteration %i: q = %.6f, B = %.6f', ...
          iter, qM, BM);
    disp(s)
end

q = (qH+qL)/2;
s = sprintf( 'equilibrium q = %.6f, annualized r = %.6f', ...
      q, 1/q-1);
disp(s)


%%
%(e)
Nbd = 1000; Nz = 2;
C = zeros(Nbd, Nz);
bL = -4.0; bH = 4.0; Nb1 = 100;
b1grid = logspace(log(bL - bL + 1)/log(10), log(bH - bL + 1)/log(10.0), Nb1)';
b1grid = b1grid + ones(size(b1grid))*(bL-1);
wgrid = [0.1, 1.0]; 

Nbd = 1000; mu = zeros(Nbd, Nz); 
bdgrid = linspace(b1grid(1), 4, Nbd);
% Compute the distribution of consumption
for iz = 1:Nz
    
    for ibd = 1:Nbd
        C(ibd, iz) = bdgrid(ibd) + wgrid(iz) - GdM(ibd, iz);

    end

end
% Compute the mean and std of consumption
CMean = sum(sum(C.*muM));
CVar = sum(sum(((C-CMean).^2).*muM));

s = sprintf( 'Mean(C) = %.4f, Std(C) = %.4f', ...
      CMean, sqrt(CVar));
disp(s)


%%
%%Golden Search Function
function [x, fx] = goldenSearch(a, b, f, display)



    iter = 0;
    tol = 10e-5; distance = 10*tol;
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




%%Piecewise Linear Interpolation Function
function [v] = interpolation(v0, bNext, bgrid)

    if(bNext <= bgrid(1))
        v = v0(1);
    elseif(bNext >= bgrid(end))
        v = v0(end);
    else
        condMet = bNext <= bgrid;
        j = find(condMet, 1, "first") - 1;
        w = (bgrid(j+1) - bNext) / (bgrid(j+1) - bgrid(j));
        v = w*v0(j) + (1 - w)*v0(j+1);
    end

end

%%Step 1
function [v, G] = Step1(q, display)


    
    Nz = 2;
    wgrid = [0.1, 1.0]; Pi = [0.5, 0.5; 0.075, 0.925];
    beta = 0.99322; sigma = 1.5;
    u = @(c) (c^(1-sigma) - 1) / (1 - sigma);
    
    bL = -4.0; bH = 4.0; Nb1 = 100;
    b1grid = logspace(log(bL - bL + 1)/log(10), log(bH - bL + 1)/log(10.0), Nb1)';
    b1grid = b1grid + ones(size(b1grid))*(bL-1); b1grid = b1grid';

    tol1 = 1e-4; distance1 = 10*tol1; eps = 1e-5; iter = 0;
    
    v = zeros(Nb1, Nz); Tv = v; G = v; v0 = v;
    
    
    while(distance1 > tol1)
        
    
        for iz = 1:Nz
            
            v0(:, iz) = Pi(iz, :)*(v(:, :)');
            for ib = 1:Nb1
                obj = @(bNext) - u(b1grid(ib) + wgrid(iz) - q*bNext) ...
                       - beta*interpolation(v0(:, iz), bNext, b1grid);
                bLowerBar = b1grid(1);
                bUpperBar = (b1grid(ib) + wgrid(iz))/q - eps;
                [G(ib, iz), Tv(ib, iz)] = goldenSearch(bLowerBar, bUpperBar,...
                    obj, false);
                Tv(ib, iz) = -Tv(ib, iz);
        
            end
    
    
        end
    
        distance1 = max(max(abs(Tv-v)));
        v = Tv;
        iter = iter + 1;
        if(mod(iter, 100) == 0 && display == true)
            s1 = sprintf( 'Iteration %i: ||Tv-v|| = %.5f', ...
                 iter, distance1);
            disp(s1)
    
            s2 = sprintf( '(bmin, bmax) =  (%.5f, %.5f)', ...
               min(min(G)), max(max(G)));
            disp(s2) 
        end
    
    
    end
    if(display == true)
        s3 = sprintf( 'q = %.4f', q);
        disp(s3)     
    end

end

%%Step2

function [mu, Gd] = Step2(G, display)
    
    Nz = 2;
    Pi = [0.5, 0.5; 0.075, 0.925];
    
    
    
    bL = -4.0; bH = 4.0; Nb1 = 100;
    b1grid = logspace(log(bL - bL + 1)/log(10), log(bH - bL + 1)/log(10.0), Nb1)';
    b1grid = b1grid + ones(size(b1grid))*(bL-1); b1grid = b1grid';


    Nbd = 1000; mu = zeros(Nbd, Nz); 
    bdgrid = linspace(b1grid(1), 4, Nbd);
    mu(1, 1) = 1; Gd = zeros(Nbd, Nz);
    
    tol2 = 1e-5; distance2 = 10*tol2; iter = 0;
    
    while(distance2 > tol2)
    
        Tmu = zeros(Nbd, Nz);
        for iz = 1:Nz
            for ib = 1:Nbd
                
                if(mu(ib, iz)>0)
                    if(Gd(ib, iz) ~= 0)
                        bNext = Gd(ib, iz);
                    else
                        Gd(ib, iz) = interpolation(G(:, iz), bdgrid(ib), b1grid);
                        bNext = Gd(ib, iz);
                    end
                    condMet = bNext < bdgrid;
                    j = find(condMet, 1, "first") - 1;
                    w = (bdgrid(j+1) - bNext) / (bdgrid(j+1) - bdgrid(j));
                end
                for lz = 1:2
                    Tmu(j, lz) = Tmu(j, lz) + w*Pi(iz, lz)*mu(ib, iz);
                    Tmu(j+1, lz) = Tmu(j+1, lz) + (1-w)*Pi(iz, lz)*mu(ib, iz);
                end
        
            end
        
        
        end
        
        distance2 = max(max(abs(Tmu-mu)));
        
        
        mu = Tmu;
        iter = iter + 1;
        if(mod(iter, 50) == 0 && display == true)
            s = sprintf( 'Iteration %i: ||Tmu-mu|| = %.5f', ...
                 iter, distance2);
            disp(s)
        end
    
    end



end

%%Step3

function [B] = Step3(mu, Gd, display)
    
    Nz = 2; Nbd = 1000;

    B = 0;
    for iz = 1:Nz
        for ib = 1:Nbd
            B = B + mu(ib, iz)*Gd(ib, iz);
    
        end
    end
    if(display)
        s = sprintf( 'Aggregate Bond Demand at q = 1: %.4f', ...
             B);
        disp(s)
    end
end
