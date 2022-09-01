
%%Q1

%(a)
sigmaEps = 0.00712; rho = 0.956; Nz = 9;
medInd = floor((Nz+1)/2);

sigma = sqrt(sigmaEps^2/(1-rho^2));

x = linspace(-2.575*sigma, 2.575*sigma, Nz);
z = exp(x);

s = sprintf( 'Median of x:%.4f, Median of z:%.4f', ...
    median(x), median(z));
disp(s)


%(b)

w = x(2) - x(1); Pi = zeros(Nz, Nz);
pd = makedist('Normal', 'mu', 0, 'sigma', 0.00712);

for iz = 1:Nz
    Pi(iz, 1) = cdf(pd, x(1) - rho*x(iz) + w/2);
    Pi(iz, Nz) = 1 - cdf(pd, x(Nz-1) - rho*x(iz) + w/2);

    for j = 2:Nz-1
        Pi(iz, j) = cdf(pd, x(j) - rho*x(iz) + w/2)...
            - cdf(pd, x(j-1) - rho*x(iz) + w/2);
    end
end

s = sprintf('Transition Matrix:');
disp(s)
disp(Pi)
s = sprintf('Transition Probability: %.4f', Pi(medInd, medInd));
disp(s)


%(c) 

T = 1000; zSim = zeros(1, T); izSim = zSim;

cumPi = cumsum(Pi, 2);
rng(1234567);
efSim = rand(1, T);

izSim(1) = medInd;

for t = 1:T-1
    cSumVec  = cumPi(izSim(t), 1:Nz);
    condMet = efSim(t+1) <= cSumVec;
    izSim(t+1) = find(condMet, 1, "first");
end

zSim = z(izSim);
s = sprintf('Persistence of Process: %.4f', corrcoef(zSim));
disp(s)

%% 

%%Q2


%(a)


alpha = 0.261; beta = 0.99; delta = 0.0176; Nk = 275;

kStar = (alpha / (1/beta - 1 + delta))^(1/(1-alpha));

k = linspace(0.5*kStar, 1.5*kStar, Nk);
s = sprintf( 'Median of k: %.4f', ...
    median(k));
disp(s)


%(b)

tol = 0.001; IKHat = zeros(Nk, Nz); v = zeros(Nk, Nz); Tv = v;
iter = 1; distance = 10*tol; G = v;


for ik = 1:Nk
    for iz = 1:Nz
        IKHat(ik, iz) = sum(z(iz)*k(ik)^alpha + (1-delta)*k(ik) >= k);
    end
end




while(distance > tol)
    for iz = 1:Nz
        for ik = 1:Nk
            X = zeros(1, Nk);
            for jk = 1:IKHat(ik, iz)            
                R = log(z(iz)*k(ik)^alpha + (1-delta)*k(ik) - k(jk));
                Ev = Pi(iz, :)*v(jk, :).';

                X(jk) = R + beta*Ev;
            end
            [Tv(ik, iz), G(ik, iz)] = max(X);

        end


    end
    distance = max(max(abs(Tv-v)));
    v = Tv;
    iter = iter + 1;
    s = sprintf( 'Iteration %i: ||Tv-v|| = %.4f', ...
    iter, distance);
    disp(s)

    
end

kNext = k(G);

s = sprintf( ['Max of decision rule: %.4f, ' ...
    'Min of decision rule: %.4f'], ...
    max(max(kNext)), min(min(kNext)));

disp(s)

figure()

[kgrid, zgrid] = ndgrid(k, z);

sp(1) = subplot(2,1,1); 
surf(kgrid, zgrid, v)
title('Value Function')
xlabel('k') 
ylabel('z')
zlabel('v')

sp(2) = subplot(2,1,2); 
surf(kgrid, zgrid, kNext)
title('Decision Rule')
xlabel('k') 
ylabel('z')
zlabel('g(k, z)')


%(c)

kSim = zeros(1, T); ikSim = kSim; cSim = kSim; iSim = kSim; ySim = kSim;
kSim(1) = k(10); ikSim(1) = 10;

for t = 1:T-1
    ikSim(t+1) = G(ikSim(t), izSim(t));
    kSim(t+1) = k(ikSim(t+1));
    ySim(t) = zSim(t)*kSim(t)^alpha;
    cSim(t) = ySim(t) + (1-delta)*kSim(t) - kSim(t+1);
    iSim(t) = kSim(t+1) - (1-delta)*kSim(t);
end

figure()

sp(1) = subplot(2,2,1); 
plot(900:T-1, zSim(900:T-1))
legend({'TFP'}, 'FontSize',6)

sp(2) = subplot(2,2,2); 
plot(900:T-1, ySim(900:T-1))
legend({'Output'}, 'FontSize',6)

sp(3) = subplot(2,2,3); 
plot(900:T-1, cSim(900:T-1))
legend({'Consumption'}, 'FontSize',6)

sp(4) = subplot(2,2,4); 
plot(900:T-1, iSim(900:T-1))
legend({'Investment'}, 'FontSize',6)

%% 
%%Q3

theta = 0.5; delta = 0.08; beta = 1/1.046; tol = 1e-9;
diff = tol*10;

kAStar = (theta/(1/beta-1+delta))^(1/(1-theta));
kL = 0.01; kH = 100;

f = @(x) theta*(x)^(theta-1) - 1/beta + 1 - delta;

while(diff>tol)
    kMed = (kL+kH)/2;
    if f(kL)*f(kMed) < 0
       kH = kMed;
    else
       kL = kMed;
    end

    diff = kH - kL;
end

s = sprintf( ['Analytical sol. of k: %.4f, ' ...
    'Numerical sol. of k: %.4f, ', 'Error: %.4f '], ...
    kAStar, (kH+kL)/2, diff/(1e-13));

disp(s)

%% 

%%Q4

%(b)

Nz = 2; Pi = [0.95, 0.05; 0.05, 0.95]; z = [0.96, 1.04];
alpha = 0.20; nu = 0.60; p = 1; w = 1;
delta = 0.10; beta = 0.96; phi = 0.50;
Nk = 200; 

k = linspace(0.01, 0.60, Nk);
tol = 1e-5; v = zeros(Nk, Nz); Tv = v;
iter = 1; distance = 10*tol; G = v;

R = zeros(Nk, Nz);

for ik = 1:Nk
    for iz = 1:Nz
        R(ik, iz) = (1-nu)*(z(iz)*((nu/w)^nu)*(k(ik)^alpha))^(1/(1-nu));
    end
end


psi = @(k, kp) phi/2*(((kp-(1-delta)*k)/k)^2)*k;

while(distance > tol)
    
    for iz = 1:Nz
        for ik = 1:Nk

            X = zeros(1, Nk);
            for jk = 1:Nk          
                capitalCost = p*(k(jk) - (1-delta)*k(ik)) +...
                    psi(k(ik), k(jk));
                Ev = Pi(iz, :)*v(jk, :).';

                X(jk) = R(ik, iz) - capitalCost + beta*Ev;

            
            end
            [Tv(ik, iz), G(ik, iz)] = max(X);

        end


    end
    distance = max(max(abs(Tv-v)));
    v = Tv;
    iter = iter + 1;
    s = sprintf( 'Iteration %i: ||Tv-v|| = %.6f', ...
    iter, distance);
    disp(s)
end

kNext = k(G);

s = sprintf( ['Max of decision rule: %.4f, ' ...
    'Min of decision rule: %.4f'], ...
    max(max(kNext)), min(min(kNext)));

disp(s)

%(c)

Dv = zeros(Nk-1, Nz); q = zeros(Nk-1, Nz); Q = zeros(Nk, Nz);

for l = 1:Nk-1
    Dv(l,:) = (v(l+1, :) - v(l, :))/(k(l+1) - k(l));
end

for iz = 1:Nz
    for lk = 1:Nk
        if lk ~= Nk
            q(lk, iz) = beta*Pi(iz, :)*Dv(lk, :).';
        end
        Q(lk, iz) = beta*(Pi(iz, :)*v(lk, :).')/k(lk);
    end
end



figure()

[kgrid, zgrid] = ndgrid(k, z);

sp(1) = subplot(2,1,1); 
surf(kgrid, zgrid, Q)
title('Average Q')
xlabel('k') 
ylabel('z')
zlabel('Q')


[kgrid, zgrid] = ndgrid(k(1:Nk-1), z);


sp(2) = subplot(2,1,2); 
surf(kgrid, zgrid, q)
title('Marginal q')
xlabel('k') 
ylabel('z')
zlabel('q')


%(d)

T = 8000; kSim = zeros(1, T); iSim = kSim; 
ikSim = kSim; izSim = kSim;

ikSim(1) = 100;

cumPi = cumsum(Pi, 2);
rng(123456);
efSim = rand(1, T);

izSim(1) = 1; i_k_ratio = zeros(1, T-1); 
QSim = zeros(1, T-1);

for t = 1:T-1
    cSumVec  = cumPi(izSim(t), 1:Nz);
    condMet = efSim(t+1) <= cSumVec;
    izSim(t+1) = find(condMet, 1, "first");
    ikSim(t+1) = G(ikSim(t), izSim(t));
    kSim(t+1) = k(ikSim(t+1));
    iSim(t) = kSim(t+1) - (1-delta)*kSim(t);
    i_k_ratio(t) = iSim(t)/kSim(t);

    QSim(t) = beta*Pi(izSim(t), :)*v(ikSim(t+1), :)'/kSim(t+1);
end

fitlm(QSim', i_k_ratio')


%(e) 

profitRate = zeros(1, T-1);
for t = 1:T-1
    profitRate(t) = R(ikSim(t), izSim(t))/k(ikSim(t));
end

X = [QSim; profitRate];
fitlm(X', i_k_ratio')



%(f)


Pi = [0.57, 0.43; 0.57, 0.43];
v = zeros(Nk, Nz); Tv = v;
iter = 1; distance = 10*tol; G = v;

R = zeros(Nk, Nz);

for ik = 1:Nk
    for iz = 1:Nz
        R(ik, iz) = (1-nu)*(z(iz)*((nu/w)^nu)*(k(ik)^alpha))^(1/(1-nu));
    end
end




while(distance > tol)
    
    for iz = 1:Nz
        for ik = 1:Nk

            X = zeros(1, Nk);
            for jk = 1:Nk          
                capitalCost = p*(k(jk) - (1-delta)*k(ik)) +...
                    psi(k(ik), k(jk));
                Ev = Pi(iz, :)*v(jk, :).';

                X(jk) = R(ik, iz) - capitalCost + beta*Ev;

            
            end
            [Tv(ik, iz), G(ik, iz)] = max(X);

        end


    end
    distance = max(max(abs(Tv-v)));
    v = Tv;
    iter = iter + 1;
    s = sprintf( 'Iteration %i: ||Tv-v|| = %.6f', ...
    iter, distance);
    disp(s)
end

kNext = k(G);


for t = 1:T-1
    cSumVec  = cumPi(izSim(t), 1:Nz);
    condMet = efSim(t+1) <= cSumVec;
    izSim(t+1) = find(condMet, 1, "first");
    ikSim(t+1) = G(ikSim(t), izSim(t));
    kSim(t+1) = k(ikSim(t+1));
    iSim(t) = kSim(t+1) - (1-delta)*kSim(t);
    i_k_ratio(t) = iSim(t)/kSim(t);

    QSim(t) = beta*Pi(izSim(t), :)*v(ikSim(t+1), :)'/kSim(t+1);
end


profitRate = zeros(1, T-1);


for t = 1:T-1
    profitRate(t) = R(ikSim(t), izSim(t))/k(ikSim(t));
end



X = [QSim; profitRate];
fitlm(X', i_k_ratio')