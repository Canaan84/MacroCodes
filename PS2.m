%%Q3

f = @(x) 20 + x^2 - 8*x -log(x);

tol = 10e-5; distance = 10*tol;
iter = 1;
r = (3- sqrt(5)) / 2;

a = 0; b = 10;
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
    
    s = sprintf( 'Iteration %i: [a, b, c, d] = [%.4f, %.4f, %.4f, %.4f], |b-a| = %.4f', ...
    iter, a, b, c, d, distance);
    disp(s)

    iter = iter + 1;
end

xStar = (a+b) / 2;

s = sprintf( '(x^*, f(x^*)) = (%.4f, %.4f)', ...
    xStar, f(xStar));
disp(s)

%%
%%Q4

Nz = 2;
wgrid = [0.1, 1.0]; Pi = [0.5, 0.5; 0.075, 0.925];
beta = 0.99322; sigma = 1.5;
u = @(c) (c^(1-sigma) - 1) / (1 - sigma);

bL = -4.0; bH = 4.0; Nb1 = 100;
b1grid = logspace(log(bL - bL + 1)/log(10), log(bH - bL + 1)/log(10.0), Nb1)';
b1grid = b1grid + ones(size(b1grid))*(bL-1);

%(a)

tol = 1e-4; distance = 10*tol;

v = zeros(Nb1, Nz); Tv = v; G = v; v0 = v;


q = 1;

while(distance > tol)

    for iz = 1:Nz
        v0(:, iz) = Pi(iz, :)*(v(:, :)');
    
        for ib = 1:Nb1
            obj = @(bNext) u(b1grid(ib) + wgrid(iz) - q*bNext) ...
                   + beta*interpolation(v0, bNext, iz, b1grid);
            bLowerBar = b1grid(1);
            bUpperBar = (b1grid(ib) + wgrid(iz))/q;
            [G(ib, iz), Tv(ib, iz)] = goldenSearch(bLowerBar, bUpperBar, obj);
    
        end


    end

    distance = max(max(abs(Tv-v)));
    v = Tv;
    iter = iter + 1;
    s = sprintf( 'Iteration %i: ||Tv-v|| = %.6f', ...
    iter, distance);
    disp(s)


end





%%
%%Golden Search Function
function [x, fx] = goldenSearch(a, b, f)


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
    end
    
    x = (a+b) / 2;
    fx = f(x);
end

%%Piecewise Linear Interpolation
function [v0] = interpolation(v, bNext, iz, bgrid)
    if(bNext <= bgrid(1))
        v0 = v(1, iz);
    elseif(bNext >= bgrid(end))
        v0 = v(end, iz);
    else
        condMet = bNext <= bgrid;
        j = find(condMet, 1, "first") - 1;
        w = (bgrid(j+1) - bNext) / (bgrid(j+1) - bgrid(j));
        v0 = w*v(j, iz) + (1 - w)*v(j+1, iz);
    end

end
