clear

rng('default');
v = @(q,k) 0 + q.*randn(1,k);        % ~N(0,q)
u = zeros(1,50)+1;
w = @(r,k) 0 + r.*randn(1,k);        % ~N(0,r)

%%
f = 1;
h = 1;
q = 0.01;
r = 1;
n = 50;

v = v(q,n);
w = w(r,n);
x(1) = 0;   % x0
z(1) = w(1);

result(1) = x(1);
for k = 2:n
    x(k) = f*x(k-1) + u(k-1) + v(k-1); % x_{k+1}
    z(k) = h*x(k) + w(k);
end

%%
