function grad = objective_fd(fun, T, isnormalized)

if nargin < 3
    isnormalized = [];
end

if isempty(isnormalized)
    isnormalized = false;
end

grad = zeros(numel(T),1);
for i = 1:numel(T)
    T_d = T;
    if isnormalized
        T_d(i) = T_d(i) + T_d(i)*1e-8i;
    else
        T_d(i) = T_d(i) + 1e-8i;
    end
    fun_d = fun(T_d);
    grad(i) = imag(fun_d)/1e-8;
end