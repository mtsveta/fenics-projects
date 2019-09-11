function biot_symbolic_computations()
    syms x_0 x_1 t
    syms mu_f bet alph mu lmbda
    example_num = 2;
    
    switch example_num
        case 1
            % --------------------------------------------------------------------%
            % Example 1:
            % --------------------------------------------------------------------%
            K = [1 0; 0 1];
            p = t * x_0 * (1 - x_0) * x_1 * (1 - x_1);
            u = t * [x_0 * (1 - x_0) * x_1 * (1 - x_1); 
                     x_0 * (1 - x_0) * x_1 * (1 - x_1)];
            dim = 2;
        case 2
            % --------------------------------------------------------------------%
            % Example 2:
            % --------------------------------------------------------------------%
            K = [1 0; 0 1];
            p = t * x_0 * (1 - x_0) * x_1 * (1 - x_1);
            u = t * [x_0*x_0 + x_1*x_1; 
                     x_0 + x_1];  
            dim = 2;
    end
    
    [s_f, f] = calculate_source_and_load(K, p, u, dim);
    
    disp(['p = ', latex(p)]);
    disp(['u = ', latex(u)]);
    disp(['s_f = ', latex(s_f)]);
    disp(['f   = ', latex(f)]);
    
    %q = - 1 / mu_f * nabla_p;
end

function [s_f, f] = calculate_source_and_load(K, p, u, dim)
    
    syms x_0 x_1 t
    syms mu_f bet alph mu lmbda
    
    args = [x_0, x_1];
    
    nabla_p = simplify(gradient(p, args));
    p_t = simplify(diff(p , t, 1));
    div_K_nabla_p = simplify(divergence(K * nabla_p, args));
    
    %nadla_u = [diff(u(1), x0, 1) diff(u(1), x1, 1); diff(u(2), x0, 1) diff(u(2), x1, 1)];
    nabla_u = grad_tensor(u, args, dim); %simplify([transpose(gradient(u(1), [x0, x1])); transpose(gradient(u(2), [x0, x1]))]);
    eps_u = simplify(0.5 * (nabla_u + transpose(nabla_u)));
    div_u = simplify(divergence(u, args));
    div_u_t = simplify(diff(div_u, t, 1));
    I = unit(dim);
    
    sigma_por = simplify(2 * mu * eps_u + lmbda * div_u * I - alph * p * I);
    f = simplify( - div_tensor(sigma_por, args, dim)); %divergence(sigma_por, [x0, x1]));
    s_f = simplify(bet * p_t - 1 / mu_f * div_K_nabla_p + alph * div_u_t);
end 

function I = unit(dim)
    switch dim
        case 2
            I = [1 0; 0 1];
        case 3
            I = [1 0 0; 0 1 0; 0 0 1];
    end
end

function grad_u = grad_tensor(u, args, dim)
    switch dim
        case 2
            grad_u = simplify([transpose(gradient(u(1), args)); transpose(gradient(u(2), args))]);
        case 3
            grad_u = simplify([transpose(gradient(u(1), args)); transpose(gradient(u(2), args)); transpose(gradient(u(3), args))]);
    end
        
end

function grad_u = div_tensor(u, args, dim)
    switch dim
        case 2
            grad_u = simplify([divergence(u(1, :), args); divergence(u(2, :), args)]);
        case 3
            grad_u = simplify(divergence([u(1, 1), u(2, 2), u(3, 3)], args));
    end
        
end
