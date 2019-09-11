function biot_examples()
    
    format long;
    clc; 
    clf;
   
    syms x y z t u p alph bet mu lmbda g f
    
    % beta * p_t + alpha * div (u_t) - div ( K*grad_p) = g
    % -div( lambda * div(u) * I + 2 * mu * eps(u) - alpha * p * I) = f
    
    examples = [5];
    
    for i = 1 : length(examples)
        example_num = examples(i);
        if example_num == 1 % manual_example_2d_t
            d = 2;
            
            p = x*(1 - x)*y*(1 - y)*t;
            
            u = [x*(1 - x)*y*(1 - y)*t;
                 x*(1 - x)*y*(1 - y)*t];
            
            % parameters
            K = [1 0;
                 0 1];
            %{
            alph = 1;
            bet  = 1;
            
            mu = 1;
            lmbda = 1;
            %}
        elseif example_num == 2 % simple_example_2d_t
            d = 2;
            p = x*(1 - x)*y*(1 - y)*t;
            
            u = [(x^2 + y^2)*t;
                 (x + y)*t];
            
            % parameters
            K = [1 0;
                 0 1];
            
            %{
            alpha = 1;
            beta  = 1;
            
            mu = 1;
            lmbda = 1;
            %}
        elseif example_num == 3
            d = 2;
            p = sin(pi*x)*sin(3*pi*y)*(t^2+t+1);
            
            u = [(x^2 + y^2)*sin(t);
                 (x + y)*exp(t)];
            
            % parameters
            K = [1 0;
                 0 1];
        elseif example_num == 4
            d = 3;
            p = x*(1 - x)*y*(1 - y)*z*(1 - z)*t;
            
            u = [(x^2 + y^2)*t;
                 (x + y)*t;
                 sin(t)];
            
            % parameters
            K = [1 0 0;
                 0 1 0;
                 0 0 1];
             
        elseif example_num == 5
            d = 2;
            p = x*(1 - x)*y*(1 - y);
            
            u = [x*(1 - x)*y*(1 - y);
                 x*(1 - x)*y*(1 - y)];
            
            % parameters
            K = [1 0;
                 0 1];
            alph = 1;
            lmbda = 1;
            mu = 1;
            bet = 1;
        elseif example_num == 6 % 
            d = 2;
            p = 10^8*x*(1 - x)*y*(1 - y)*t;
            
            u = [(x^2 + y^2)*t;
                 (x + y)*t];
            
            % parameters
            K = [100 0;
                 0 100];
            
            %{
            alpha = 1;
            beta  = 1;
            
            mu = 1;
            lmbda = 1;
            %}
        end
        
        I = eye(d);
    
        fprintf('p = \n'); disp(p);
        fprintf('u = \n'); disp(u);
        
        % p
        if d == 2
            vars = [x, y];
            grad_p  = [simplify(diff(p, x)); 
                        simplify(diff(p, y))];
        elseif d == 3
            vars = [x, y, z];
            grad_p  = [simplify(diff(p, x)); 
                        simplify(diff(p, y));
                        simplify(diff(p, z))];
        end
                
        p_t = diff(p, t);
        
        fprintf('grad_p = \n'); disp(grad_p);
        fprintf('p_t = \n'); disp(p_t);
        
        % u
        grad_u = simplify(jacobian(u, vars));
        eps_u = simplify(1/2 * (grad_u + transpose(grad_u)));
        div_u = simplify(divergence(u, vars));
        div_u_t = simplify(diff(div_u, t));

        fprintf('grad_u = \n'); disp(grad_u);
        fprintf('eps_u = \n'); disp(eps_u);
        fprintf('div_u = \n'); disp(div_u);
        fprintf('div_u_t = \n'); disp(div_u_t);
        
        div_K_gradp = simplify(divergence(K * grad_p, vars));                     % diffusion part of the flow equation
        sigma_por = simplify(lmbda * div_u * I + 2 * mu * eps_u - alph * p * I); % lambda * div(u) * I + 2 * mu * eps(u) - alpha * p * I
        
        fprintf('div_K_gradp = \n'); disp(div_K_gradp);
        fprintf('sigma_por = \n'); disp(sigma_por);
        
        g = simplify(bet * p_t + alph * div_u_t - div_K_gradp);
        f = simplify(- Divergence(sigma_por, vars));
        fprintf('g = \n'); disp(g);
        fprintf('f = \n'); disp(f);
        
    end
end

function result = Divergence(tensor, vars)
    dim = length(vars);
    %result = zeros(length(vars), 1);
    result = sym(zeros(length(vars), 1));
    for i = 1 : dim
        result(i, 1) = divergence(tensor(i, :), vars);
    end
end

