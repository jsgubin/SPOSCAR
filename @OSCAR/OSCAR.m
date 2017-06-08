classdef OSCAR
    %SVOR Summary of this class goes here
    %Detailed explanation goes here
    
    properties
        beta = [];   %path
        Groups = {};   % number of iterations
        pms = [];
        Time = [];
        Steps =[];
    end
    methods(Static = true)
        %-----------  add J
        function obj = OSCAR(x,y,pms)
            global fake_zero epsilon_active max_active epsilon_SP epsilon
            fake_zero=10^-10;
            max_active = 40
            epsilon_active =0.01;
%             L = 1000000; % Lipschitz constant
%             L = 300000; % Lipschitz constant
            L = 50000; % Lipschitz constant
            running_flag = 1;
            zeta = 1.1;
            old_tao = 1;
            t = 1;
            old_beta = zeros(length(x(1,:)),1);
            yyy = old_beta;
            tic;
            objective_function = [];
            [objective_function(t),b] = OSCAR.ComputeObjective(x,y,old_beta,old_beta,pms,L);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            while running_flag == 1
                L_flag =1;
                while L_flag == 1
                    [new_beta,groups] = OSCAR.SolveSubProglem(x,y,yyy,L,pms);
                    [Fobject,Qobject]=OSCAR.ComputeObjective(x,y,new_beta,yyy,pms,L);
                    if Fobject <= Qobject
                        L_flag = 0;
                    else
                        L = L*zeta;
                    end
                end
                new_tao = (1+sqrt(1+4*old_tao^2))/2;
                yyy = new_beta+ ((old_tao-1)/new_tao)*(new_beta-old_beta);
                t = t+1;
                if t>5000
                    a=1;
                end
                objective_function(t) = Fobject;
                [running_flag,difference] = OSCAR.CheckDetermine(x,y,new_beta,objective_function(t),epsilon,pms);
            end
            toc;
            OSCAR.PlotEnergy(objective_function);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Steps = t;
            obj.beta = new_beta;
            obj.Groups = groups;
            obj.Time = toc;
            obj.pms = pms;
        end
        function PlotEnergy(objective_function)
            mylength = length(objective_function);
            xx=[1:mylength];
            plot(xx,objective_function);
        end
        function [F_objective, Qobject]= ComputeObjective(x,y,new_beta,yyy,pms,L)
            mylength = length(x(1,:));
            primal_1 = norm(x*new_beta-y)^2;
            beta_abs = abs(new_beta);
            [order_beta_abs, index_order_gamma] = sort(beta_abs,'descend');
            order = [1:mylength];
            out_second_term = sum((pms(1) + pms(2)*(mylength - order))'.*order_beta_abs);
            F_objective = primal_1 + out_second_term;
            deriv_f = x'* (2*(x*yyy - y));
            primal_21 = norm(x*yyy-y)^2;
            primal_22 = (new_beta - yyy)'*deriv_f + (L/2)*norm(yyy - new_beta)^2;
            Qobject = primal_21+primal_22+out_second_term;
        end
        function [beta,Group] = SolveSubProglem(x,y,beta_input,L,pms)
            global fake_zero epsilon_active max_active
            mylength = length(x(1,:));
            deriv_f = x'* (2*(x*beta_input - y));
            aa = L*beta_input - deriv_f;
            a=(2/L)*abs(aa);
            [order_a, index_order_a] = sort(a,'descend');
            for i=1:mylength
                InitialGroup{i} = [i index_order_a(i)];
            end
            mystack{1} = InitialGroup{1};
            mytop = 1;
            for i = 2:mylength
                g = InitialGroup{i};
                v_top = OSCAR.GetV(mystack{mytop}(:,1),order_a,mylength,L,pms);
                v_g = OSCAR.GetV(g(:,1),order_a,mylength,L,pms);
                while mytop > 0 & v_g >= v_top-fake_zero
                    g = [mystack{mytop};g]
                    mytop = mytop-1;
                    if mytop > 0
                        v_g = OSCAR.GetV(g(:,1),order_a,mylength,L,pms);
                        v_top = OSCAR.GetV(mystack{mytop}(:,1),order_a,mylength,L,pms);
                    end
                end
                mystack{mytop+1} = g;
                mytop = mytop+1;
            end
            beta_group_abs = [];
            out_second_term = 0;
            active = [];
            for i=1:mytop
                tmp_z = OSCAR.GetV(mystack{i}(:,1),order_a,mylength,L,pms);
                beta_group_abs(i) = max([tmp_z 0]);
                    if length(active)<max_active &  beta_group_abs(i)>epsilon_active
                        active = [active;i];
                    end
            end
            beta = zeros(mylength,1);
            for i=1:mytop
                local_index = mystack{i}(:,2);
                beta(local_index) = beta_group_abs(i);
            end
            local_label = (1*(aa>0)+(-1)*(aa<=0));
            beta = beta.*local_label;
            Group.groups = mystack;
            Group.num_groups = mytop;
            Group.active = active;
            Group.beta_group_abs = beta_group_abs;
        end
        function out_v = GetV(input,a,dim,L,pms)
            local_a = a(input);
            local_w = (2/L)*(pms(1) + pms(2)*(dim - input));
            out_v = sum(local_a - local_w)/(2*(input(end)-input(1)+1));
        end
        function  [running_flag,difference] = CheckDetermine(x,y,new_beta,objective_function,epsilon,pms)
            running_flag = 1;
            dim = length(x(1,:)); 
            primal = objective_function;
            deriv_f = 2*(x*new_beta - y);
            gamma_h = x'* deriv_f;
            abs_gamma_h = abs(gamma_h);
            [order_gamma_h, index_order_gamma] = sort(abs_gamma_h,'descend');
            order = [1:dim]
            lamda = pms(1) + (dim-order)*pms(2);
            sum_gama = zeros(dim,1);
            sum_lamda = zeros(dim,1);
            tmp1 = 0;
            tmp2 = 0;
            for i=1:dim
                sum_gama(i) = tmp1 + order_gamma_h(i);
                tmp1 = sum_gama(i);
                sum_lamda(i) = tmp2 + lamda(i);
                tmp2 = sum_lamda(i);
            end
            tmp = sum_gama./sum_lamda;
            r_star = max(tmp);
            first_term = min([1,1/r_star]);
            alpha =  first_term*deriv_f;
            dual = 0.25 * norm(alpha)^2 + alpha'* y;
            difference = primal+ dual;
            if difference<epsilon
                running_flag = 0;
            end
        end
    end
end

