classdef SPOSCAR
    %SVOR Summary of this class goes here
    %Detailed explanation goes here
    
    properties
        Path=[];   %path
        Steps=[];   % number of iterations
        NumSubProblem = [];
        Direction=[];
        Time = [];
        range_zeta = [];
    end
    methods(Static = true)
        %-----------  add J
        function obj=SPOSCAR(x,y,d,range_zeta)
            global KTYPE KSCALE fake_zero size_training  NumSingular epsilon_SP
            % KTYPE = 6;
            fake_zero=10^-8;
            fake_zero2=10^-8;
            NumSingular=0;
%             epsilon_SP = 1.1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            out=SPOSCAR.SolutionPath(x,y,d,range_zeta);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Steps=out.Steps;
            obj.Path = out.Path;
            obj.Time = out.Time;
            obj.Direction =d;
            obj.range_zeta= range_zeta;
            obj.NumSubProblem = out.NumSubProblem;
        end
        function out_sp=SolutionPath(x,y,d,range_zeta)
            global fake_zero time_batch  epsilon epsilon_SP
            tic;
            path=[];
            NumSubProblem = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%core code
            current_zeta = range_zeta(2);
            epsilon = 0.05;
            epsilon_SP = 0.1;
            old_beta = zeros(length(x(1,:)),1);
            L = 50000;
            [objective_function,b] = OSCAR.ComputeObjective(x,y,old_beta,old_beta,d,L);
            epsilon = epsilon*objective_function;
            epsilon_SP = epsilon_SP*objective_function;
            tic;
            while current_zeta > range_zeta(1)+fake_zero
%                 local_flag=1;
%                 while local_flag==1
                    [initial_solution,out_time]=SPOSCAR.Initial(x,y,d,current_zeta);
%                     if sum(abs(initial_solution.beta))<0.5 | initial_solution.Groups.num_groups<2
%                         current_zeta = 0.9*current_zeta;
%                     else
%                         local_flag=0;
%                     end
%                 end
                time_batch = out_time;
                [local_path]=SPOSCAR.Run(initial_solution,x,y,d);
                path=[path;local_path];
                current_zeta = local_path(end).zeta;
                NumSubProblem = NumSubProblem+1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc;
            local_time=toc;
            steps=length(path);
            out_sp.Steps=steps;
            out_sp.Path=path;
            out_sp.Time=local_time;
            out_sp.NumSubProblem = NumSubProblem;
        end
        function [out_path] = Run(initial_solution,x,y,d)
            global fake_zero max_iterations
            dim = length(x(1,:));
            path=[initial_solution];
            mystruct=initial_solution;
            [Xh, partition, label] = SPOSCAR.PreparXh(x,initial_solution,dim);
%             out_KKT1=SPOSCAR.CheckOptimal(mystruct,partition,Xh,y,dim,d);
            steps=0;
            flag_running = 1;
            while mystruct.zeta > fake_zero*100 & flag_running==1
                xi=SPOSCAR.GetXi(Xh,y,mystruct,partition,dim,d);
                [max_delta,index,flag_max]=SPOSCAR.GetMaxChange(mystruct,xi,Xh,y,d,partition);
                [mystruct,partition]=SPOSCAR.Updata(mystruct,max_delta,index,flag_max,xi,Xh,y,d,partition);
%                 out_KKT1=SPOSCAR.CheckOptimal(mystruct,partition,Xh,y,dim,d);
%                 mystruct.KKT=out_KKT1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%accessorial code
                [flag_running,difference,beta] = SPOSCAR.CheckDetermine(x,y,mystruct,d,label);
                mystruct.beta = beta;
                path=[path;mystruct];
                steps=steps+1;
                if steps>max_iterations
                    error('a');
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            out_path = path;
        end
        function [running_flag,difference,beta] = CheckDetermine(x,y,mystruct,d,label)
            global epsilon_SP
            running_flag = 0;
            beta = SPOSCAR.RecoverBeta(mystruct,length(x(1,:)),label);
            pms = d*mystruct.zeta;
            [F_objective, Qobject]= OSCAR.ComputeObjective(x,y,beta,beta,pms,0.1)
            [running_flag,difference] = OSCAR.CheckDetermine(x,y,beta,F_objective,epsilon_SP,pms);
            if difference<epsilon_SP
                running_flag = 1;
            else
                running_flag = 0;
            end
        end
        function beta = RecoverBeta(mystruct,dim,label)
            Groups = mystruct.Groups;
            groups = Groups.groups;
            num_groups = Groups.num_groups;
            beta_group_abs = Groups.beta_group_abs';
            beta = zeros(dim,1);
            for i=1:num_groups
                local_group = groups{i};
                index_tmp = local_group(:,2);
                beta(index_tmp) = beta_group_abs(i);
            end
            beta = beta.*label;
        end
        function [Xh, partition,label] = PreparXh(x,initial_solution,dim)
            beta = initial_solution.beta;
            Groups = initial_solution.Groups;
            groups = Groups.groups;
            num_groups = Groups.num_groups;
            Xh = [];
            label = zeros(dim,1);
            for i=1:num_groups
               index = groups{i}(:,2);
               local_beta = beta(index);
               local_label = (1*(local_beta>0)+(-1)*(local_beta<=0));
               label(index) = local_label;
               local_x = x(:,index)*local_label;
               Xh = [Xh local_x];
            end
            active = Groups.active;
            total_index = [1:num_groups];
            total_index(active) = [];
            partition{1}=active;
            partition{2}=total_index;
        end
        function [out_mystruct,partition] = Updata(mystruct,max_delta,index,flag_max,xi,Xh,y,d,partition)
            global epsilon fake_zero 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%update zeta, beta and v
            zeta = mystruct.zeta;
            groups = mystruct.Groups;
            zeta = zeta + max_delta;
            index_A = partition{1};
            old_beta = mystruct.Groups.beta_group_abs';
            beta = mystruct.Groups.beta_group_abs';
            beta(index_A) = beta(index_A) + xi*max_delta;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%update partition and order
            if flag_max==1
                local_index = find(partition{1}==index);
                partition{1}(local_index) = [];
                partition{2}=sort([partition{2} index]);
            end
            if flag_max==3
                tmp1 = groups.groups{index(1)};
                groups.groups{index(1)} = groups.groups{index(2)};
                groups.groups{index(2)} = tmp1;
                
            end
            mystruct.Groups.beta_group_abs = beta';
            out_mystruct.zeta = zeta;
            out_mystruct.Groups = groups;
            out_mystruct.pms = d*zeta;
            out_mystruct.beta = mystruct.beta;
        end
        function index_maxzero = FindIndexMaxZero(index_A_hat,order)
            local_order = order(index_A_hat);
            [max_order, index_max_order] = max(local_order);
            index_maxzero = index_A_hat(index_max_order);
        end
        function [max_delta,index,flag_max] = GetMaxChange(mystruct,xi,x,y,d,partition)
            [change_zeta1 index_change_zeta1]=SPOSCAR.FindChangeZeta1(mystruct,xi,partition);
            [change_zeta3 index_change_zeta3]=SPOSCAR.FindChangeZeta3(mystruct,xi,partition);
            [max_delta index,flag_max]=SPOSCAR.FindMaxChange(change_zeta1, index_change_zeta1,...
                change_zeta3,index_change_zeta3,mystruct);
        end
        function [change_zeta1 index_change_zeta1] = FindChangeZeta1(mystruct,xi,partition)   
            global fake_zero
            change_zeta1=inf;
            index_change_zeta1 = inf;
            index_A = partition{1};
            beta = mystruct.Groups.beta_group_abs';
            local_beta = beta(index_A);
            local_xi = xi;
            max_change_zeta = -inf;
            index_max_change_zeta = -inf;
            local_index_A=index_A;
            if length(index_A)>0
                local_index=find(local_xi<fake_zero & local_xi>-fake_zero);
                local_beta(local_index)=[];
                local_xi(local_index)=[];
                local_index_A(local_index)=[];
                local_index = find(local_xi.*local_beta<=0);
                local_beta(local_index)=[];
                local_xi(local_index)=[];
                local_index_A(local_index)=[];
                local_change_al = - local_beta./local_xi;
                [max_alpha_change max_alpha_change_index]=max(local_change_al);
                if length(local_change_al)==0
                    change_zeta1 = -inf;
                    index_change_zeta1 = -inf;
                else
                    change_zeta1=max_alpha_change;
                    index_change_zeta1 = local_index_A(max_alpha_change_index);
                end
            end
        end
        function [change_zeta3 index_change_zeta3] = FindChangeZeta3(mystruct,xi,partition)   
            global fake_zero
            index_A = partition{1};
            beta = mystruct.Groups.beta_group_abs';
            local_beta = beta(index_A);
            mylength = length(index_A);
            beta1 = local_beta(1:mylength-1);
            beta2 = local_beta(2:mylength);
            xi1 = xi(1:mylength-1);
            xi2 = xi(2:mylength);
            A = [xi2- xi1];
            B= [beta1-beta2];
%             [max2,fval,exitflag,output]=linprog(1,A,B,[],[],-inf,0);
            [max_zeta_change max_zeta_change_index] =...
                SPOSCAR.MyIneqaulity( A,B,mylength-1);
            if max_zeta_change==-inf
            	change_zeta3 = -inf;
            	index_change_zeta3 = -inf;
            else
                change_zeta3 = max_zeta_change;
                index_change_zeta3 = [max_zeta_change_index max_zeta_change_index+1];
            end
        end
        function [max_zeta max_zeta_index] = MyIneqaulity( A,B,mylength)
            index = find(A<0 & B>0);
            if length(index)>0
                [max_value index_value]=max(B(index)./A(index));
                max_zeta = max_value;
                max_zeta_index = index(index_value);
                max_zeta_index=mod(max_zeta_index,mylength);
                if max_zeta_index==0
                    max_zeta_index = mylength;
                end
            else
                max_zeta = -inf;
                max_zeta_index = [];
            end
        end
        function [max_delta index,flag_max] = FindMaxChange(change_zeta1, index_change_zeta1,...
                change_zeta3,index_change_zeta3,mystruct)
            max_delta=-inf;
            index=0;
            flag_max=0;
            if max_delta < change_zeta1
                max_delta=change_zeta1;
                index=index_change_zeta1;
                flag_max = 1;
            end
            if max_delta < change_zeta3
                max_delta = change_zeta3;
                index = index_change_zeta3;
                flag_max = 3;
            end
            if max_delta + mystruct.zeta<0
                max_delta = -mystruct.zeta;
                index = 0;
                flag_max = 0;
            end
        end
        function [xi] = GetXi(Xh,y,mystruct,partition,dim,d)
            global fake_zero NumSingular
            active_set = partition{1};
            H_AA = Xh(:,active_set)'*Xh(:,active_set);
            groups = mystruct.Groups;
            mylength = length(active_set);
            wg = SPOSCAR.ComputeWg(mylength,active_set,groups,dim,d);
            xi=linsolve(2*H_AA,-wg);
            tmp=det(H_AA);
            if tmp<fake_zero & tmp>-fake_zero
                NumSingular=NumSingular+1;
            end
        end
        function wg = ComputeWg(mylength,active_set,groups,dim,d)
            wg =[];
            totoal_index =0;
            for i=1:mylength
                index_group = active_set(i);
                leng_one_group = length(groups.groups{i}(:,1));
                tmp1 = (totoal_index+1+totoal_index+leng_one_group)*leng_one_group/2;
                tmp2 = d(2)*(leng_one_group*dim-tmp1) + d(1)*leng_one_group;
                totoal_index = totoal_index+leng_one_group;
                wg = [wg;tmp2];
            end
        end
        function [initial_solution,out_time]=Initial(x,y,d,zeta)
            pms = d*zeta;
            os = OSCAR(x,y,pms)
            initial_solution.beta = os.beta;
            initial_solution.Groups = os.Groups;
            initial_solution.pms = os.pms;
            initial_solution.zeta = zeta;
            out_time=os.Time;
        end
        function upper = ComputUpper(training_x,training_y,dir)
            d= dir;
            local_tmp1 = abs(training_x'*training_y);
            [ranked_out,local_rank]=sort(local_tmp1,'ascend'); 
            [a,local_rank] = sort(local_rank,'ascend');
%             local_rank = length(training_x(1,:))+1 - local_rank;
            local_tmp2 = d(1)+d(2)*(local_rank-1);
            local_tmp3 = local_tmp2.^-1;
            local_tmp4 = local_tmp3.*local_tmp1;
            [upper] = max(local_tmp4);
        end
        function  KKTflag = CheckOptimal(mystruct,partition,Xh,y,dim,d)
            global fake_zero
            index_A = partition{1};
            index_A_hat = partition{2};
            beta = mystruct.Groups.beta_group_abs';
            local_beta = beta(index_A);
            H_AA = Xh(:,index_A)'*Xh(:,index_A);
            groups = mystruct.Groups;
            mylength = length(index_A);
            wg = SPOSCAR.ComputeWg(mylength,index_A,groups,dim,d);
            equ1 = 2*H_AA*local_beta- 2*Xh(:,index_A)'*y + mystruct.zeta*wg;
            tmp1 = equ1 > fake_zero | equ1 < -fake_zero;
            tmp2 = sum(tmp1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check the order
            if tmp2>0 
                KKTflag = 0;
            else
                KKTflag = 1;
            end
        end
    end
end

