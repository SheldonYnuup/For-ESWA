clear all;

clc;

Function_name='F1'; % Name of function

SearchAgents_no=30; % Number of search agents

Max_iteration=1000; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=function_detail(Function_name); 

% main
[Best_AOBLMOA,AOBLMOA_cg_curve,AOBLMOA_pos,position_history,fitness_history,Trajectories]=AOBLMOA(SearchAgents_no,SearchAgents_no,SearchAgents_no,Max_iteration,lb,ub,dim,fobj);  

% figure plot
figure('Position',  [680   400   560   560])
semilogy(1:Max_iteration,AOBLMOA_cg_curve,'DisplayName','AOBLMOA','Color',[0.18 0.54 0.34],'LineWidth',1.5)
hold on
legend('AOBLMOA');
title('Convergence curve')
xlabel('Iteration');
ylabel('Best score obtained so far');
box on
axis tight
hold on

function [lb,ub,dim,fobj] = function_detail(Function_name)
switch Function_name
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim = 10;      
end

end
function o = F1(x)
o=sum(x.^2);
end

