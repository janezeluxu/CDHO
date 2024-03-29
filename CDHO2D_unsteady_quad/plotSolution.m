function plotSolution(fileName,solution,plotnum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misc. Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
global node_Size;
 
[OrderList,meshData,vertexData,Area,...
    IBC,IBC_vertex,BCval,BCb,solution_ien,order_list,edge_usage,...
    iper,uHBCnode,masternodes,slavenodes] ...
    = buildMeshStruct(fileName);
 writeMesh(meshData,vertexData,'boundary')
 %'./testvtk/cases/boundary/solution/p2dynamic.txt'
 u=load(solution);
 %[x,y,ien] = boundary();
 %sp = trimesh(ien,x,y,u(1:node_Size));
%  [X,Y] = meshgrid(0:0.5:1, 0:0.5:1);
%  Z = zeros(3,3);
%  sp = mesh(X,Y,Z);

%% mesh plot settings
% set(sp,'EdgeColor',[0 0 0],'LineWidth',2);
%   axis([0 1 0 1])
%   view(-45,90)
%   xL = xlim;
% yL = ylim;
% line([0 0], yL);  %x-axis
% line(xL, [0 0]);  %y-axis
% set(gca,'visible','off')
% set(gcf,'units','points','position',[100,100,500,500])
% 
% figure1 = figure(1)
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.168000000000001 0.04 0.0759999999999997 0.046],'String',{'x'},...
%     'FontSize',36,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% 
% % Create arrow
% annotation(figure1,'arrow',[0.13 0.222],[0.109 0.108]);
% 
% % Create arrow
% annotation(figure1,'arrow',[0.134 0.226],...
%     [0.519000000000001 0.518000000000001]);

% % Create textbox
% annotation(figure1,'textbox',...
%     [0.0760000000000007 0.164 0.0879999999999993 0.078000000000001],...
%     'String','y',...
%     'FontSize',36,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);

% % Create textbox
% annotation(figure1,'textbox',...
%     [0.186000000000001 0.404 0.0879999999999994 0.078000000000001],...
%     'String','Flow',...
%     'FontSize',36,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);


% annotation(figure1,'textbox',...
%     [0.588000000000001 0.046 0.0759999999999996 0.046],'String',{'x'},...
%     'FontSize',30,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.456000000000001 0.202 0.0879999999999992 0.0780000000000011],...
%     'String','y',...
%     'FontSize',30,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% 
% % Create arrow
% annotation(figure1,'arrow',[0.518000000000001 0.610000000000002],...
%     [0.111 0.11]);
% 
% % Create arrow
% annotation(figure1,'arrow',[0.516 0.516],[0.111 0.22]);
% 
% % Create arrow
% annotation(figure1,'arrow',[0.328 0.424],[0.317000000000001 0.42]);
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.278000000000003 0.4 0.0999999999999969 0.0580000000000015],...
%     'String','Flow',...
%     'FontSize',30,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);

% %  %ylabel('Exact','FontSize',30)
% %ylabel('\tau_{alg}','FontSize',40)
%  ylabel('\tau_{dyn2}','FontSize',40)
% %ylabel('Inclined \tau_{dyn2}','FontSize',40)
%  set(sp,'EdgeColor',[0 0 0]);
%  %set(sp,'EdgeColor',[0 0 0],'LineWidth',2);
%  %set(gca,'visible','off')
 hold on
 polyplot(meshData,solution_ien,order_list,edge_usage,vertexData,u,plotnum);
 %set(sp,'EdgeColor',[1 1 1]);
 
 %axis off
 %h = colorbar('eastoutside','Ticks',[0,0.2,0.4,0.6,0.8,1.0]);
 %h.Limits = [0 1.0];
 %set(gca, 'FontSize', 20)

%  axis([0 1 0 1])
%  view(0,90)
%  xticks([0.0 0.5 1.0])
%  yticks([ 0.0 0.5 1.0])
%  zticks([ 0.0 0.5 1.0])
%  xlim([0 1.0])
%  ylim([0 1.0])
%  zlim([-0.1 1.0])
%  set(gca, 'XTickLabel',{'0.0', '0.5', '1.0'},...
%      'YTickLabel',{ '0.0', '0.5', '1.0'},...
%      'ZTickLabel',{ '0.0', '0.5', '1.0'})
% set(gca, 'FontSize', 30)
 %% coutour plot
end