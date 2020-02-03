function plotSolution(fileName,solution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Misc. Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
global node_Size;
 
 [OrderList,meshData,vertexData,Area,IBC,IBC_vertex,BCval,BCb,uHBCE,solution_ien,...
    order_list,edge_usage,iper] = preprocess(fileName);
 writeMesh(meshData,vertexData,'boundary')
 %'./testvtk/cases/boundary/solution/p2dynamic.txt'
 u=load(solution);
 %[x,y,ien] = boundary();
 %sp = trimesh(ien,x,y,u(1:node_Size));
 %set(sp,'EdgeColor',[0 0 0]);
%hold on 
polyplot(meshData,solution_ien,order_list,edge_usage,vertexData,u);
 view(0,90)
 axis([0 1 0 1 0 1.2])
% x1 = 0.15; x2 = 0.25;
% y1 = 0.15; y2 = 0.225;
% 
% % P1 = [x1,y1];
% % P2 = [x2,y1];
% % pts = [P1; P2];
% % plot(pts(:,1), pts(:,2),'r-','MarkerSize',8,'MarkerFaceColor','r','LineWidth',2)
% 
%  axis([0 1 0 0.2])
%  xticks([0.0 0.2 0.4 0.6 0.8 1.0])
%  yticks([0.0 0.1 0.2])
%  a = get(gca,'XTickLabel');
%  set(gca,'fontsize',20)
% %grid on
% %grid minor
% set(gca,'units','points','Position',[50,50,900,300])
% figure1 = figure(1);
% % Create textbox
% annotation(figure1,'textbox',[0.212 0.8975 0.044 0.0640000000000001],...
%     'String',{'0.2'},...
%     'FontSize',20,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% ax1 = gca();
% ax1.TickDir = 'out';
% %% Create axes 2
% figure1 = figure(1);
% axes2 = axes('Parent',figure1,...
%     'Position',[0.6 0.5375 0.348 0.46]);
% hold(axes2,'on');
% %% Create plot
% polyplot(meshData,solution_ien,order_list,edge_usage,vertexData,u);
% %axis ([0.15 0.25 0.15 0.225])
% axis ([x1 x2 y1 y2])
% xticks([0.2])
% %yticks([0.15  0.2])
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% box(axes2,'on');
% 
% P1 = [0.15,0.2];
% P2 = [0.25,0.2];
% pts = [P1; P2];
% plot(pts(:,1), pts(:,2),'k-','MarkerSize',18,'LineWidth',0.5)
% hold on
% 
% P1 = [0.2,0.2+0.003];
% P2 = [0.2,0.2];
% pts = [P1; P2];
% plot(pts(:,1), pts(:,2),'k-','MarkerSize',18,'LineWidth',0.5)
% hold on
% 
% % Create textbox right
% annotation(figure1,'textbox',[0.756 0.8775 0.044 0.0640000000000001],...
%     'String',{'0.2'},...
%     'FontSize',20,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% % Create arrow
% annotation(figure1,'arrow',[0.343 0.534],[0.924 0.9225],'Color',[1 0 0],...
%     'LineWidth',5,'HeadLength',30,'HeadWidth',30);
% % Create rectangle
% annotation(figure1,'rectangle',[0.188 0.685 0.086 0.2825],'Color',[1 0 0],...
%     'LineWidth',5);
% % Create rectangle on right
% annotation(figure1,'rectangle',[0.6 0.5375 0.349 0.46],'Color',[1 0 0],...
%     'LineWidth',5);
% set(gca,'fontsize',20)
% set(groot,'defaultFigurePaperPositionMode','auto')
% set(gcf,'units','points','position',[100,100,1000,400])

 %% coutour plot
end