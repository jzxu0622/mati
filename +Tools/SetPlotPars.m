%% This script is for the convenience of plotting
% 1. Defines some MATLAB visualization parameters for convenience
% 2. Show examples how to plot a few different kinds of figures
% -------------------------------------------------------------------------------------------------------------------------
% NOTES: 
%       1. Please use this script after a figure opened
%       2. This is a script so all defined parameters will go into main script
% -------------------------------------------------------------------------------------------------------------------------


%% -------------------- figure size -------------------------
% MRM author guidelines
% 1 column = 3.42 inches or (8.67-cm) wide
% 1.5 column = 5.125 inches or (13.02-cm) wide
% 2 columns = 6.92 inches or (17.15-cm) wide. 
fig_width = 3.42 ;  % inches
fig_height = 3 ;  % inches


%% -------------------- Parameter names ------------------------
LW = 'LineWidth' ; 
lw.line = 2 ; 
FS = 'FontSize' ; fs.label = 12 ; fs.axis = 9 ; fs.legend=16;
FN = 'FontName' ; fn = 'Times New Roman' ; fn1 = 'Arial' ; 
FW = 'FontWeight' ; 
MS = 'MarkerSize' ; ms = 7 ; 
MFC = 'MarkerFaceColor' ; mfc = [1 1 1]*0.4 ; 
MEC = 'MarkerEdgeColor' ; mec = [1 1 1]*0.4 ; 

%  USAGE: plot(x, y, char(fig_marker(n)), MFC,char(fig_color(n)))
fig_marker_line = {'bo-','r>-','gd-','c^-','mp-','ks-'} ; 
fig_marker = {'bo','r>','kd','c^','mp','gs'} ; 
fig_color = {'b','r','k','c','m','g'} ; 

%% color
% transpancy
% col = [0 1 0] ; color_saturation = 0.3 ; color_new = col + (1-col)*(1-color_saturation) ; 

%% --------- set Nrow X Ncol ------------    
% figure(1) ;  clf ; set(gcf, 'Units', 'inches', 'Position', [9 6 6 5]) ; %colormap gray
% Nrow = 5 ; Ncol = 4 ; 
% xmin = 0.02 ; ymin = 0.05 ; lx = 0.18 ; ly = 0.10 ; dx = 0.02 ; dy = 0.08 ; [nx,ny]=meshgrid(1:Ncol,1:Nrow) ; 
% for n=1:length(nx(:))
%     hax1(n) = axes('Position',[xmin+(nx(n)-1)*(lx+dx), ymin+(Nrow-ny(n))*(ly+dy), lx, ly]) ; 
% end


%% set up figure  for presentation
% set(gcf, 'Units','inches', 'PaperUnits', 'inches', 'Position', [9 6 fig_width fig_height])
% set(gcf, 'color','none'); % 'white','none'
set(gcf,'renderer','painters')      % vector
set(groot,'defaultAxesFontName', 'Times New Roman') ; set(0,'DefaultTextFontName', 'Times New Roman') ; 

set(groot, ...
    'DefaultFigurePaperPositionMode','manual', ...                                              % --------- Paper -------
    'DefaultAxesFontSize', fs.axis, ...                                                                        % -------Axes -------
    'DefaultAxesFontName', 'Times New Roman' , ... % 'New Times', 'Arial', 'Helvetica' , 'New Times'% ...
    'DefaultAxesXColor','k', ...
    'DefaultAxesYColor','k', ...
    'DefaultAxesZColor','k', ...
    'DefaultAxesLineWidth', 0.5, ...
    'DefaultAxesFontSize', 12, ...
    'DefaultAxesFontName', 'Times New Roman', ... % 'New Times', 'Arial', 'Helvetica' 
    'DefaultAxesFontWeight', 'normal', ... % 'bold','normal'
    'DefaultLineMarkerSize', ms, ...                                                                    % -------Line -------
    'DefaultLineMarkerFaceColor','none', ... % 'none','r'
    'DefaultTextFontName','Times New Roman', ... % 'New Times', 'Arial', 'Helvetica'                    % -------Text -------
    'DefaultTextFontWeight','normal', ... % 'bold','normal'     
    'DefaultTextMargin', 1, ...
    'DefaultTextFontSize', fs.label, ...
    'DefaultUicontrolFontName', 'Arial', ...
    'DefaultUitableFontName', 'Arial', ...
    'DefaultUipanelFontName', 'Arial' ...
    ); 


% ...%'DefaultFigurePaperPositionMode','manual', ...
% ...
% ----- Axes ------
%     'DefaultAxesLineWidth', 0.5, ...
%     'DefaultAxesFontSize', 9, ...
%     'DefaultAxesFontName', 'Arial', ... % 'New Times', 'Arial', 'Helvetica' 
%     'DefaultAxesFontWeight', 'normal', ... % 'bold','normal'
% ...%    'DefaultAxesOuterPosition',[0.5 0.5 0.9 0.9], ...
% ...%    'DefaultAxesOuterPosition',[0.05 0.05 1 1], ...
% ...%    'DefaultAxesPosition',[0.15 0.15 0.8 0.8], ...
% ...%     'DefaultAxesPosition',[0.15 0.15 0.8 0.8], ...
% ...%'DefaultAxesTickLength',[0.01 0.025], ...     
% ...%'DefaultAxesTickDir','in', ...  
% ...%'DefaultAxesColorOrder',[0 0 1; 0 1 1; 1 1 0; 1 0 0; 0.75 0 0.75; 0.75 0.75 0; 0.25 0.25 0.25], ... 
% ...
% ... % line 
%     'DefaultLineLineWidth', lw, ...
% ...%'DefaultLineColor', cl, ...
%     'DefaultLineMarkerSize', ms, ...
%     'DefaultLineMarkerFaceColor','none', ... % 'none','r'
% ...%'DefaultLineMarkerEdgeColor','k', ...     
% ...    
% ... % text
%     'DefaultTextFontName','Arial', ... % 'New Times', 'Arial', 'Helvetica'
%     'DefaultTextFontWeight','normal', ... % 'bold','normal'     
%     'DefaultTextMargin', 1 ...
%   'DefaultTextFontSize', fs...

%% ---------------------------------------- errorbar --------------------------------
% ------ adjust separate lines in errorbars ------------
% hebxy=My_errorbarxy(mean(kin_highq), mean(AXR_all), std(kin_highq), std(AXR_all),{'ko-', 'b', 'r'}); 
% set(hebxy.hErrorbar(1,1), 'linewidth',1) 
% set(hebxy.hMain, 'markerfacecolor','g','markersize',6)

%% ---------------------------------------- colorbar ----------------------------------
% hcb1 = colorbar('location','EastOutside') ; 
% title('T_{1} [ms]') ; set(hcb1,'xlim',[1500 2550]) ; 
% length_cb = 1 ; width_cb = 0.3 ; position_tmp  = get(hcb1,'position') ; position_tmp(3) = position_tmp(3)*length_cb ; position_tmp(4) = position_tmp(4).*width_cb ; 
% set(hcb1, 'position', position_tmp) ; 


%% ---------------------------------- smooth parametric maps ------------------
% h=1/9*ones(3) ; % smoothing using a 3X3X1 kernel
% AxDmap = filter2(h, AxDmap) ; 


%% ------------------------ legend -------------------------------------
% edgecolor
% hlg1 = legend('intra-','extra-','both-') ; set(hlg1,'PlotBoxAspectRatio',[1 0.4 1]) ; 
% set(hlg1, FS,fs,'color','none','edgecolor','none', 'NumColumns',2, 'location','EastOutside')

% Greek letters
% legend(sprintf('\\alpha=\\beta'))

% add '\n' in legend
% ylabel(['DDR{\perp}',sprintf('\n'),'[\mum^2]'],FS,fs)

%% ----------------------------------- set up gray area ----------------------------
% % indicate zoom image position
% xzoom1 = 0 ; xzoom2 = 2 ; 
% yzoom1 = 0 ; yzoom2 = 2 ; 
% FCzoom = [1 1 1].*0.9 ; 
% % Create rectangle
% verts = [xzoom1, yzoom1 ; xzoom2, yzoom1 ; xzoom2, yzoom2 ; xzoom1, yzoom2] ; 
% faces = [1 2 3 4] ; 
% hph1 = patch('Faces',faces, 'Vertices',verts,'FaceColor', FCzoom,'edgecolor','none');


%% --------------------------------   set title and labels ----------------------------
% set(gca,'TickLabelInterpreter', 'tex'); % works for R2015a

%ylabel('$\bar{d}$fit [$\mu$m]',FS,fs.label,'Interpreter','Latex') ;  
% ylabel('R= $\displaystyle\frac{dADC}{df}$ [$\mu$m]','interpreter','latex')
% title('{\itvin} [%]')  %italic
% title('\itvin\rm [%]')  %italic

%% --------------------------- shaded error bars ----------------------
% clf ; clear ; hold on  ; 
% x = linspace(0,2*pi,50)'; y = sin(x); dy = .1*(1+rand(size(y))).*y;  % made-up error values
% hf = fill([x;flipud(x)],[y-dy;flipud(y+dy)],[1 1 1]*0.8,'linestyle','none');
% plot(x,y,'r-','linewidth',2)


%% -------------------------------- line with varying color -----------------------------
% clf
% xlim([1,5]);
% h2b = plot(2:6,  15:-1:11, '.-r', 'LineWidth',8, 'DisplayName',' 0.7'); h2b.Color(4)=0.7;  % 30% transparent
% legend('show','Location','west') ; 
% cd=uint8([255,200,250,50,0; 0,50,250,150,200; 0,0,0,100,150; 179,150,200,70,50]) ; 
% set(h2b.Edge, 'ColorBinding','interpolated', 'ColorData',cd) ; 

%% ------------------------------ export figures --------------------------------
% output fig1.eps with color, 300 dpi and tiff preview
%print -depsc -r300 -tiff -r300 fig2
% print('-depsc','-tiff','-r600','fig3')


%% ------------------------------ export gif animation ---------------------------
% load mri
% imwrite(D,map,'gif_example.gif','DelayTime',0.2,'LoopCount',inf) 

%% ---------------------------- export movie --------------------------------
%         nz = 24 ; 
%         file_DCE_video = fullfile(pars.dir.Output, sprintf('%s_DCE_video_nz=%d.avi', pars.scan.studyID, nz)); 
%         v = VideoWriter(file_DCE_video) ; open(v) ; 
%         figure(1) ; clf ; colormap gray ; set(gcf,'color','k') ; set(gcf, 'Units', 'inches', 'Position', [9 6 4.5 4])
%         for nt = 1:size(DCE_nii.img,4)
%             imagesc(fliplr(rot90(DCE_nii.img(:,:,nz,nt),-1))) ; axis image ; axis off ; caxis([0 2000])
%             frame = getframe(gcf) ; 
%             writeVideo(v,frame)
%         end
%         close(v)
