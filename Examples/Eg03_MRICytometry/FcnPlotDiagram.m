function FcnPlotDiagram(sim, vdist, fitopts)
    %% plot projection 
%     fcolor_true = [0.89 0.94 0.90]*0.96 ; ecolor_true = [0.3 0.5 1] ; line_true = 'r--' ;  lw_true = 2 ; 
    col = [0.07 0.6 1]; %[0 0 1] ; 
    color_saturation = 0.6 ; color_new = col + (1-col)*(1-color_saturation) ; 
    fcolor_true = color_new ; ecolor_true = color_new ; line_true = 'r--' ;  lw_true = 2 ;  line_fit = 'k-' ; lw_fit = 2 ; 

    matrix_invw = vdist{11} ; matrix_ex = vdist{12} ; matrix_in = vdist{13} ; 
    
    %%
    figure ; clf ; mati.Tools.SetPlotPars() ; 
    set(gcf,'Units','inches', 'PaperUnits', 'inches', 'position',[6 6 6.9 2.4]) ; 
    fs.axis=6 ; fs.label=8 ; fs.legend=9 ; set(gcf,'DefaultAxesFontSize',fs.axis, 'DefaultAxesFontName', 'Times New Roman' )
    Nrow = 2 ; Ncol = 6 ; 
    xmin = 0.05 ; ymin = 0.14 ; lx = 0.09 ; ly = 0.3 ; dx = 0.06 ; dy = 0.12 ; [nx,ny]=meshgrid(1:Ncol,1:Nrow) ; 
    for n=1:length(nx(:))
        hax(n) = axes('Position',[xmin+(nx(n)-1)*(lx+dx), ymin+(Nrow-ny(n))*(ly+dy), lx, ly]) ; 
    end
    
    %% betaex
    nax=1 ; axes(hax(nax)) ; hold on ;  set(gca,FN,fn,FS,fs.axis) ;
%     hae1(4) = area(sim.betaexs, sim.vbetaex*sim.vex*100) ; set(hae1(4), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle',':',LW,lw_true) ; 
    hbar(4) = bar(sim.betaexs, sim.vbetaex*sim.vex*100) ; set(hbar(4), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle','-',LW,lw_true,'BarWidth',0.2) ; 
%     plot(sim.betaexs, sim.vbetaex*sim.vex*100,line_true,LW,lw_true) ; 
    plot(fitopts.betaexs, vdist{4}*100,line_fit,LW,lw_fit) ; 
    box on ; ylim([0 20]) ; xlim([0 10]) ; xlabel('\it\betaex\rm [\mum^2]',FN,fn, FS,fs.label) ; ylabel('probability [%]',FN,fn, FS,fs.label) ; set(get(gca,'XLabel'),'Position',[5 -5 -1])
    set(gca,'xtick',0:2:10,'xticklabel',0:2:10, 'ytick',0:5:20,'yticklabel',0:5:20, 'XAxisLocation','bottom','YAxisLocation', 'left') ; camroll(-90) ; hold off ; 

    %% 
    nax=2 ; delete(hax(nax)) ; 
    
    %% matrix_ex
    nax=3 ; axes(hax(nax)) ; 
    imagesc(matrix_ex) ; colormap hot ; set(gca,FN,fn,FS,fs.axis) ; 
    xlabel('\itDex\rm [\mum^2/ms]',FN,fn, FS,fs.label) ;   set(gca,'xtick',[1 5 10 15 20],'xticklabel',0:1:3, 'XAxisLocation', 'top') ; 
    ylabel('\it\betaex\rm [\mum^2]',FN,fn, FS,fs.label) ;   set(gca,'ytick',[1:4:22],'yticklabel',0:2:10, 'YAxisLocation', 'right') ; set(get(gca,'YLabel'),'Rotation',270,'Position',[-2 12 1]) ; 

    %% Dex
    nax=4 ; axes(hax(nax)) ; hold on ;  set(gca,FN,fn,FS,fs.axis) ; 
%     hae1(7) = area(sim.Dexs, sim.vDex*sim.vex*100) ; set(hae1(7), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle',':',LW,lw_true) ; 
%     plot(sim.Dexs, sim.vDex*sim.vex*100,line_true,LW,lw_true) ; 
%     plot([1 1]*pars.structure.Dfree ,[0 sim.vfree]*100,line_true,LW,lw_true) ; 
    hbar(7) = bar(sim.Dexs, sim.vDex*sim.vex*100) ; set(hbar(7), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle','-',LW,lw_true,'BarWidth',1) ; 
%     hbar(8) = bar(pars.structure.Dfree ,sim.vfree*100) ;  set(hbar(8), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle','-',LW,lw_true,'BarWidth',0.1) ; 
    plot(fitopts.Dexs, vdist{3}*100,line_fit,LW,lw_fit) ; 
    box on ; ylim([0 10]) ; xlim([0 3.3]) ; xlabel('\itDex\rm [\mum^2/ms]',FN,fn, FS,fs.label) ; ylabel('probability [%]',FN,fn, FS,fs.label) ; %set(get(gca,'YLabel'),'Position',[4.3 5.1 -1])
    set(gca,'xtick',0:1:3,'xticklabel',0:1:3, 'ytick',0:2:10,'yticklabel',0:2:10,'YAxisLocation','left')
    
    %% matrix_invw
    nax=5 ; axes(hax(nax)) ; 
    imagesc(matrix_invw) ; colormap hot ; set(gca,FN,fn,FS,fs.axis ) ; 
    ylabel('\itDin\rm [\mum^2/ms]',FN,fn, FS,fs.label) ;  set(gca,'ytick',[1 5 10 15 20],'yticklabel',0:1:3, 'YAxisLocation', 'left') ; %set(get(gca,'YLabel'),'Rotation',270,'Position',[64 9 1]) ; 
    xlabel('\itd\rm [\mum]',FN,fn, FS,fs.label) ;   set(gca,'xtick',[1 10:10:50],'xticklabel',0:5:25, 'XAxisLocation', 'top') ; 

    %% dvw
    nax=6 ; axes(hax(nax)) ;  hold on ;  
%     hae1(6) = area(sim.ds, sim.vd*sim.vin*100) ; set(hae1(6), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle',':',LW,lw_true) ; 
%     plot(sim.ds, sim.vd*sim.vin*100, line_true,LW,lw_true)
    hbar(6) = bar(sim.ds, sim.vdvw*sim.vin*100) ; set(hbar(6), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle','-',LW,lw_true,'BarWidth',0.1) ; 
    plot(fitopts.ds, vdist{1}*100,line_fit,LW,lw_fit)
%     plot(fitopts.ds, vdist{5}*100,'b-',LW,lw_fit)
     box on ; set(gca,FN,fn,FS,fs.axis, 'ylim',[0 10],'xlim',[-1 25]) ;  xlabel('\itdvw\rm [\mum]',FN,fn, FS,fs.label) ; ylabel('probability [%]',FN,fn, FS,fs.label) ;  set(gca,'xtick',0:5:25,'xticklabel',0:5:25, 'ytick',0:2:10,'yticklabel',0:2:10)

     %% Din
    nax=7 ; axes(hax(nax)) ; hold on ; 
    if sim.flag.DinDist  == 'y'
%         hae1(1) = area(sim.Dins, sim.vDin*sim.vin*100) ; set(hae1(1), 'facecolor',fcolor_true,'edgecolor',fcolor_true) ; 
%         plot(sim.Dins, sim.vDin*sim.vin*100, line_true,LW,lw_true)
        hbar(1) = bar(sim.Dins, sim.vDin*sim.vin*100) ; set(hbar(1), 'facecolor',fcolor_true,'edgecolor',fcolor_true,'BarWidth',1) ; 
    else
        plot(sim.Din*[1 1], [0 1]*sim.vin*100, line_true,LW,lw_true)
    end
    plot(fitopts.Dins, vdist{2}*100,line_fit,LW,lw_fit)
    box on ; set(gca,FN,fn,FS,fs.axis, 'ylim',[0 10],'xlim',[0 3.3]) ; xlabel('{\itDin} [\mum^2/ms]',FN,fn, FS,fs.label) ;   ylabel('probability [%] ',FN,fn, FS,fs.label) ; set(gca,'xtick',0:1:3,'xticklabel',0:1:3, 'ytick',0:2:10,'yticklabel',0:2:10,'Xdir','reverse',  'XAxisLocation', 'bottom', 'YAxisLocation', 'right') ; set(get(gca,'XLabel'),'Rotation',270,'Position',[1.5 11 1],FS,fs.label) ; camroll(90) ; hold off ; 

    %% 
    nax=8 ; delete(hax(nax)) ; 

    %% matrix_in
    nax=9 ; axes(hax(nax)) ; 
    imagesc(matrix_in) ; colormap hot ; set(gca,FN,fn,FS,fs.axis) ; 
    ylabel('\itDin\rm [\mum^2/ms]',FN,fn, FS,fs.label) ;  set(gca,'ytick',[1 5 10 15 20],'yticklabel',0:1:3, 'YAxisLocation', 'left') ; %set(get(gca,'YLabel'),'Rotation',270,'Position',[64 9 1]) ; 
    xlabel('\itd\rm [\mum]',FN,fn, FS,fs.label) ;   set(gca,'xtick',[1 10:10:50],'xticklabel',0:5:25, 'XAxisLocation', 'top') ; 

    %% d
    nax=10 ; axes(hax(nax)) ;  hold on ;  
%     hae1(6) = area(sim.ds, sim.vd*sim.vin*100) ; set(hae1(6), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle',':',LW,lw_true) ; 
%     plot(sim.ds, sim.vd*sim.vin*100, line_true,LW,lw_true)
    hbar(6) = bar(sim.ds, sim.vd*sim.vin*100) ; set(hbar(6), 'facecolor',fcolor_true, 'edgecolor',ecolor_true,'linestyle','-',LW,lw_true,'BarWidth',0.1) ; 
    plot(fitopts.ds, vdist{5}*100,line_fit,LW,lw_fit)
    box on ; set(gca,FN,fn,FS,fs.axis, 'ylim',[0 10],'xlim',[-1 25]) ; xlabel('\itd\rm [\mum]',FN,fn,FS,fs.label) ; ylabel('probability [%]',FN,fn,FS,fs.label) ;  set(gca,'xtick',0:5:25,'xticklabel',0:5:25, 'ytick',0:2:10,'yticklabel',0:2:10)
    hlg = legend('truth','MRI cytometry')  ; set(hlg, FS, fs.legend,'location','EastOutside','color','none','edgecolor','none') ; hlg.Position = hlg.Position + [lx*1.2 0 0 0] ; 

     %% Din
    nax=11 ; axes(hax(nax)) ; hold on ; 
    if sim.flag.DinDist  == 'y'
%         hae1(1) = area(sim.Dins, sim.vDin*sim.vin*100) ; set(hae1(1), 'facecolor',fcolor_true,'edgecolor',fcolor_true) ; 
%         plot(sim.Dins, sim.vDin*sim.vin*100, line_true,LW,lw_true)
        hbar(1) = bar(sim.Dins, sim.vDin*sim.vin*100) ; set(hbar(1), 'facecolor',fcolor_true,'edgecolor',fcolor_true,'BarWidth',1) ; 
    else
        plot(sim.Din*[1 1], [0 1]*sim.vin*100, line_true,LW,lw_true)
    end
    plot(fitopts.Dins, vdist{6}*100,line_fit,LW,lw_fit)
    box on ; set(gca,FN,fn,FS,fs.axis, 'ylim',[0 10],'xlim',[0 3.3]) ; xlabel('{\itDin} [\mum^2/ms]',FN,fn,FS,fs.label) ;   ylabel('probability [%] ',FN,fn,FS,fs.label) ; set(gca,'xtick',0:1:3,'xticklabel',0:1:3, 'ytick',0:2:10,'yticklabel',0:2:10,'Xdir','reverse',  'XAxisLocation', 'bottom', 'YAxisLocation', 'right') ; set(get(gca,'XLabel'),'Rotation',270,'Position',[1.5 -4 1]) ; camroll(90) ; hold off ; 

    %% 
    nax=12 ; delete(hax(nax)) ; 

end
