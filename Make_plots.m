%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT THE RESULTS

load(Filename_parameters);
load(Filename_results);

if(Single_Point_source == 0)
if(Sampling_type == 0 || Sampling_type == 1)
    
    if(To_save(1))
    if(Coronagraph == 0)    
    figure('name', 'Plane B')
    plot(vp*s_b, log10(I_B), 'b-', 'Linewidth', 1.5)
    xlabel('Radius in plane B ($R_{\odot}$)','Interpreter', 'LaTex');
    ylabel('Intensity $I_{B}(r)$ (log)','Interpreter', 'LaTex');
    title('Plane B', 'Interpreter', 'latex');
    xlim([0 3]);
    ylim([-6 1]);
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end

    if(Coronagraph == 1)
    figure('name', 'Plane O''')
    plot(vp*s_op/(z1*tan(Rsun)), log10(I_Op), 'r-', 'Linewidth', 1.5)
    hold on
    plot([Rio_ext, Rio_ext]/(z1*tan(Rsun)), [-10, 1], 'k-')
    xlabel('Radius in plane O'' ($R_{\odot}$)','Interpreter', 'LaTex');
    ylabel('Intensity $I_{O''}(r)$ (log)','Interpreter', 'LaTex');
    title('Plane O''', 'Interpreter', 'latex');
    xlim([0 3]);
    ylim([-10 0]);
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end
    end
    
    if(To_save(2))
    figure('name', 'Plane C')
    plot(vp*s_a*1000, log10(I_C), 'r-', 'Linewidth', 1.5)
    hold on
    plot([Rp*Fls*1000, Rp*Fls*1000], [-10, 1], 'k-')
    xlabel('Radius in plane C (mm)','Interpreter', 'LaTex');
    ylabel('Intensity $I_{C}(r)$ (log)','Interpreter', 'LaTex');
    title('Plane C', 'Interpreter', 'latex');
    xlim([0 30]);
    ylim([-10 0]);
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end

    if(To_save(3))
    figure('name', 'Plane D')
    plot(vp*s_d, log10(I_D), 'r-', 'Linewidth', 1.5)
    xlabel('Radius in plane D ($R_{\odot}$)','Interpreter', 'LaTex');
    ylabel('Intensity $I_{D}(r)$ (log)','Interpreter', 'LaTex');
    title('Plane D', 'Interpreter', 'latex');
    xlim([0 3]);
    ylim([-12 0]);
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end

elseif( Sampling_type == 2)
    
    if(To_save(1))
    if(Coronagraph == 0)    
    figure('name', 'Plane B')
    imagesc(v*s_b, v*s_b, log10(I_B))
    xlabel('x ($R_{\odot}$)','Interpreter', 'LaTex');
    ylabel('y ($R_{\odot}$)','Interpreter', 'LaTex');
    title('Plane B', 'Interpreter', 'latex');
    colorbar
    colorAxis = colorbar;
    colorAxis.Label.String='Intensity $I_{B}(x,y)$ (log))'; 
    colorAxis.Label.Interpreter='latex';
    colorAxis.Label.FontSize=12;
    colorAxis.TickLabelInterpreter='latex'; 
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end

    if(Coronagraph == 1)
    figure('name', 'Plane O''')
    imagesc(v*s_op*1000, v*s_op*1000, log10(I_Op))
    xlabel('x (mm)','Interpreter', 'LaTex');
    ylabel('y (mm)','Interpreter', 'LaTex');
    title('Plane O''', 'Interpreter', 'latex');
    colorbar
    colorAxis = colorbar;
    colorAxis.Label.String='Intensity $I_{O''}(x,y)$ (log))'; 
    colorAxis.Label.Interpreter='latex';
    colorAxis.Label.FontSize=12;
    colorAxis.TickLabelInterpreter='latex'; 
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end
    end
    
    if(To_save(2))
    figure('name', 'Plane C')
    imagesc(v*s_a*1000, v*s_a*1000, log10(I_C))
    xlabel('x (mm)','Interpreter', 'LaTex');
    ylabel('y (mm)','Interpreter', 'LaTex');
    title('Plane C', 'Interpreter', 'latex');
    colorbar
    colorAxis = colorbar;
    colorAxis.Label.String='Intensity $I_{C}(x,y)$ (log))'; 
    colorAxis.Label.Interpreter='latex';
    colorAxis.Label.FontSize=12;
    colorAxis.TickLabelInterpreter='latex'; 
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end

    if(To_save(3))
    figure('name', 'Plane D')
    imagesc(v*s_b, v*s_b, log10(I_D))
    xlabel('x ($R_{\odot}$)','Interpreter', 'LaTex');
    ylabel('y ($R_{\odot}$)','Interpreter', 'LaTex');
    title('Plane D', 'Interpreter', 'latex');
    colorbar
    colorAxis = colorbar;
    colorAxis.Label.String='Intensity $I_{D}(x,y)$ (log))'; 
    colorAxis.Label.Interpreter='latex';
    colorAxis.Label.FontSize=12;
    colorAxis.TickLabelInterpreter='latex'; 
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end
end

elseif( Single_Point_source == 1)
    
    if(To_save(1))
    if(Coronagraph == 0)    
    figure('name', 'Plane B')
    imagesc(v*s_b, v*s_b, log10(I_B))
    xlabel('x ($R_{\odot}$)','Interpreter', 'LaTex');
    ylabel('y ($R_{\odot}$)','Interpreter', 'LaTex');
    title('Point Source in Plane B', 'Interpreter', 'latex');
    colorbar
    colorAxis = colorbar;
    colorAxis.Label.String='Intensity $I_{B}(x,y)$ (log))'; 
    colorAxis.Label.Interpreter='latex';
    colorAxis.Label.FontSize=12;
    colorAxis.TickLabelInterpreter='latex'; 
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end

    if(Coronagraph == 1)
    figure('name', 'Plane O''')
    imagesc(v*s_op*1000, v*s_op*1000, log10(I_Op))
    xlabel('x (mm)','Interpreter', 'LaTex');
    ylabel('y (mm)','Interpreter', 'LaTex');
    title('Point Source in Plane O''', 'Interpreter', 'latex');
    colorbar
    colorAxis = colorbar;
    colorAxis.Label.String='Intensity $I_{O''}(x,y)$ (log))'; 
    colorAxis.Label.Interpreter='latex';
    colorAxis.Label.FontSize=12;
    colorAxis.TickLabelInterpreter='latex'; 
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end
    end

    if(To_save(2))
    figure('name', 'Plane C')
    imagesc(v*s_a*1000, v*s_a*1000, log10(I_C))
    xlabel('x (mm)','Interpreter', 'LaTex');
    ylabel('y (mm)','Interpreter', 'LaTex');
    title('Point Source in Plane C', 'Interpreter', 'latex');
    colorbar
    colorAxis = colorbar;
    colorAxis.Label.String='Intensity $I_{C}(x,y)$ (log))'; 
    colorAxis.Label.Interpreter='latex';
    colorAxis.Label.FontSize=12;
    colorAxis.TickLabelInterpreter='latex'; 
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end

    if(To_save(3))
    figure('name', 'Plane D')
    imagesc(v*s_b, v*s_b, log10(I_D))
    xlabel('x ($R_{\odot}$)','Interpreter', 'LaTex');
    ylabel('y ($R_{\odot}$)','Interpreter', 'LaTex');
    title('Point Source in Plane D', 'Interpreter', 'latex');
    colorbar
    colorAxis = colorbar;
    colorAxis.Label.String='Intensity $I_{D}(x,y)$ (log))'; 
    colorAxis.Label.Interpreter='latex';
    colorAxis.Label.FontSize=12;
    colorAxis.TickLabelInterpreter='latex'; 
    % l=legend();
    % l.Location = 'eastoutside';
    % l.Box='off';
    % l.Interpreter='LaTex';
    Fig=gca;
    set(Fig,'TickLabelInterpreter', 'LaTex');
    Fig.XMinorTick='on';
    Fig.YMinorTick='on';
    Fig.TickLength=[0.03 0.02];
    Fig.FontSize = 12;
    end

end