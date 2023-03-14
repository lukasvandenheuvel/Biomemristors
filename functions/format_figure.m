function [fig,fig_property] = format_figure(width,height,fontsize,format)
    clear figure_property;
    set(gcf,'color','w');
    fig_property.units = 'inches';
    fig_property.format = format; 
    fig_property.Preview= 'none';
    fig_property.Width= num2str(width); 
    fig_property.Height= num2str(height); 
    fig_property.Units= 'inches';
    fig_property.Color= 'rgb';
    fig_property.Background= 'w';
    fig_property.FixedfontSize= num2str(fontsize);
    fig_property.ScaledfontSize= 'auto';
    fig_property.FontMode= 'scaled';
    fig_property.FontSizeMin= num2str(fontsize);
    fig_property.FixedLineWidth= '1';
    fig_property.ScaledLineWidth= 'auto';
    fig_property.LineMode= 'none';
    fig_property.LineWidthMin= '0.1';
    fig_property.FontName= 'Arial';
    fig_property.FontWeight= 'auto';
    fig_property.FontAngle= 'auto';
    fig_property.FontEncoding= 'latin1';
    fig_property.PSLevel= '3';
    fig_property.Renderer= 'painters';
    fig_property.Resolution= '600';
    fig_property.LineStyleMap= 'none';
    fig_property.ApplyStyle= '0';
    fig_property.LockAxes= 'off';
    fig_property.LockAxesTicks= 'off';
    fig_property.ShowUI= 'off';
    fig_property.SeparateText= 'off';
    fig=gcf;
    set(fig,'PaperUnits','inches');
    set(fig,'PaperPositionMode','auto');
    set(fig,'PaperSize',[str2double(fig_property.Width) str2double(fig_property.Height)]); % Canvas Size
    set(fig,'Units','inches');
    set(gca,'FontSize',fontsize);
end