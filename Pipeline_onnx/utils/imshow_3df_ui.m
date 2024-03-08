function [c] = imshow_3df_ui(v)

    set(gcf, 'InnerPosition', [300 300 860 618])
    clf            
    var = v.a;

    c.image     = axes('position',[.00 .25 .60 .75]);
    if var.sno > 1                
        c.slider    = uicontrol('Style','slider','Units','normalized','Position', [0 0.12 0.6 0.04],'Min',1,'Max',var.sno,'Value',var.s,'SliderStep',[1/(var.sno-1) 10/(var.sno-1)]);
        c.slidertxt = uicontrol('Style','text'  ,'Units','normalized','Position', [0 0.16 0.6 0.04],'String',sprintf('image size = (%d, %d),  Slice# %d / %d ', size(var.IMGLv1,1), size(var.IMGLv1,2), var.s, var.sno), 'FontSize', 9);
    else
        c.slider    = uicontrol('Style','slider','Units','normalized','Position', [0 0.12 0.6 0.04],'Min',1,'Max',10,'Value',5,'SliderStep',[1/9 10/9]);
        c.slider.Visible = 'off';
        c.slidertxt = uicontrol('Style','text'  ,'Units','normalized','Position', [0 0.16 0.6 0.04],'String',sprintf('image size = (%d, %d),  2D image ',size(var.IMGLv1,1), size(var.IMGLv1,2)), 'FontSize', 9);
    end             

    c.panel1.name = uipanel('Title','info',         'FontSize',10,'Position',[.65 .70 .15 .30]);
    c.panel2.name = uipanel('Title','image control','FontSize',10,'Position',[.65 .05 .15 .65]);
    c.panel3.name = uipanel('Title','view',         'FontSize',10,'Position',[.82 .70 .15 .30]);
    c.panel4.name = uipanel('Title','image range',  'FontSize',10,'Position',[.82 .05 .15 .65]);

    % panel 1
    arrowcontents = sprintf([strjust(sprintf('%12s',var.arrow(1)), 'center'),'\n',...
                             strjust(sprintf('%12s',' ^'), 'center'),'\n',...
                             strjust(sprintf('%10s',[var.arrow(3),'  <-    ->  ',var.arrow(4)]), 'center'),'\n',...
                             strjust(sprintf('%12s','v'), 'center'),'\n',...
                             strjust(sprintf('%12s',var.arrow(2)), 'center')]);  
    c.panel1.arrow = uicontrol('Style','text','Units','normalized','Position',[.68 .84 .09 .12],'String',arrowcontents);
    c.panel1.imdim = uicontrol('Style','text','Units','normalized','Position',[.68 .77 .09 .04],'String',sprintf('input = %dD',v.IMGdim));
    if v.IMGdim == 2
        imsizcontents = sprintf('(%d,%d)', size(var.IMGLv0,1), size(var.IMGLv0,2));
    elseif v.IMGdim == 3
        imsizcontents = sprintf('(%d,%d,%d)', size(var.IMGLv0,1), size(var.IMGLv0,2), size(var.IMGLv0,3));
    elseif v.IMGdim == 4
        imsizcontents = sprintf('(%d,%d,%d,%d)', size(var.IMGLv0,1), size(var.IMGLv0,2), size(var.IMGLv0,3), size(var.IMGLv0,4));
    elseif v.IMGdim == 5
        imsizcontents = sprintf('(%d,%d,%d,%d,%d)', size(var.IMGLv0,1), size(var.IMGLv0,2), size(var.IMGLv0,3), size(var.IMGLv0,4), size(var.IMGLv0,5));
    end            
    c.panel1.imsiz = uicontrol('Style','text','Units','normalized','Position',[.66 .73 .13 .04],'String',imsizcontents);

    % panel 2            
    c.panel2.orign  = uicontrol('Style','pushbutton','Units','normalized','Position',[.66 .60 .13 .04],'String','original');
    c.panel2.trans  = uicontrol('Style','pushbutton','Units','normalized','Position',[.66 .54 .13 .04],'String','transpose');
    c.panel2.flipud = uicontrol('Style','pushbutton','Units','normalized','Position',[.66 .48 .13 .04],'String','flipud');
    c.panel2.fliplr = uicontrol('Style','pushbutton','Units','normalized','Position',[.66 .42 .13 .04],'String','fliplr');
    c.panel2.rotl   = uicontrol('Style','pushbutton','Units','normalized','Position',[.66 .36 .06 .04],'String','L');
    c.panel2.rotr   = uicontrol('Style','pushbutton','Units','normalized','Position',[.73 .36 .06 .04],'String','R');
    c.panel2.flat   = uicontrol('Style','pushbutton','Units','normalized','Position',[.66 .30 .13 .04],'String','flattenaxes');
    c.panel2.cropbx = uicontrol('Style','checkbox'  ,'Units','normalized','Position',[.66 .24 .13 .04],'String','crop');
    c.panel2.cropln = uicontrol('Style','checkbox'  ,'Units','normalized','Position',[.66 .20 .13 .04],'String','show line');
    c.panel2.cropxt = uicontrol('Style','text'      ,'Units','normalized','Position',[.66 .13 .13 .04],'String','x','HorizontalAlignment','left');
    c.panel2.cropyt = uicontrol('Style','text'      ,'Units','normalized','Position',[.66 .08 .13 .04],'String','y','HorizontalAlignment','left');
    c.panel2.cropx1 = uicontrol('Style','edit'      ,'Units','normalized','Position',[.68 .13 .05 .04],'String',sprintf('%3d',var.CROP(1)));
    c.panel2.cropx2 = uicontrol('Style','edit'      ,'Units','normalized','Position',[.74 .13 .05 .04],'String',sprintf('%3d',var.CROP(2)));
    c.panel2.cropy1 = uicontrol('Style','edit'      ,'Units','normalized','Position',[.68 .08 .05 .04],'String',sprintf('%3d',var.CROP(3)));
    c.panel2.cropy2 = uicontrol('Style','edit'      ,'Units','normalized','Position',[.74 .08 .05 .04],'String',sprintf('%3d',var.CROP(4)));

    % panel 3                 
    c.panel3.ax = uicontrol('Style','togglebutton','Units','normalized','Position',[.83 .93 .13 .04],'String','Axial', 'Value', true);
    c.panel3.sg = uicontrol('Style','togglebutton','Units','normalized','Position',[.83 .87 .13 .04],'String','Sagittal');
    c.panel3.cr = uicontrol('Style','togglebutton','Units','normalized','Position',[.83 .81 .13 .04],'String','Coronal');
    c.panel3.d4 = uicontrol('Style','checkbox'    ,'Units','normalized','Position',[.83 .75 .13 .04],'String','4th dim');
    c.panel3.d5 = uicontrol('Style','checkbox'    ,'Units','normalized','Position',[.83 .71 .13 .04],'String','5th dim');            

    % panel 4
    imgnumcontents = {};
    for i=1:v.IMGnum
        imgnumcontents{i} = num2str(i);
    end
    c.panel4.numtxt  = uicontrol('Style','text'        ,'Units','normalized','Position', [.83 .61 .13 .04],'String',sprintf('# image = %d',v.IMGnum));
    c.panel4.imgnum  = uicontrol('Style','listbox'     ,'Units','normalized','Position', [.83 .55 .13 .06],'String',imgnumcontents);
    c.panel4.mintxt  = uicontrol('Style','text'        ,'Units','normalized','Position', [.83 .49 .06 .04],'String','min');
    c.panel4.maxtxt  = uicontrol('Style','text'        ,'Units','normalized','Position', [.90 .49 .06 .04],'String','max');
    c.panel4.minval  = uicontrol('Style','edit'        ,'Units','normalized','Position', [.83 .46 .06 .04],'String',sprintf('%3.2f',v.disprange(1,1)));
    c.panel4.maxval  = uicontrol('Style','edit'        ,'Units','normalized','Position', [.90 .46 .06 .04],'String',sprintf('%3.2f',v.disprange(2,1)));
    c.panel4.auto    = uicontrol('Style','pushbutton'  ,'Units','normalized','Position', [.83 .39 .13 .04],'String','Auto');
    c.panel4.cmaptxt = uicontrol('Style','text'        ,'Units','normalized','Position', [.83 .32 .13 .04],'String','colormap');
    c.panel4.cmappop = uicontrol('Style','popupmenu'   ,'Units','normalized','Position', [.83 .278 .13 .04],'String',{'gray','jet'});
%             c.panel4.datatxt = uicontrol('Style','text'        ,'Units','normalized','Position', [.83 .20 .13 .04],'String','Datatip');
%             c.panel4.datatip = uicontrol('Style','togglebutton','Units','normalized','Position', [.83 .17 .13 .04],'String','on/off');

    fn1 = fieldnames(c);
    for i1=4:length(fn1)
        fn2 = fieldnames(c.(fn1{i1}));
        for i2=2:length(fn2)
            if strcmp(c.(fn1{i1}).(fn2{i2}).Style,'edit') || strcmp(c.(fn1{i1}).(fn2{i2}).Style,'text')
                set(c.(fn1{i1}).(fn2{i2}),'FontSize',9)
            else
                set(c.(fn1{i1}).(fn2{i2}),'FontSize',10)
            end
        end
    end
end

