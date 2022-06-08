function  imshow_3df( varargin )
    % input type 1 : (IMG, IMG, IMG, ...)
    % input type 2 : (IMG, disprange, IMG, IMG, disprange, ...)
    % input type 3 : (..., 'slice', 1:5:120)
    % input type 4 : (..., 'range', [0 1])
    
    for i=1:length(varargin)
        if ischar(varargin{i})
            break
        end        
        if ~ isequal(size(varargin{i}), [1 2])
            varargin{i} = permute(varargin{i},[2 1 3 4 5]);
        end               
    end
    
    v = imshow_3df_parse_input(varargin);
    c = imshow_3df_ui(v);
    IMGLv3 = 0;
    IMGLv3idx = 0;
    IMGLv1process
    IMGLv2process
    IMGView

    function [var] = load_var()
        if v.IMGview == 'A'
            var = v.a;
        elseif v.IMGview == 'S'
            var = v.s;
        elseif v.IMGview == 'C'
            var = v.c;
        end   

        if v.flag.view4thdim
            var.sno = v.sno4;
            var.s   = v.s4;
        elseif v.flag.view5thdim
            var.sno = v.sno5;
            var.s   = v.s5;
        end            
    end
    function save_var(var)
        if v.flag.view4thdim
            v.sno4 = var.sno;
            v.s4   = var.s;
        elseif v.flag.view5thdim
            v.sno5 = var.sno;
            v.s5   = var.s;
        end
        if v.IMGview == 'A'
            if v.flag.view4thdim || v.flag.view5thdim
                var.sno = v.a.sno;
                var.s   = v.a.s;
            end
            v.a   = var;
        elseif v.IMGview == 'S'
            if v.flag.view4thdim || v.flag.view5thdim
                var.sno = v.a.sno;
                var.s   = v.a.s;
            end          
            v.s = var;
        elseif v.IMGview == 'C'
            if v.flag.view4thdim || v.flag.view5thdim
                var.sno = v.a.sno;
                var.s   = v.a.s;
            end                    
            v.c = var;
        end            
    end 

    function IMGLv1process(object)
        var = load_var();
        set(c.panel2.cropx1,'String',sprintf('%3d',var.CROP(1)));
        set(c.panel2.cropx2,'String',sprintf('%3d',var.CROP(2)));
        set(c.panel2.cropy1,'String',sprintf('%3d',var.CROP(3)));
        set(c.panel2.cropy2,'String',sprintf('%3d',var.CROP(4)));

        arrowcontents = sprintf([strjust(sprintf('%12s',var.arrow(1)), 'center'),'\n',...
                                 strjust(sprintf('%12s',' ^'), 'center'),'\n',...
                                 strjust(sprintf('%10s',[var.arrow(3),'  <-    ->  ',var.arrow(4)]), 'center'),'\n',...
                                 strjust(sprintf('%12s','v'), 'center'),'\n',...
                                 strjust(sprintf('%12s',var.arrow(2)), 'center')]);   
        set(c.panel1.arrow, 'String',arrowcontents);
    end
    function IMGLv2process(object)
        var = load_var();
        var.IMGLv2 = var.IMGLv1;
        if get(c.panel2.cropln, 'Value')
            var.IMGLv2(var.CROP(1),:,:,:,:,:) = repmat(cos(pi*(1:size(var.IMGLv1,2)) )*inf,[1,1,size(var.IMGLv1,3),size(var.IMGLv1,4),size(var.IMGLv1,5),size(var.IMGLv1,6)]);
            var.IMGLv2(var.CROP(2),:,:,:,:,:) = repmat(cos(pi*(1:size(var.IMGLv1,2)) )*inf,[1,1,size(var.IMGLv1,3),size(var.IMGLv1,4),size(var.IMGLv1,5),size(var.IMGLv1,6)]);
            var.IMGLv2(:,var.CROP(3),:,:,:,:) = repmat(cos(pi*(1:size(var.IMGLv1,1))')*inf,[1,1,size(var.IMGLv1,3),size(var.IMGLv1,4),size(var.IMGLv1,5),size(var.IMGLv1,6)]);
            var.IMGLv2(:,var.CROP(4),:,:,:,:) = repmat(cos(pi*(1:size(var.IMGLv1,1))')*inf,[1,1,size(var.IMGLv1,3),size(var.IMGLv1,4),size(var.IMGLv1,5),size(var.IMGLv1,6)]);
        end
        if get(c.panel2.cropbx, 'Value')
            var.IMGLv2 = var.IMGLv2(var.CROP(1):var.CROP(2),var.CROP(3):var.CROP(4),:,:,:,:);
        end            
        if v.flag.view4thdim
            if v.IMGview == 'A'
                var.IMGLv2 = permute(var.IMGLv2(:,:,v.a.s,:,:,:),[1 2 4 5 3 6]);    
            elseif v.IMGview == 'S'
                var.IMGLv2 = permute(var.IMGLv2(:,:,v.s.s,:,:,:),[1 2 4 5 3 6]);    
            elseif v.IMGview == 'C'
                var.IMGLv2 = permute(var.IMGLv2(:,:,v.c.s,:,:,:),[1 2 4 5 3 6]);    
            end
        elseif v.flag.view5thdim
            if v.IMGview == 'A'
                var.IMGLv2 = permute(var.IMGLv2(:,:,v.a.s,:,:,:),[1 2 5 4 3 6]);    
            elseif v.IMGview == 'S'
                var.IMGLv2 = permute(var.IMGLv2(:,:,v.s.s,:,:,:),[1 2 5 4 3 6]);    
            elseif v.IMGview == 'C'
                var.IMGLv2 = permute(var.IMGLv2(:,:,v.c.s,:,:,:),[1 2 5 4 3 6]);    
            end            
        end
        save_var(var);
    end
    function IMGView(object)
        var = load_var();
        IMGLv3 = var.IMGLv2;
        IMGLv3idx = zeros(size(IMGLv3));
        set(c.panel4.minval, 'Value', double(v.disprange(1,get(c.panel4.imgnum,'Value'))));
        set(c.panel4.maxval, 'Value', double(v.disprange(2,get(c.panel4.imgnum,'Value'))));
        for i=1:v.IMGnum
            IMGLv3(:,:,:,:,:,i) = (IMGLv3(:,:,:,:,:,i)-v.disprange(1,i))/(v.disprange(2,i)-v.disprange(1,i));
            IMGLv3idx(:,:,:,:,:,i) = i;
        end
        if v.IMGdim < 5
            IMGLv3 = permute(IMGLv3, [1 2 3 4 6 5]);
            IMGLv3 = IMGLv3(:,:,:,:,:,1);
            IMGLv3idx = permute(IMGLv3idx, [1 2 3 4 6 5]);
            IMGLv3idx = IMGLv3idx(:,:,:,:,:,1);
        end
        N = [size(IMGLv3,1),size(IMGLv3,2),size(IMGLv3,3),size(IMGLv3,4),size(IMGLv3,5)];
        if xor(v.flag.flatten, v.IMGnum > 1)
            IMGLv3 = permute(reshape(permute(IMGLv3,[1 4 3 2 5]),[N(1),N(4),N(3),N(2)*N(5)]),[1 4 3 2]);
            IMGLv3 = permute(reshape(permute(IMGLv3,[3 2 1 4]),[N(3),N(2)*N(5),N(1)*N(4)]),[3 2 1]);
            IMGLv3idx = permute(reshape(permute(IMGLv3idx,[1 4 3 2 5]),[N(1),N(4),N(3),N(2)*N(5)]),[1 4 3 2]);
            IMGLv3idx = permute(reshape(permute(IMGLv3idx,[3 2 1 4]),[N(3),N(2)*N(5),N(1)*N(4)]),[3 2 1]);            
        else
            IMGLv3 = permute(reshape(permute(IMGLv3,[4 2 3 1 5]),[N(4),N(2),N(3),N(1)*N(5)]),[4 2 3 1]);
            IMGLv3 = permute(reshape(permute(IMGLv3,[1 3 2 4]),[N(1)*N(5),N(3),N(2)*N(4)]),[1 3 2]);  
            IMGLv3idx = permute(reshape(permute(IMGLv3idx,[4 2 3 1 5]),[N(4),N(2),N(3),N(1)*N(5)]),[4 2 3 1]);
            IMGLv3idx = permute(reshape(permute(IMGLv3idx,[1 3 2 4]),[N(1)*N(5),N(3),N(2)*N(4)]),[1 3 2]);              
        end

        if var.sno > 1                
            set(c.slidertxt,'String',sprintf('image size = (%d, %d),  Slice# %d / %d ', size(var.IMGLv1,1), size(var.IMGLv1,2), var.s, var.sno));
        else
            set(c.slidertxt,'String',sprintf('image size = (%d, %d),  2D image ',size(var.IMGLv1,1), size(var.IMGLv1,2)));
        end                     
        cla(c.image);
        imshow(IMGLv3(:,:,var.s), [0 1])
        set(get(gca,'children'),'cdata',IMGLv3(:,:,var.s))
        colormap(gca, c.panel4.cmappop.String{get(c.panel4.cmappop,'Value')})
    end

    set(c.slider, 'Callback', @SliceSlider);   
    function SliceSlider(object,event)
        var = load_var();
        var.s = round(get(c.slider,'Value'));
        if var.sno > 1
            set(c.slider, 'Min',1,'Max',var.sno,'Value',var.s,'SliderStep',[1/(var.sno-1) 10/(var.sno-1)]);
            set(c.slidertxt,'String',sprintf('image size = (%d, %d),  Slice# %d / %d ', size(var.IMGLv1,1), size(var.IMGLv1,2), var.s, var.sno));
        else
            set(c.slidertxt,'String',sprintf('image size = (%d, %d),  2D image ',size(var.IMGLv1,1), size(var.IMGLv1,2)));
        end    
        set(get(gca,'children'),'cdata',IMGLv3(:,:,var.s))
        save_var(var);
    end    

    set(gcf,'WindowScrollWheelFcn', @mouseScroll);
    function mouseScroll(object,event)       
        var = load_var();
        UPDN = event.VerticalScrollCount;
        var.s = var.s - UPDN;
        if ( var.s < 1)
            var.s = 1;
        elseif ( var.s > var.sno )
            var.s = var.sno;
        end
        if var.sno > 1
            set(c.slider, 'Min',1,'Max',var.sno,'Value',var.s,'SliderStep',[1/(var.sno-1) 10/(var.sno-1)]);
            set(c.slidertxt,'String',sprintf('image size = (%d, %d),  Slice# %d / %d ', size(var.IMGLv1,1), size(var.IMGLv1,2), var.s, var.sno));
        else
            set(c.slidertxt,'String',sprintf('image size = (%d, %d),  2D image ',size(var.IMGLv1,1), size(var.IMGLv1,2)));
        end        
        set(get(gca,'children'),'cdata',IMGLv3(:,:,var.s))
        save_var(var);     
    end

    set(c.panel2.orign, 'Callback', @orignftn);
    function orignftn(object,event)
        var = load_var();
        var.IMGLv1 = var.IMGLv0;
        var.CROP = [1,size(var.IMGLv1,1),1,size(var.IMGLv1,2)];
        if v.IMGview == 'A'
            var.arrow = 'APRL';
        elseif v.IMGview == 'S'
            var.arrow = 'SIAP';
        elseif v.IMGview == 'C'
            var.arrow = 'SIRL';
        end
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView           
    end

    set(c.panel2.trans, 'Callback', @transftn);
    function transftn(object,event)
        var = load_var();
        var.IMGLv1 = permute(var.IMGLv1, [2 1 3 4 5 6]);
        var.CROP = [var.CROP(3),var.CROP(4),var.CROP(1),var.CROP(2)];
        var.arrow = [var.arrow(3),var.arrow(4),var.arrow(1),var.arrow(2)];
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView        
    end

    set(c.panel2.flipud, 'Callback', @flipudftn);
    function flipudftn(object,event)
        var = load_var();
        var.IMGLv1 = flip(var.IMGLv1,1);
        var.CROP = [size(var.IMGLv1,1)+1-var.CROP(2),size(var.IMGLv1,1)+1-var.CROP(1),var.CROP(3),var.CROP(4)];
        var.arrow = [var.arrow(2),var.arrow(1),var.arrow(3),var.arrow(4)];
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView        
    end

    set(c.panel2.fliplr, 'Callback', @fliplrftn);
    function fliplrftn(object,event)
        var = load_var();
        var.IMGLv1 = flip(var.IMGLv1,2);
        var.CROP = [var.CROP(1),var.CROP(2),size(var.IMGLv1,2)+1-var.CROP(4),size(var.IMGLv1,2)+1-var.CROP(3)];
        var.arrow = [var.arrow(1),var.arrow(2),var.arrow(4),var.arrow(3)];
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView        
    end

    set(c.panel2.rotl, 'Callback', @rotlftn);
    function rotlftn(object,event)
        var = load_var();
        var.IMGLv1 = rot90(var.IMGLv1,1);
        var.CROP = [size(var.IMGLv1,2)+1-var.CROP(4),size(var.IMGLv1,2)+1-var.CROP(3),var.CROP(1),var.CROP(2)];
        var.arrow = [var.arrow(4),var.arrow(3),var.arrow(1),var.arrow(2)];
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView        
    end

    set(c.panel2.rotr, 'Callback', @rotrftn);
    function rotrftn(object,event)
        var = load_var();
        var.IMGLv1 = rot90(var.IMGLv1,-1);
        var.CROP = [var.CROP(3),var.CROP(4),size(var.IMGLv1,1)+1-var.CROP(2),size(var.IMGLv1,1)+1-var.CROP(1)];
        var.arrow = [var.arrow(3),var.arrow(4),var.arrow(2),var.arrow(1)];
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView        
    end

    set(c.panel2.flat, 'Callback', @flatftn);
    function flatftn(object,event)
        v.flag.flatten = ~v.flag.flatten;
        IMGView        
    end

    set(c.panel2.cropbx, 'Callback', @cropbxftn);
    function cropbxftn(object,event)
        IMGLv2process
        IMGView              
    end

    set(c.panel2.cropln, 'Callback', @croplnftn);
    function croplnftn(object,event)
        IMGLv2process
        IMGView              
    end

    set(c.panel2.cropx1, 'Callback', @cropx1ftn);
    function cropx1ftn(object,event)
        var = load_var();
        x1 = str2num(get(c.panel2.cropx1, 'string'));
        if x1 < 1
            x1 = 1;
        elseif x1 >= var.CROP(2)
            x1 = var.CROP(2)-1;
        end
        var.CROP(1) = x1;
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView              
    end

    set(c.panel2.cropx2, 'Callback', @cropx2ftn);
    function cropx2ftn(object,event)
        var = load_var();
        x2 = str2num(get(c.panel2.cropx2, 'string'));
        if x2 <= var.CROP(1)
            x2 = var.CROP(1)+1;
        elseif x2 > size(var.IMGLv1,1)
            x2 = size(var.IMGLv1,1);
        end
        var.CROP(2) = x2;
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView              
    end

    set(c.panel2.cropy1, 'Callback', @cropy1ftn);
    function cropy1ftn(object,event)
        var = load_var();
        y1 = str2num(get(c.panel2.cropy1, 'string'));
        if y1 < 1
            y1 = 1;
        elseif y1 >= var.CROP(4)
            y1 = var.CROP(4)-1;
        end
        var.CROP(3) = y1;
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView              
    end

    set(c.panel2.cropy2, 'Callback', @cropy2ftn);
    function cropy2ftn(object,event)
        var = load_var();
        y2 = str2num(get(c.panel2.cropy2, 'string'));
        if y2 <= var.CROP(3)
            y2 = var.CROP(3)+1;
        elseif y2 > size(var.IMGLv1,2)
            y2 = size(var.IMGLv1,2);
        end
        var.CROP(4) = y2;
        save_var(var);
        IMGLv1process
        IMGLv2process
        IMGView              
    end

    set(c.panel3.ax, 'Callback', @axftn);
    function axftn(object,event)
        v.IMGview = 'A';
        set(c.panel3.ax, 'Value', true);
        set(c.panel3.sg, 'Value', false);
        set(c.panel3.cr, 'Value', false);
        IMGLv1process
        IMGLv2process            
        IMGView   
    end

    set(c.panel3.sg, 'Callback', @sgftn);
    function sgftn(object,event)
        v.IMGview = 'S';
        set(c.panel3.ax, 'Value', false);
        set(c.panel3.sg, 'Value', true);
        set(c.panel3.cr, 'Value', false);
        IMGLv1process
        IMGLv2process            
        IMGView   
    end

    set(c.panel3.cr, 'Callback', @crftn);
    function crftn(object,event)
        v.IMGview = 'C';
        set(c.panel3.ax, 'Value', false);
        set(c.panel3.sg, 'Value', false);
        set(c.panel3.cr, 'Value', true);
        IMGLv1process
        IMGLv2process            
        IMGView   
    end

    set(c.panel3.d4, 'Callback', @d4ftn);
    function d4ftn(object,event)
        v.flag.view4thdim = get(c.panel3.d4,'Value');
        set(c.panel3.d5,'Value',false);
        v.flag.view5thdim = false;
        if v.IMGdim <= 3
            v.flag.view4thdim = false;
            set(c.panel3.d4,'Value',false);
        end
        IMGLv2process
        IMGView        
    end

    set(c.panel3.d5, 'Callback', @d5ftn);
    function d5ftn(object,event)
        v.flag.view5thdim = get(c.panel3.d5,'Value');
        set(c.panel3.d4,'Value',false);
        v.flag.view4thdim = false;
        if v.IMGdim <= 4
            v.flag.view5thdim = false;
            set(c.panel3.d5,'Value',false);
        end
        IMGLv2process
        IMGView        
    end

    set(c.panel4.imgnum, 'Callback', @imgnumftn);
    function imgnumftn(object,event)
        dispr = v.disprange(:,get(c.panel4.imgnum,'Value'));
        set(c.panel4.minval, 'String', sprintf('%3.2f',dispr(1)));
        set(c.panel4.maxval, 'String', sprintf('%3.2f',dispr(2)));
    end

    set(c.panel4.minval, 'Callback', @minvalftn);
    function minvalftn(object,event)
        imgidx = get(c.panel4.imgnum,'Value');
        v.disprange(1,imgidx) = str2double(get(c.panel4.minval, 'String'));
        if v.disprange(1,imgidx) >= v.disprange(2,imgidx)
            v.disprange(1,imgidx) = v.disprange(2,imgidx) - 1e-5;
        end
        IMGView
    end

    set(c.panel4.maxval, 'Callback', @maxvalftn);
    function maxvalftn(object,event)
        imgidx = get(c.panel4.imgnum,'Value');
        v.disprange(2,imgidx) = str2double(get(c.panel4.maxval, 'String'));
        if v.disprange(1,imgidx) >= v.disprange(2,imgidx)
            v.disprange(2,imgidx) = v.disprange(1,imgidx) + 1e-5;
        end
        IMGView
    end

    set(c.panel4.auto, 'Callback', @autoftn);
    function autoftn(object,event)
        var = load_var();
        imgidx = get(c.panel4.imgnum,'Value');
        v.disprange(1,imgidx) = min(min(min(min(min(var.IMGLv2(:,:,:,:,:,imgidx))))));
        v.disprange(2,imgidx) = max(max(max(max(max(var.IMGLv2(:,:,:,:,:,imgidx))))));
        set(c.panel4.minval, 'String', sprintf('%3.2f',v.disprange(1,imgidx)));
        set(c.panel4.maxval, 'String', sprintf('%3.2f',v.disprange(2,imgidx)));        
        IMGView
    end

    set(c.panel4.cmappop, 'Callback', @cmappopftn);
    function cmappopftn(object,event)      
        IMGView
    end

%     set(c.panel4.datatip, 'Callback', @datatipftn);
%     function datatipftn(object,event)      
%         if get(c.panel4.datatip, 'Value')
%             dcm = datacursormode;
%             dcm.Enable = 'on';
%             dcm.UpdateFcn = @datatipdispftn;
%         else
%             dcm.Enable = 'off';
%         end
%     end
    function txt = datatipdispftn(~, info)
        x = info.Position(1);
        y = info.Position(2);
        var = load_var();
        value = IMGLv3(y,x,var.s);
        imgidx = IMGLv3idx(y,x,var.s);
        value = value*(v.disprange(2,imgidx)-v.disprange(1,imgidx))+v.disprange(1,imgidx);
        if value < 0.001
            txt = sprintf("value = %0.3e", value);
        else
            txt = sprintf("value = %4.3f", value);
        end
    end
    dcm = datacursormode;
    dcm.UpdateFcn = @datatipdispftn;
    datacursormode off

end
