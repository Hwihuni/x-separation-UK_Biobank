function [v] = imshow_3df_parse_input(varargin)
    varargin = varargin{1};
    
    v.IMGdim = length(size(varargin{1}));
    v.IMG(:,:,:,:,:,1) = varargin{1}; v.IMGnum = 1; 
    v.disprange(:,1) = [min(varargin{1}(:))/2, max(varargin{1}(:))/2];
    for i=2:length(varargin)
        if ischar(varargin{i})
            break
        end
        if ~isequal(size(varargin{i}), [1 2]) 
            v.IMGnum = v.IMGnum+1; v.IMG(:,:,:,:,:,v.IMGnum) = varargin{i};
            v.disprange(:,v.IMGnum) = [min(varargin{i}(:))/2, max(varargin{i}(:))/2];
        else
            v.disprange(:,v.IMGnum) = varargin{i};
        end
    end
    
    for k=i:length(varargin)
        if strcmpi(varargin{k},'range')
            v.disprange = repmat(varargin{k+1}', [1,v.IMGnum]);
        end       
        if strcmpi(varargin{k},'slice')
            v.IMG = v.IMG(:,:,varargin{k+1},:,:,:);
        end               
    end
    
    if ~ isreal(v.IMG)
        v.IMG = abs(v.IMG);
        v.disprange = abs(v.disprange);
    end
    
   v.a.IMGLv0 = v.IMG;
%     v.a.IMGLv0 = rot90(v.IMG,2);
%     v.a.IMGLv0 = flipud(v.IMG);
    v.s.IMGLv0 = flip(permute(v.IMG, [3 1 2 4 5 6]),1);
    v.c.IMGLv0 = flip(permute(v.IMG, [3 2 1 4 5 6]),1);
    v.a.IMGLv1 = v.a.IMGLv0;
    v.s.IMGLv1 = v.s.IMGLv0;
    v.c.IMGLv1 = v.c.IMGLv0;

    v.a.sno    = size(v.a.IMGLv0,3);  v.a.s = round(v.a.sno/2);
    v.s.sno    = size(v.s.IMGLv0,3);  v.s.s = round(v.s.sno/2);
    v.c.sno    = size(v.c.IMGLv0,3);  v.c.s = round(v.c.sno/2);
    v.sno4     = size(v.a.IMGLv0,4);  v.s4  = round(v.sno4 /2);
    v.sno5     = size(v.a.IMGLv0,5);  v.s5  = round(v.sno5 /2);

    v.a.CROP  = [1,size(v.a.IMGLv0,1),1,size(v.a.IMGLv0,2)];
    v.s.CROP  = [1,size(v.s.IMGLv0,1),1,size(v.s.IMGLv0,2)];
    v.c.CROP  = [1,size(v.c.IMGLv0,1),1,size(v.c.IMGLv0,2)];

    v.a.arrow = ['A','P','R','L'];
    v.s.arrow = ['S','I','A','P'];
    v.c.arrow = ['S','I','R','L'];

    v.flag.flatten = false;  
    v.flag.view4thdim = false;  
    v.flag.view5thdim = false;       

    v.IMGview = 'A';              
end

