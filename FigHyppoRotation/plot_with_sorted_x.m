function varargout = plot_with_sorted_x(varargin)
[arg1,ii]=sort(varargin{1});
arg2 = varargin{2}(ii);

 [varargout{1:nargout}]  = plot(arg1,arg2,varargin{3:end});