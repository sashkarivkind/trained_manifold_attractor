function opt=nom_opt_assigner(opt,nom)
%assignes nom values fields to the fields that are missing in opt 
%example:
% nom_opt_assigner(struct('puki',8,'chaku','zumba'),struct('puki',80,'muki',9))
% ans = 
%      puki: 8
%     chaku: 'zumba'
%      muki: 9

nomfields=fields(nom);
for ff=1:length(nomfields)
    this_field=nomfields{ff};
    if ~isfield(opt,this_field)
        opt.(this_field)=nom.(this_field);
    end
end