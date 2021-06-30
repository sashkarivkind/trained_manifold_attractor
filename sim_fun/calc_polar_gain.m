function G_polar = calc_polar_gain(G_cartesian,theta)

G_polar.cell = {};
G_polar.flat = {};
% G_polar.flat = zeros(length(theta_vec,4));
for ii=1:length(G_cartesian) %equals also to length theta
    G_polar.flat{ii} = cell(2,2);
    for oo=1:2
        for uu=1:2
            G_polar.flat{ii}{oo,uu}=...
                zeros(length(G_cartesian{ii}),1);
        end
    end
    R=[cos(theta(ii)), -sin(theta(ii));sin(theta(ii)), cos(theta(ii))];
    % G_polar.flat = zeros(length(theta_vec,4));
    %     R=inv(R);
    for jj=1:length(G_cartesian{ii})
        G_polar.cell{ii}{jj} = inv(R)*G_cartesian{ii}{jj}*R;
        for oo=1:2
            for uu=1:2
                G_polar.flat{ii}{oo,uu}(jj)=G_polar.cell{ii}{jj}(oo,uu);
            end
        end
        % G_polar.flat(ii,:) = G_polar.cell{ii}(:)';
    end
end