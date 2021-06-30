function net_out=update_net_g(net,g)
net.g=g;
net.W=net.g*randn(net.N)/sqrt(net.N);
net_out=net;
end