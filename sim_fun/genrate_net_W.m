function net_out=genrate_net_W(net)
net.W=randn(net.N)/sqrt(net.N);
net_out=net;
end