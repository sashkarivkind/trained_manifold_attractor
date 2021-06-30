function u_out=create_step_input(u_on, nsteps,switch_step)

u_out=zeros(2,nsteps); 
u_out(:,switch_step:nsteps)=diag(u_on)*ones(size(u_out(:,switch_step:nsteps)));
