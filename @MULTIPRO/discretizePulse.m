function pulseArray=discretizePulse(pulse,delta_t)
w_d=delta_t;
pulseArray.delta_t=delta_t;
for i_n=1:pulse.Nacq
    delta=pulse.delta(i_n); % gradient width
    tau=pulse.Delta(i_n); % gradient separation
    n=pulse.n(i_n);     
    w_t=w_d:w_d:(delta+tau);
    n_impulse=length(w_t);
    g_max=pulse.G(i_n);
    g=zeros(1,n_impulse);
    switch pulse.shape(i_n)
        case "tpgse" 
            trise=pulse.trise(i_n);
            tp=pulse.tp(i_n);
            dg=g_max/(trise/w_d);
            t_0=floor(trise/w_d);
            g(1:t_0)=w_t(1:t_0)*dg/w_d;
            t_1=floor((trise+tp)/w_d);
            g(t_0:t_1)=g_max;
            t_2=floor((trise*2+tp)/w_d)-2;
            g(t_1:t_2)=g_max-w_t(1:(t_2-t_1+1))*dg/w_d;
            g((length(w_t)-t_2+1):length(w_t))=-1.*g(1:t_2);    
        case "pgse"
            t_1=floor(pulse.delta(i_n)/w_d);
            g(1:t_1)=g_max;
            g((length(w_t)-t_1+1):length(w_t))=-1.*g(1:t_1);            
        case "tcos"
            trise=pulse.trise(i_n);
            tp=pulse.tp(i_n);
            t3=tp+trise/2;
            dg=g_max/(trise/w_d);
            %% start of the diffusion gradient
            t_0=floor(trise/w_d);
            g(1:t_0)=w_t(1:t_0)*dg/w_d;
            t_1=floor((trise+tp)/w_d);
            g(t_0:t_1)=g_max;
            t_2=floor((trise*3+tp)/w_d);
            g(t_1:t_2)=g_max-w_t(1:(t_2-t_1+1))*dg/w_d;
            t_3=floor((trise*3+tp+t3)/w_d);
            g(t_2:t_3)=-g_max;
            %% middle of the diffusion gradient
            if n>1
                tt_0=floor(t3/w_d);
                gg(1:tt_0)=-g_max;
                tt_1=floor((trise*2+t3)/w_d);
                gg(tt_0:tt_1)=-g_max+w_t(1:(tt_1-tt_0+1))*dg/w_d;
                tt_2=floor((trise*2+t3*3)/w_d);
                gg(tt_1:tt_2)=g_max;
                tt_3=floor((trise*4+t3*3)/w_d);
                gg(tt_2:tt_3)=g_max-w_t(1:(tt_3-tt_2+1))*dg/w_d;
                tt_4=floor((trise*4+t3*4)/w_d);
                gg(tt_3:tt_4)=-g_max;
                g((t_3+1):(t_3+tt_4*(n-1)))=repmat(gg(1:tt_4),[1 n-1]);
                t_5=t_3*2+tt_4*(n-1);
                g((t_3+tt_4*(n-1)+1):t_5)=g(t_3:-1:1);
            else
                t_5=t_3*2;
                g((t_3+1):t_5)=g(t_3:-1:1);
            end
            %% end of the diffusion gradient
            g((length(w_t)-t_5+1):length(w_t))=-1.*g(1:t_5);
        case "cos" % cos 
            t_0=floor(delta/n/w_d);
            gg(1:t_0)=g_max.*cos(w_t(1:t_0)/w_t(t_0)*2*pi);
            g(1:t_0*n)=repmat(gg(1:t_0),[1 n]);
            g((length(w_t)-t_0*n+1):length(w_t))=-1.*g(1:t_0*n);
    end
    pulseArray.q(i_n,1:n_impulse)=g*w_d;
    pulseArray.n(i_n)=n_impulse;
    pulseArray.delta_q(i_n,1:n_impulse)=g./(g_max/(trise/w_d));
end
pulseArray.nq=trise/w_d;
end