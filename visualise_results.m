function visualise_results(Meas,Wall,Pos,En,Truth,t_m)
%%function to visualise results as described below. 
%Inputs: 
%%% Wall: structure with the geometry of the wall and
%%% time-stepping for the heat transfer model.
%%% Meas: structure that contains synthetic measurements of surface heat 
%%% flux and near air temperatures, standard deviations of surface heat 
%%% flux measurements, number of measurements within each assimilation subinterval
%%% Pos: structure with statistics of the posterior ensemble of (k, c,
%%% T_0, R_I, R_E) computed via the ensemble Kalman inversion algorithm
%%% described in Appendix B. 
%%% En: structure with the prior ensemble (of size J) of (k, c, T_0, R_I, R_E).
%%% Truth: structure that contains the true thermophysical
%%% properties used to generate synthetic data.
%%% t_m: assimilation time at which quantities of interest are visualised


n=Pos.En.n;

%%%plot synthetic data (surface heat fluxes and near air temperatures)
plot_data;
%%%plot posterior and true c, k and T_{0} at given time step t_m
plot_spatial_thermal(t_m,'');
%%%plot running mean and credible intervals of R_I, R_E, C-value and
%%%U-value within the assimilation interval [0,t_m]. This also includes the
%%%running estimate computed via the average method
plot_running_stats(t_m);
%%%plot densities of R_I, R_E, C-value and U-value computed at time t_m
plot_densities(t_m);
%%%plot predictive surface heat flux priors on the interval [T,2T]
plot_predictive_priors();
%%%plot predictive surface heat flux posteriors on the interval [T,2T]
plot_predictive_posterior()


    function plot_data()
        %%%%%%%%%%%%%%%%%%plot synthetic data%%%%%%%%%%%%%%%%%%
        M=length(Meas.HF.ext);
        mi=min(min(Meas.HF.ext),min(Meas.HF.int));
        ma=max(max(Meas.HF.ext),max(Meas.HF.int));
        lab='Synthetic heat flux measurements';
        xData = linspace(0,M,M)*5/60/24;
        xData2 = linspace(-M/2,M,3/2*M)*5/60/24;
        figure('Position', [100, 100, 1049, 295*2]);
        subplot(2,1,2);
        plot(xData,Meas.HF.int(1:M),'r','linewidth',1.1)
        hold on
        plot(xData,Meas.HF.ext(1:M),'-.k','linewidth',1.1)
        axis([xData(1),xData(end),0.9*mi,1.2*ma])
        h=legend('Internal','External');
        set(h,'fontsize',22,'interpreter','latex','location','northeast')
        xlabel('Time [days]','fontsize',15)
        ylabel('Heat flux [W m^{-2}]','fontsize',15)
        h=title(lab);
        set(h,'fontsize',20,'interpreter','latex');
        
        subplot(2,1,1);
        plot(xData2,Meas.T_all.int,'-r','linewidth',1.1)
        hold on
        plot(xData2,Meas.T_all.ext,'-.k','linewidth',1.1)
        axis([xData2(1),xData2(end),0,30])
        line([0 0], [0 30],'color',[0,0,0],'linestyle', '-.','linewidth',1.5);
        h=legend('Internal','External');
        set(h,'fontsize',22,'interpreter','latex','location','northeast')
        xlabel('Time [days]','fontsize',15)
        ylabel('Temperature [C]','fontsize',15)
        h=title('Synthetic temperature measurements');
        set(h,'fontsize',20,'interpreter','latex');
    end


    function plot_spatial_thermal(time,lab)        
        mi{1}=0.1e6;
        ma{1}=6e6;
        mi{2}=0;
        ma{2}=5.0;
        mi{3}=4;
        ma{3}=23;
        Nx=2^n;
        Nx_truth=length(Truth.thermal{1});
        Nx_fine=Nx_truth*8;
        hx_truth=Wall.L/Nx_truth;
        hx=Wall.L/Nx;
        x=linspace(0,Wall.L,Nx+1)';
        x_fine=linspace(0,Wall.L,Nx_fine)';
        x_truth=linspace(0,Wall.L,Nx_truth)';
        x_t=linspace(0,Wall.L,Nx_truth+1)';
        figure('Position', [0, 0, 1200, 300]);
        for j=1:3
            if (j<3)
                truth=Truth.exp_thermal{j}(1).*(x_fine<=x_truth(2));
                for i=2:Nx_truth-1
                    truth=truth+Truth.exp_thermal{j}(i).*((x_fine>x_truth(i))&(x_fine<=x_truth(i+1)));
                end
                truth=truth+Truth.exp_thermal{j}(Nx_truth).*(x_fine>x_truth(Nx_truth));
                for ii=1:3
                    stat(:,ii)=Pos.thermal{j}.stat{ii}(1,time).*(x_fine<=x(2));
                    for i=2:Nx-1
                        stat(:,ii)=stat(:,ii)+Pos.thermal{j}.stat{ii}(i,time).*((x_fine>x(i))&(x_fine<=x(i+1)));
                    end
                    stat(:,ii)=stat(:,ii)+Pos.thermal{j}.stat{ii}(Nx,time).*(x_fine>x(Nx));
                end
            else
                truth=zeros(length(x_fine),1);
                truth(1)=Truth.thermal{j}(1);
                
                for i=1:Nx_truth
                    truth=truth+Truth.thermal{j}(i)*(x_t(i+1)-x_fine)/hx_truth.*((x_fine>x_t(i))&(x_fine<=x_t(i+1)))+Truth.thermal{j}(i+1)*(x_fine-x_t(i))/hx_truth.*((x_fine>x_t(i))&(x_fine<=x_t(i+1)));
                end
                
                for ii=1:3
                    stat(:,ii)=zeros(length(x_fine),1);
                    stat(1,ii)=Pos.thermal{j}.stat{ii}(1,time);
                    for i=1:Nx
                        stat(:,ii)=stat(:,ii)+Pos.thermal{j}.stat{ii}(i,time)*(x(i+1)-x_fine)/hx.*((x_fine>x(i))&(x_fine<=x(i+1)))+Pos.thermal{j}.stat{ii}(i+1,time)*(x_fine-x(i))/hx.*((x_fine>x(i))&(x_fine<=x(i+1)));
                    end
                end
            end
            subplot(1,3,j);
            plot(x_fine,truth,'r','linewidth',1.5)
            hold on
            plot(x_fine,stat(:,3),'--k','linewidth',1.5)
            X=[x_fine',fliplr(x_fine')];
            Y=[stat(:,1)',fliplr(stat(:,2)')];
            f=fill(X,Y,[0.8 0.8 0.8]);
            set(f,'EdgeColor','none')
            plot(x_fine,stat(:,3),'--k','linewidth',1.5)
            plot(x_fine,truth,'r','linewidth',1.5)
            h1=xlabel('$$x$$  wall thickness [m]');
            set(h1,'fontsize',12,'interpreter','latex');
            axis([0,Wall.L,mi{j},ma{j}]);
            leg2=legend('truth','mean','CI');
            set(leg2,'fontsize',12,'interpreter','latex','location','northeast');
            if (j==1)
                y_lab=ylabel('Volumetric heat capacity [J m$^{-3}$ K$^{-1}$]','fontsize',10);
                if (time~=1)
                    ti=title(strcat('~~~~   ~~Posterior of $$c(x)$$,',lab));
                else
                    ti=title(strcat('Prior of $$c(x)$$'));
                end
            elseif (j==2)
                y_lab=ylabel('Thermal conductivity [W m$^{-1}$ K$^{-1}$]','fontsize',10);
                if (time~=1)
                    ti=title(strcat('Posterior of $$\kappa(x)$$,',lab));
                else
                    ti=title(strcat('Prior of $$\kappa(x)$$'));
                end
            elseif (j==3)
                if (time~=1)
                    ti=title(strcat( 'Posterior of $$T_{0}(x)$$,',lab));
                else
                    ti=title(strcat('Prior of $$T_{0}(x)$$'));
                end
                y_lab=ylabel('Initial temperature[C]','fontsize',10);
            end
            set(y_lab,'fontsize',10,'interpreter','latex');
            set(ti,'fontsize',15,'interpreter','latex');
        end
    end

    function plot_running_stats(time)
        tt = linspace(0,time,time)*5/60/24*Meas.DeltaT;
        %%% compute the running estimate computed via the average method.
        s(1)=Meas.HF.int(1)/(Meas.T_int(1)-Meas.T_ext(1));
        for i=1:round(Wall.Nt/Meas.DeltaT)-1
            s(i+1)=sum(Meas.HF.int(1:i*Meas.DeltaT))/sum(Meas.T_int(1:i*Meas.DeltaT)-Meas.T_ext(1:i*Meas.DeltaT));
        end
        figure('Position', [0, 0, 1600, 250]);
        yl{2}='$R_{E}$ [W$^{-1}$m$^2$ K)]';
        yl{1}='$R_{I}$ [W$^{-1}$m$^2$ K)]';
        ti{1}='$$R_{I}$$';
        ti{2}='$$R_{E}$$';
        
        yl{4}='U-value [W/(m$^2$ K)]';
        yl{3}='C-value [J/(m$^2$ K)]';
        ti{3}='$$\mathcal{C}$$';
        ti{4}='$$\mathcal{U}$$';
        
        for i=1:4
            ti{i}=strcat(ti{i},',~~$$h=2^{-',num2str(n),'}L$$');
        end
        for i=1:4
            subplot(1,4,i);
            low=prctile(Pos.thermal{i+3}(:,1:time),2.5);
            high=prctile(Pos.thermal{i+3}(:,1:time),97.5);
            X=[tt,fliplr(tt)];
            Y=[low,fliplr(high)];
            fill(X,Y,[0.8 0.8 0.8]);
            hold on
            plot(tt, mean(Pos.thermal{i+3}(:,1:time)),'--k','linewidth',1.0)
            plot(tt, Truth.thermal{i+3}*ones(time,1),'r','linewidth',1.0)
            y_lab=ylabel(yl{i});
            h=legend('CI','mean','truth');
            if (i==3)
                ylim([2.4e5,5.0e5])
            elseif(i==4)
                plot(tt, s(1:time),'b','linewidth',1.0);
                hold on
                ylim([1.2,2.0])
                h=legend('CI','mean','truth', 'AV');
            end
            if (i==1)
                ylim([0.05,0.25])
            elseif(i==2)
                ylim([0.02,0.12])
            end
            set(y_lab,'fontsize',15,'interpreter','latex');
            xlim([0,tt(end)])
            box on
            if (i==1)||(i==2)
                set(h,'fontsize',15,'interpreter','latex','location','northeast');
            else
                set(h,'fontsize',15,'interpreter','latex','location','southeast');
            end
            h_ti=title(ti{i});
            h1=xlabel('Time [days]');
            set(h1,'fontsize',15,'interpreter','latex');
            set(h_ti,'fontsize',20,'interpreter','latex');
        end
    end

    function plot_densities(time)
        
        x_l{4}=0.9;
        x_u{4}=2.2;
        x_l{3}=1.7e5;
        x_u{3}=4.5e5;
        x_l{2}=0.01;
        x_u{2}=0.14;
        x_l{1}=0.01;
        x_u{1}=0.3;
        
        lab{4}='$$\mathcal{U}$$ [W/(m$^2$ K)]';
        lab{3}='$$\mathcal{C}$$ [J/(m$^2$ K)]';
        lab{1}='$$R_{I}$$ [W$^{-1}$m$^2$ K)]';
        lab{2}='$$R_{E}$$ [W$^{-1}$m$^2$ K)]';
        
        leg{4}='$$\mathcal{U}^{\dagger}$$';
        leg{3}='$$\mathcal{C}^{\dagger}$$';
        leg{2}='$$R_{E}^{\dagger}$$';
        leg{1}='$$R_{I}^{\dagger}$$';
        
        ti{4}='$$\mathcal{U}$$';
        ti{3}='$$\mathcal{C}$$';
        ti{2}='$$R_{E}$$';
        ti{1}='$$R_{I}$$';
        
        
        left_color = [1 0 0];
        right_color = [0 0 1];
        loc{1}='northeast';
        loc{2}='northeast';
        loc{3}='northwest';
        loc{4}='northwest';
        c=0;
        fig = figure('Position', [0, 0, 1600, 250]);
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        for i=1:4
            c=c+1;
            subplot(1,4,i);
            [f_prior,xi_prior] = ksdensity(Pos.thermal{c+3}(:,1),'NumPoints',500,'Kernel','epanechnikov');%,'Bandwidth',200)
            yyaxis left
            AX1=plot(xi_prior,f_prior,'-.r','linewidth', 1.5);
            ylabel('Prior density','fontsize',12)
            ylim([0 max(f_prior)*1.5])
            yyaxis right
            [f_pos,xi_pos] = ksdensity(Pos.thermal{c+3}(:,time),'NumPoints',500,'Kernel','epanechnikov');%,'Bandwidth',200)
            hold on
            plot(xi_pos,f_pos,'b','linewidth', 1.5);
            ylabel('Posterior density','fontsize',12)
            line([Truth.thermal{c+3} Truth.thermal{c+3}], [0 max(f_pos)*1.5],'Color','k','linewidth',1.0,'linestyle','--');
            ylim([0 max(f_pos)*1.5])
            AX1.Color='r';
            box on
            xlim([x_l{c},x_u{c}])
            h=legend('Prior','Posterior',leg{c});
            set(h,'fontsize',12,'interpreter','latex','location',loc{c});
            h1=xlabel(lab{c});
            tit=title(strcat('~~Distribution of~',ti{c}));
            set(h1,'fontsize',12,'interpreter','latex');
            set(tit,'fontsize',17,'interpreter','latex');
        end
    end

    function plot_predictive_priors()
        Grid.Nx=2^n;
        Grid.hx=Wall.L/Grid.Nx;
        tt = linspace(Wall.Nt+1,2*Wall.Nt,Wall.Nt)*5/60/24;
        for en=1:En.J
            Fluxes=Heat_FEM(Grid.hx,Grid.Nx,2*Wall.Nt,Wall,exp(En.thermal{2}(:,en)),...'
                exp(En.thermal{1}(:,en)),En.thermal{3}(:,en),exp(En.thermal{4}(:,en)),exp(En.thermal{5}(:,en)), Meas.T_int,Meas.T_ext,1);
            HF.int(:,en)=Fluxes.int(Wall.Nt+1:2*Wall.Nt);
            HF.ext(:,en)=Fluxes.ext(Wall.Nt+1:2*Wall.Nt);
        end
        low =[ prctile(HF.int',2.5); prctile(HF.ext',2.5)];
        high =[ prctile(HF.int',97.5); prctile(HF.ext',97.5)];
        me =[ mean(HF.int'); mean(HF.ext')];
        
        HF_meas=[Meas.HF.int(Wall.Nt+1:2*Wall.Nt),Meas.HF.ext(Wall.Nt+1:2*Wall.Nt)];
        tit{1}='Prior predictions of internal heat flux';
        tit{2}='Prior predictions of external heat flux';
        
        loc{1}='northeast';
        loc{2}='northeast';
        
        
        for i=1:2
            figure('Position', [0, 0, 2600, 280]);
            mi(i)=1.05*min(min(low(i,:)),min(HF_meas(i,:)));
            ma(i)=1.05*max(max(high(i,:)),max(HF_meas(i,:)));
            plot(tt,HF_meas(:,i),'*r','markersize', 1.0);
            hold on
            plot(tt,me(i,:),'--k','linewidth', 0.75);
            X=[tt,fliplr(tt)];                
            Y=[low(i,:),fliplr(high(i,:))];   
            fill(X,Y,[0.8 0.8 0.8]);            
            plot(tt,HF_meas(:,i),'*r','markersize', 1.0);
            plot(tt,me(i,:),'--k','linewidth', 0.75);
            if (i==1)
                axis([tt(1),tt(end),0,40])
            else
                axis([tt(1),tt(end),-20,60])
            end
            h=legend('Measurements','predictions (mean)','predictions (CI)');
            set(h,'fontsize',12,'interpreter','latex','location',loc{i});
            xlabel('Time [days]','fontsize',12)
            ylabel('Heat flux [W m^{-2}]','fontsize',12)
            ti=title(tit{i});
            set(ti,'fontsize',17,'interpreter','latex');%,'location','northwest')
            
        end
        
    end

    function plot_predictive_posterior()
        
        tt = linspace(Wall.Nt+1,2*Wall.Nt,Wall.Nt)*5/60/24;
        HF=[Meas.HF.int(1+Wall.Nt:2*Wall.Nt),Meas.HF.ext(1+Wall.Nt:2*Wall.Nt)];        
        tit{1}=strcat( 'Posterior predictions of internal heat flux,~~',strcat('~~$h=L 2^{-',num2str(n),'}$'));
        tit{2}=strcat( 'Posterior predictions of external heat flux,~~',strcat('~~$h=L 2^{-',num2str(n),'}$'));
        
        for i=1:2
            figure('Position', [0, 0, 2600, 280]);
            plot(tt,HF(:,i),'*r','markersize', 1.0);
            hold on
            plot(tt,Pos.Predictive.mean(i,:),'--k','linewidth', 0.75);
            X=[tt,fliplr(tt)];                
            Y=[Pos.Predictive.low(i,:),fliplr(Pos.Predictive.high(i,:))];     
            fill(X,Y,[0.8 0.8 0.8]);                
            plot(tt,HF(:,i),'*r','markersize', 1.0);
            plot(tt,Pos.Predictive.mean(i,:),'--k','linewidth', 0.75);
            axis([tt(1),tt(end),-20,60])
            h=legend('Measurements','predictions (mean)','predictions (CI)');
            set(h,'fontsize',12,'interpreter','latex','location','northeast');
            xlabel('Time [days]','fontsize',12)
            ylabel('Heat flux [W m^{-2}]','fontsize',12)
            ti=title(tit{i});
            set(ti,'fontsize',17,'interpreter','latex');
            
        end
    end

end
