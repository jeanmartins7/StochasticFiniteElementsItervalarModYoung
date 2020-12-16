%Start script

function intrval_solution
global A Amp m M3 freq N
clear all

clc

close all

global omega

%%
%Intervalos dos Parametros

%------------------------------

%  - modulo de Young
p_1_m=7e10;

p_1_l=p_1_m-.05*p_1_m;

p_1_r=p_1_m+.05*p_1_m;

p_1_I=[p_1_l,p_1_r];


%Numero de intervalos

n_I=3;




%Options GA

options = gaoptimset('PopulationSize',n_I*10,'PopulationType','doubleVector');
% Analysis Intervalar

omega_f=0.5;% Frequencia final

delta_omega=0.001;

omega_i=0;

iiii=1;
for omega=omega_i:delta_omega:omega_f
    
    lb = [p_1_l];
    
    ub = [p_1_r];
    
    
    %Interval Analysis for all Parameters
    
    [x_min,fval_min] = ga(@func1,n_I,[],[],[],[],lb,ub,[],options);
    
    [x_max,fval_max] = ga(@func1_max,n_I,[],[],[],[],lb,ub,[],options);
    
    y(iiii)=func1([p_1_m]);
    
    y_min(iiii)=fval_min;
    
    y_max(iiii)=-fval_max;
    
    p_min(iiii,:)=x_min;
    
    p_max(iiii,:)=x_max;
    
    
    
    iiii=iiii+1
    
end
%Plot of interval Analysis of all parameters

plot([omega_i:delta_omega:omega_f],y,'r');

hold on

plot([omega_i:delta_omega:omega_f],y_min,'b');

hold on

plot([omega_i:delta_omega:omega_f],y_max,'g ');

title('Total interval Analysis','Interpreter','latex','fontsize',16)

ylabel('Deslocamento','Interpreter','latex','fontsize',22)

xlabel('Tempo','Interpreter','latex','fontsize',22)

box on

saveas(gcf,'interval_solution2gdl_massa2.tif');
saveas(gcf,'interval_solution2gdl_massa2.eps');
save('y_max', 'y_min', 'y');

end

function z=func1(GG)
global A Amp m freq omega M
N = 5;
rho = 2700;
Le = 5;
h = 0.5;
b = 0.5;
Ae = b*h;
I = (b*h*(b^2 + h^2))/24;

o=omega;

if o == 0
    z = 0;
    
else
    E = GG(1);
    uo = 0.001;
    [K,M] = beam_ef_creating(N,rho,E,Ae,I,Le);
    alfa=1e-2;
    beta=1e-5;
    C=(alfa*M+beta*K);
    z = time_response(M,K,C,o,uo);
end


end

function z=func1_max(GG)
global A Amp m freq omega M

N = 5;
rho = 2700;
Le = 5;
h = 0.5;
b = 0.5;
Ae = b*h;
I = (b*h*(b^2 + h^2))/24;

o=omega;
if o == 0
    z = 0;
    
else
    E = GG(1);
    uo = 0.001;
    [K,M] = beam_ef_creating(N,rho,E,Ae,I,Le);
    alfa=1e-2;
    beta=1e-5;
    C=(alfa*M+beta*K);
    z = time_response(M,K,C,o,uo);
    z=-z;
end


end

