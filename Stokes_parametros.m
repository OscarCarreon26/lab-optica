% Calculo de parametros de Stokes
% Oscar Alberto Carreon Gastelum - A01741876
% Ultima Modificacion: 14/04/2023 - 15:16

clc; clear all; close all;

n_p = 2^8;     %Numero de puntos para la discretizacion
aux = (-n_p/2:n_p/2-1)  * 2/n_p ;
[X,Y] = meshgrid(aux,aux);
[phi,r] = cart2pol(X,Y);

% Representacion del campo optico
m = 0;
w0 = 1;
E = (r./w0).^abs(m) .* exp(-r.^2/w0^2) .* exp(1i * m .*phi);

% Polarizacion del haz incidente

Pol0 = 1/sqrt(2) * [1;-1i];

E0x = Pol0(1).*E;
E0y = Pol0(2).*E;

E0 = sqrt(abs(E0x).^2 + abs(E0y).^2);
I0 = sum(sum(E0(:).^2)) * (aux(2)-aux(1)).^2;

I0_n = I0/I0; % Intensidad del haz polarizado circularmente a la derecha nomalizado

% Retardador de media onda, con el eje rapido orientado de manera vertical
Retardador = [cos(2*pi/2), sin(2*pi/2); sin(2*pi/2),-cos(2*pi/2)];    % Retardadro 1/2 onda con el eje rapido orientado de forma vertical

n = 100;                    % Numero de divisiones para el angulo
theta = linspace(0,2*pi,n); % Angulo para rotar la posicion original del retardador
I_ret_theta = zeros(1,n);

for ii = 1:1:n
    Ret_theta = [cos(-theta(ii)), sin(-theta(ii)); -sin(-theta(ii)), cos(-theta(ii))] * Retardador *[cos(theta(ii)), sin(theta(ii)); -sin(theta(ii)), cos(theta(ii))];
    
    E_ret_theta_x = Ret_theta(1,1)*E0x + Ret_theta(1,2)*E0y;
    E_ret_theta_y = Ret_theta(2,1)*E0x + Ret_theta(2,2)*E0y;

    E_ret_theta = sqrt(abs(E_ret_theta_x).^2 + abs(E_ret_theta_y).^2);

    I_ret_theta(1,ii) = (sum(sum(E(:).^2)) * (aux(2)-aux(1)).^2);

    
    
end









