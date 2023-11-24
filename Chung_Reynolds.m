%_________________________________________________________________________%
%  Fick's Law Algorithm (FLA) source codes version 1.0                    %
%                                                                         %
%  Developed in MATLAB R2021b                                             %
%                                                                         %
%  Coresponding Author:  Abdelazim G. Hussien                             %
%                                                                         %
%                                                                         %
%         e-Mail: abdelazim.hussien@liu.se                                %
%                 aga08@fayoum.edu.eg                                     %
%                                                                         %
%                                                                         %
%   Main paper: Fatma Hashim, Reham R Mostafa, Abdelazim G. Hussien,      %
%                     Seyedali Mirjalili, & Karam M. Sallam               %
%               Knowledge-based Systems                                   %
%                                                                         %
%_________________________________________________________________________%
function objval = Chung_Reynolds(x)%[-5.12,5.12]
 n = size(x,2); 
f = 0;
for i = 1:n;
    f = f+(x(i).^2).^2;
end
objval = f;
