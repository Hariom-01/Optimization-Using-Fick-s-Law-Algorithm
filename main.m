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
close all
clear all
clc
fitfun = @Chung_Reynolds;
dim=30;
Max_iteration=1000;
%SearchAgents_no=30;
NoMolecules=30;
lb=-100;
ub=100;
tlt='Chung Reynolds';
i=1;

[Xfood, Xvalue,CNVG] = FLA(NoMolecules,Max_iteration,lb,ub,dim,fitfun)
figure,
plot(CNVG,'Color', 'r')
xlim([1 1000]);
