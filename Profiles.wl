(* ::Package:: *)

BeginPackage["Profiles`"];
den::usage="den[x] Density profile";
denPlot::usage="denPlot[xl,xr] Density profile plot";
B::usage="B[x] Magnetic field profile";
Ti::usage="Ti[x] Ti profile";
Te::usage="Te[x] Te profile";


Begin["`Private`"];

R0=88;a=25;
den[x_]=1.6`*10^13 (1-(x/a)^2)^1+0.687`*10^13 Exp[-(x-8)^2/3.6`^2];
denPlot[x1_,x2_]:=Plot[den[x],{x,x1,x2}]

B0=13975;
B[x_]=B0 R0/(R0+x);

Te0=1800;
Te[x_]=Te0 (1-x^2/a^2)^2

Ti0=350;
Ti[x_]=Ti0 (1-x^2/a^2)^2

(*SetAttributes[{Ti,Te,B,den},{Protected}];*)

End[];
EndPackage[];



