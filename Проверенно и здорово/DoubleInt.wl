(* ::Package:: *)

BeginPackage["DoubleInt`",{"CutOffs`"}];
\[Omega]uh=\[Omega]27=f2\[Omega]@50.789979926671045`;
\[Omega]0=f2\[Omega]@82.5;
\[Omega]s=\[Omega]0-\[Omega]uh;
xl=leftCutOff[\[Omega]uh];
xr=rightCutOff[\[Omega]uh];
Dplus[x_]=Derivative[1,0,0][Duh][qplus[\[Omega]uh,x],\[Omega]uh,x];
Dminus[x_]=-Derivative[1,0,0][Duh][qminus[\[Omega]uh,x],\[Omega]uh,x];
(* \[CapitalLambda] *)
\[Delta]\[Psi][x_?NumericQ]:=NIntegrate[qplus[\[Omega]uh,\[Xi]]-qminus[\[Omega]uh,\[Xi]],{\[Xi],xl,x}];
func[\[Xi]_?NumericQ]:=1/Dplus[\[Xi]]+1/Dminus[\[Xi]]+(2 Sin[\[Delta]\[Psi][\[Xi]]])/Sqrt[Dplus[\[Xi]]Dminus[\[Xi]]];
\[CapitalIota]=Re@NIntegrate[func[\[Xi]],{\[Xi],xl,xr}];


EndPackage[];
