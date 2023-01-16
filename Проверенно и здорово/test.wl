(* ::Package:: *)

BeginPackage["Test`"]


intfunc;
int1;


Begin["Private`"]


func1[x_?NumericQ]:=x-RandomReal[]


func[x_?NumericQ]:=RandomReal[]x^2+func1[x]
intfunc[y_?NumericQ]:=NIntegrate[func[x],{x,func1[y],y}]
int1=intfunc[1]


End[];
EndPackage[];
