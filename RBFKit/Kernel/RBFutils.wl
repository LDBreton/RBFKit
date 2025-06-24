(* Wolfram Language Package *)

BeginPackage["RBFKit`RBFutils`"]
(* Exported symbols added here with SymbolName::usage *)  
distancematrix::usage = 
  "distancematrix[rbf, Nodes, Centers] computes a distance matrix where rbf is applied to all combinations of Nodes and Centers.";
RBFInterpolant::usage = ""
Begin["`Private`"] (* Begin Private Context *) 

(* Function to compute the distance matrix using a radial basis function (rbf) *)
distancematrix[rbf_, Nodes_, Centers_] := Outer[rbf, Nodes, Centers, 1]

RBFInterpolant /: (obj_RBFInterpolant)[args__]  := obj[[1]]["Function"][args]

RBFInterpolant /: (obj_RBFInterpolant)[s_String] := obj[[1]][s]

(* Define a simple custom format to display the key information *)
Format[interp : RBFInterpolant[assoc_Association]] := RBFInterpolant[
  Panel[
    Column[{
      "RBF Interpolation Function",
      "Variables: " <> StringRiffle[ToString /@ assoc["Variables"], ", "],
      "Polynomial Degree: " <> ToString[assoc["Degree"]],
      "Number of Centers: " <> ToString[assoc["NCenters"]]
    }],
    BaseStyle -> {FontFamily -> "Helvetica", FontSize -> 12}
  ]];


End[] (* End Private Context *)

EndPackage[]