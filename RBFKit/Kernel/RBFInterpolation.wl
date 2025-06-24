(* ::Package:: *)

(* ::Title:: *)
(*RBF Interpolation Package*)


(* ::Subtitle:: *)
(*A Mathematica package to perform Radial Basis Function interpolation with modular components*)


BeginPackage["RBFKit`RBFInterpolation`",{"RBFKit`RBFutils`"}];

RBFInterpolation::usage = 
  "RBFInterpolation[data, rbfFunc, opts] interpolates the given data using radial basis functions. \
  The `rbfFunc` can be any function of the distance r, making this function general.";

(* Exported symbols *)
RBFInterpolation;


Begin["`Private`"];

(* Generate Monomials Module *)
generateMonomials[vars_List, deg_Integer] /; deg >0 := Module[
  {numVars, exponentVectors},
  numVars = Length[vars];
  exponentVectors = Flatten[
    Table[FrobeniusSolve[ConstantArray[1, numVars], k], {k, 0, deg}],
    1
  ];
  (Times @@ Power[vars, #]) & /@ exponentVectors
];

generateMonomials[vars_List, 0] :={}


(* RBF Basis Functions Module *)
rbfBasisFunctions[centers_, rbfFunc_] := Table[
  With[{c = centers[[i]]},
    (rbfFunc[Norm[{##} - c]] &)
  ],
  {i, Length[centers]}
];

(* RBF Design Matrix Module *)
designMatrixRBF[data_, vars_, rbfFuncs_] := DesignMatrix[
  data,
  Through@*rbfFuncs @@ vars,
  vars,
  IncludeConstantBasis -> False
];

(* Polynomial Design Matrix Module *)
designMatrixPoly[data_, vars_, deg_] /; deg > 0 := DesignMatrix[
  data,
  generateMonomials[vars, deg],
  vars,
  IncludeConstantBasis -> False
];

designMatrixPoly[data_, vars_, deg_] /; deg <= 0 := {};

(* Assemble System Matrix Module *)
assembleSystemMatrix[A_, P_] := PadRight@Join[Join[A, P, 2], Transpose[P]];


(* Create Interpolation Function Module *)
createInterpolationFunction[coefficients_, data_, vars_, rbfFuncs_, monomialBasis_] /; Length[monomialBasis]>0 := Module[
  {n, polyLength, lambdaCoeffs, polyCoeffs, interpolationFunction},
  n = Length[data];
  polyLength = Length[monomialBasis];
  lambdaCoeffs = coefficients[[1 ;; n]];
  polyCoeffs = coefficients[[n + 1 ;;]];
  
  interpolationFunction = Function[
    Evaluate[vars],
    Module[{rbfSum, polySum},
      rbfSum = Sum[
        lambdaCoeffs[[i]] * rbfFuncs[[i]][Sequence @@ vars],
        {i, n}
      ];
      polySum = polyCoeffs . (monomialBasis /. Thread[vars -> vars]);
      rbfSum + polySum
    ]
  ];
  
  interpolationFunction
];


createInterpolationFunction[coefficients_, data_, vars_, rbfFuncs_, monomialBasis_]/; Length[monomialBasis]==0  := Module[
  {n,  lambdaCoeffs, interpolationFunction,rbf},
  n = Length[data];
  lambdaCoeffs = coefficients[[1 ;; n]];
  rbf = Through[rbfFuncs[Sequence @@ vars]].lambdaCoeffs;
  interpolationFunction = Function[
    Evaluate[vars],
    Module[{rbfeval},
      Evaluate@rbf
    ]
  ];
  
  interpolationFunction
];

(* Main RBF Interpolation Function *)
Options[RBFInterpolation] = {"PolynomialDegree" -> 0};

RBFInterpolation[data_, rbfFunc_, opts : OptionsPattern[]] := Module[
  {
    deg = OptionValue["PolynomialDegree"],
    vars, monomialBasis, A, P, Phi, b, coefficients,
    interpolationFunction, n, polyLength, centers, rbfFuncs
  },
  
  (* Determine the number of variables *)
  vars = Array[Symbol["x" <> ToString[#]] &, Dimensions[data][[2]] - 1];
  
  (* Generate monomial basis *)
  monomialBasis = generateMonomials[vars, deg];
  
  (* Extract input variables and function values *)
  centers = data[[All, 1 ;; -2]];
  n = Length[data];


  (* Generate RBF basis functions *)
  rbfFuncs = rbfBasisFunctions[centers, rbfFunc];
  polyLength = Length[monomialBasis];

  (* Construct design matrices *)
  A = Outer[rbfFunc[Norm[#1-#2]] &,centers,centers,1];
  P = designMatrixPoly[data, vars, deg];
  
  (* Assemble system matrix *)
  Phi = ArrayFlatten[{{A, P}, {Transpose[P],0}}];
  (* Phi = assembleSystemMatrix[A, P];*)
  
  (* Right-hand side vector *)
  b = Join[data[[All, -1]], ConstantArray[0, polyLength]];
  
  (* Solve for coefficients *)
  coefficients = LinearSolve[Phi, b];
  
  (* Create the interpolation function *)
  interpolationFunction = createInterpolationFunction[
    coefficients, data, vars, rbfFuncs, monomialBasis
  ];
  

  RBFInterpolant[
    <|
      "Function" -> interpolationFunction,
      "Variables" -> vars,
      "Degree" -> deg,
      "NCenters" -> Length@centers
    |>
  ]
];

End[];

EndPackage[];
