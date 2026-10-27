// Independent complete-electrode weak form. SI units; I is inward-positive.
// Parameters, physical regions and the voltage gauge are generated separately.
Include "parameters.pro";

Jacobian {
  { Name Volume; Case { { Region All; Jacobian Vol; } } }
  { Name Surface; Case { { Region All; Jacobian Sur; } } }
}
Integration {
  { Name Gauss;
    Case { { Type Gauss;
      Case {
        { GeoElement Tetrahedron; NumberOfPoints 4; }
        { GeoElement Triangle; NumberOfPoints 6; }
      }
    } }
  }
}
FunctionSpace {
  { Name Tissue; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Node;
        Support Domain; Entity NodesOf[All]; }
      If(TissueDegree == 2)
        { Name se; NameOfCoef ue; Function BF_Node_2E;
          Support Domain; Entity EdgesOf[All]; }
      EndIf
    }
  }
  { Name Metal; Type Scalar;
    BasisFunction {
      { Name sr; NameOfCoef Ur; Function BF_Region;
        Support Electrodes; Entity Region[Electrodes]; }
    }
    Constraint { { NameOfCoef Ur; EntityType Region; NameOfConstraint Gauge; } }
  }
}
Formulation {
  { Name CompleteElectrode; Type FemEquation;
    Quantity {
      { Name u; Type Local; NameOfSpace Tissue; }
      { Name U; Type Local; NameOfSpace Metal; }
    }
    Equation {
      Integral { [ sigma[] * Dof{d u}, {d u} ]; In Bulk; Jacobian Volume; Integration Gauss; }
      Integral { [ Dof{u} / contact[], {u} ]; In Electrodes; Jacobian Surface; Integration Gauss; }
      Integral { [ -Dof{U} / contact[], {u} ]; In Electrodes; Jacobian Surface; Integration Gauss; }
      Integral { [ -Dof{u} / contact[], {U} ]; In Electrodes; Jacobian Surface; Integration Gauss; }
      Integral { [ Dof{U} / contact[], {U} ]; In Electrodes; Jacobian Surface; Integration Gauss; }
      Integral { [ -injection[] / area[], {U} ]; In Electrodes; Jacobian Surface; Integration Gauss; }
    }
  }
}
Resolution {
  { Name CEM;
    System { { Name Sys; NameOfFormulation CompleteElectrode; } }
    Operation { Generate[Sys]; Solve[Sys]; SaveSolution[Sys]; }
  }
}
PostProcessing {
  { Name Results; NameOfFormulation CompleteElectrode;
    Quantity {
      { Name potential; Value { Local { [ {u} ]; In Bulk; } } }
      { Name electric; Value { Local { [ -{d u} ]; In Bulk; } } }
      { Name current; Value { Local { [ -sigma[] * {d u} ]; In Bulk; } } }
      { Name bulk_power; Value { Integral { [ sigma[] * {d u} * {d u} ];
        In Bulk; Jacobian Volume; Integration Gauss; } } }
      { Name contact_power; Value { Integral { [ ({U} - {u})^2 / contact[] ];
        In Electrodes; Jacobian Surface; Integration Gauss; } } }
      // @ELECTRODE_QUANTITIES@
    }
  }
}
PostOperation {
  { Name Export; NameOfPostProcessing Results;
    Operation {
      Print[potential, OnElementsOf Bulk, Depth 0, Format SimpleTable, File "potential.txt"];
      Print[electric, OnElementsOf Bulk, Depth 0, Format SimpleTable, File "electric.txt"];
      Print[current, OnElementsOf Bulk, Depth 0, Format SimpleTable, File "current.txt"];
      Print[bulk_power[Bulk], OnGlobal, Format SimpleTable, File "bulk_power.txt"];
      Print[contact_power[Electrodes], OnGlobal, Format SimpleTable, File "contact_power.txt"];
      // @ELECTRODE_OUTPUTS@
    }
  }
}
