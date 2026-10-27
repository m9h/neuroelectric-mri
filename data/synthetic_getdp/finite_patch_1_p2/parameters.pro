TissueDegree = 2;
Group {
Bulk = Region[{1}];
Electrodes = Region[{101,102}];
Domain = Region[{Bulk, Electrodes}];
E101 = Region[101];
E102 = Region[102];
}
Function {
sigma[Region[1]] = Tensor[0.5,0,0,0,0.5,0,0,0,0.5];
contact[E101] = 0.002;
injection[E101] = 0.001;
area[E101] = 0.0011999999999999999;
contact[E102] = 0.0050000000000000001;
injection[E102] = -0.001;
area[E102] = 0.0011999999999999997;
}
Constraint { { Name Gauge; Type Assign; Case { { Region E102; Value 0.; } } } }
