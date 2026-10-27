TissueDegree = 2;
Group {
Bulk = Region[{1,2}];
Electrodes = Region[{101,102}];
Domain = Region[{Bulk, Electrodes}];
E101 = Region[101];
E102 = Region[102];
}
Function {
sigma[Region[1]] = Tensor[0.20000000000000001,0,0,0,0.20000000000000001,0,0,0,0.20000000000000001];
sigma[Region[2]] = Tensor[0.80000000000000004,0,0,0,0.80000000000000004,0,0,0,0.80000000000000004];
contact[E101] = 0.002;
injection[E101] = 0.001;
area[E101] = 0.0047999999999999996;
contact[E102] = 0.0050000000000000001;
injection[E102] = -0.001;
area[E102] = 0.0047999999999999996;
}
Constraint { { Name Gauge; Type Assign; Case { { Region E102; Value 0.; } } } }
