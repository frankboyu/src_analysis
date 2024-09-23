TLorentzVector test_beam(0.0, 0.0, 9.0, 9.0);
TLorentzVector test_target(0.0, 0.0, 0.0, mass_deuteron);

TLorentzVector test1 = boost_lorentz_vector(test_beam, -(test_beam + test_target).BoostVector());
TLorentzVector test2 = boost_lorentz_vector(test_target, -(test_beam + test_target).BoostVector());

cout << "Test1: " << test1.Px() << " " << test1.Py() << " " << test1.Pz() << " " << test1.E() << endl;
cout << "Test2: " << test2.Px() << " " << test2.Py() << " " << test2.Pz() << " " << test2.E() << endl;

TLorentzVector test_cm = test_beam + test_target;
test_beam.Boost(-test_cm.BoostVector());
test_target.Boost(-test_cm.BoostVector());

cout << "Test1: " << test_beam.Px() << " " << test_beam.Py() << " " << test_beam.Pz() << " " << test_beam.E() << endl;
cout << "Test2: " << test_target.Px() << " " << test_target.Py() << " " << test_target.Pz() << " " << test_target.E() << endl;