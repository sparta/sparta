
 
  Create_Sparta_input

In order to compare continuum vs Sparta calculations, it's necessary to 
calculate the thermodynamic representation that is the basis of each representation in the
limit of thermodynamic equilibrium.


In this program, Zuzax input file for air using a Stat Mech thermo calculator is used to produce
the current Sparta species input file.

The implementation needs work. There are several unimplemented features that need fleshing out.
However, the basic premise is exercised that the two program must use the same input data if they
are to create the same long time thermodynamic results.

The input file for Zuzax is named airSM.xml. It contains a description of the species and some of
their StatMech input. This file references an internal database of StatMech quantitities that Zuzax 
maintains within its source code.
The output is contained in the file air.species





There are several important omissions within Sparta. The electronic contributions due to the
electronic partition functions are not represented within Sparta. This has large ramifications
for temperature above 5000 K. It also effects the thermodynamics for atoms and molecules which 
have closely lying excited electronic states such as the O and N atom. Thus, there are problems with
representing air at even moderate temperatures due to the lack of 

Another issue is the lack of anharmonic vibrational models for diatomic species. The N2 and O2 molecules
have significant anharmonicity that materially affect their thermodynamic functions. These effects aren't 
taken into account within Sparta. 






