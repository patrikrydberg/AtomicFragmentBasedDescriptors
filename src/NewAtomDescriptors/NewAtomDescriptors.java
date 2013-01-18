/* 
 * Copyright (C) 2012  Patrik Rydberg <patrik.rydberg@gmail.com>
 * 
 * Contact: pry@farma.ku.dk
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package NewAtomDescriptors;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormatSymbols;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Locale;
import java.lang.System;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.SMILESWriter;
import org.openscience.cdk.io.iterator.DefaultIteratingChemObjectReader;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.io.iterator.IteratingSMILESReader;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import NewAtomDescriptors.MoleculeKU.NEWDESC_PROPERTY;

public class NewAtomDescriptors {


	public static void main(String[] arguments) throws Exception{

		// Check that the arguments (molecule files) have been given
		if (arguments.length < 1){
			System.out.println("Wrong number of arguments!" + '\n' + "Usage: java -jar 2DSASA.jar <One or more moleculeFiles>");
			System.exit(0);			
		}
		
	    String[] filenames;
    	filenames = arguments;
	    
	    File inputFile;
		String infileName;
		
		//to get the decimal signs as dot, even when running on other systems
		Locale.setDefault(new Locale("en", "US"));
		//we wamt three digits after comma for SASA
		NumberFormat formatter = new DecimalFormat("###.###");

		// Iterate over all molecule infiles (it can be a single file)
		int moleculeFileNr;
		int highestMoleculeID = 1;
		int moleculeIndex = 0;
		for (moleculeFileNr = 0; moleculeFileNr < filenames.length; moleculeFileNr++) {
			
			infileName = filenames[moleculeFileNr];	
			inputFile = new File(infileName);

			//initiate csv output
			// DecimalFormat
			DecimalFormat twoDecimalFormat = new DecimalFormat("#.##");
			twoDecimalFormat.setDecimalSeparatorAlwaysShown(false);
			DecimalFormatSymbols decformat = new DecimalFormatSymbols();
			decformat.setDecimalSeparator('.');
			decformat.setGroupingSeparator(',');
			twoDecimalFormat.setMaximumFractionDigits(2);
			twoDecimalFormat.setDecimalFormatSymbols(decformat);

		    PrintWriter outfile = null;
		    try {
				outfile = new PrintWriter(new BufferedWriter(new FileWriter(infileName + "_newatomdescriptors.csv")));
			} catch (IOException e) {
				System.out.println("Could not create CSV outfile");
				e.printStackTrace();
			}
			
		    if (!inputFile.exists()) {
			    System.err.println("File not found: " + infileName);
			    System.exit(1);
			}

			DefaultIteratingChemObjectReader reader = null;
			if (infileName.endsWith(".sdf")) {  
				reader = (IteratingMDLReader) new IteratingMDLReader(new FileInputStream(infileName), DefaultChemObjectBuilder.getInstance());
			}
			else if (infileName.endsWith(".smi")){
				reader = new IteratingSMILESReader(new FileReader(infileName), DefaultChemObjectBuilder.getInstance());
			}
			else {
				System.err.println("Filetype not supported (only .sdf and .smi are supported): " + infileName);
			    System.exit(1);
			}
			Molecule mol=null;
			IAtomContainer [] atom2endofmolMols = null;
			IAtomContainer [] atom2endofbranchMols = null;
			
			int moleculePrintIndex = 0;
			//iterate molecules
			while (reader.hasNext()){
	
				mol=(Molecule)reader.next();
				moleculeIndex ++;
				moleculePrintIndex++;
				
				
				if (moleculeIndex == 1){
					outfile.println("Molecule,Atom,Mol_bonds2end,Mol_rotablebonds,Mol_AtomCount,Mol_TPSA,Mol_TPSAperAtom,Mol_Volume,Mol_HAcount,Mol_HDcount," +
					        "Mol_PIsystemSize," +
					        "Branch_bonds2end,Branch_rotablebonds,Branch_AtomCount,Branch_TPSA,Branch_TPSAperAtom,Branch_Volume,Branch_HAcount," + 
					        "Branch_HDcount,Branch_PIsystemSize");
				}

				MoleculeKU moleculeKU;
				IAtomContainer iAtomContainer;	
				CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance());				

				iAtomContainer = AtomContainerManipulator.removeHydrogens(mol);
					
				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(iAtomContainer);

				adder.addImplicitHydrogens(iAtomContainer);
				CDKHueckelAromaticityDetector.detectAromaticity(iAtomContainer); 	
						
				moleculeKU = new MoleculeKU(iAtomContainer);	
				moleculeKU.setID(Integer.toString(highestMoleculeID));
				
				System.out.println("\n ************** Molecule " + (moleculeIndex) + " **************");
				
				moleculeKU.calculateRelativeSpan();
				moleculeKU.setSymmetryNumbers();
				
				int [] EndOfMoleculeAtoms = moleculeKU.findAtomsatEndOfMolecule();
				
				int [] EndOfBranchAtoms = moleculeKU.findAtomsatEndOfBranch();
				
				String prefix = "Mol";
				atom2endofmolMols = moleculeKU.getAtoms2EndOfMolMolecules(EndOfMoleculeAtoms,prefix);
				prefix = "Branch";
				atom2endofbranchMols = moleculeKU.getAtoms2EndOfMolMolecules(EndOfBranchAtoms,prefix);
				
				//list of properties that are set for each fragment in atom2endofmolMols:
				//RotableBondCount (my own implementation, not default CDK, also excludes amide, thioamide and sulfonamide bonds)
				//BondsToEndofMol
				//AtomCount
				//TPSA (topological polar surface area)
				//TPSAperAtom
				//Volume (from atom types and bond types)
				//HAcount (hydrogen bond acceptor count)
				//HDcount (hydrogen bond donor count)
				//PISystemSize (largest pi system, using atom counts)
				
				
				//below print output in any way you'd like
				//here is a simple example of printing all molecular properties in one go
				/*
				int testmol;
				for (testmol = 0; testmol < atom2endofmolMols.length; testmol++){
					System.out.println("Fragmentatom " + (testmol+1));
					System.out.println(atom2endofmolMols[testmol].getProperties());
					System.out.println(atom2endofbranchMols[testmol].getProperties());
				}
				*/
				
				Atom currentAtom;
				String currentAtomType;

				for(int atomIndex = 0; atomIndex < moleculeKU.getAtomCount()  ; atomIndex++ ){
					
					currentAtom = (Atom) moleculeKU.getAtom(atomIndex);

					// Match atom symbol
					currentAtomType = currentAtom.getSymbol();
					
						
						int nonsymmetricatom = 0;
						if (NEWDESC_PROPERTY.IsSymmetric.get(currentAtom) != null) nonsymmetricatom = NEWDESC_PROPERTY.IsSymmetric.get(currentAtom).intValue(); 
						if(nonsymmetricatom != 1) {
						
							outfile.print((moleculeIndex) + "," + currentAtom.getSymbol() + "."+ currentAtom.getID());				
							//Atom2endofMol descriptors
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_BondsToEnd.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_RotableBondCount.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_AtomCount.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_TPSA.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_TPSAperAtom.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_Volume.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_HAcount.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_HDcount.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Mol_PIsystemSize.get(currentAtom)));
							//Atom2endofBranch descriptors
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_BondsToEnd.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_RotableBondCount.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_AtomCount.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_TPSA.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_TPSAperAtom.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_Volume.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_HAcount.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_HDcount.get(currentAtom)));
							outfile.print("," + twoDecimalFormat.format(NEWDESC_PROPERTY.Branch_PIsystemSize.get(currentAtom)));
							outfile.print("\n");
						}
				} 

			}
			outfile.close();
		}
		

	}
	
	
}



