/* 
 * Copyright (C) 2012 Patrik Rydberg <patrik.rydberg@gmail.com>
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

import java.io.InputStream;
import java.util.Arrays;
import java.util.List;
import java.util.TreeMap;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.atomtype.IAtomTypeMatcher;
import org.openscience.cdk.atomtype.SybylAtomTypeMatcher;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.graph.invariant.EquivalentClassPartitioner;
import org.openscience.cdk.graph.matrix.AdjacencyMatrix;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.manipulator.BondManipulator;

import NewAtomDescriptors.MoleculeKU.NEWDESC_PROPERTY;

public class MoleculeKU extends AtomContainer implements IMolecule {

	public enum NEWDESC_PROPERTY {
		SymmetryNumber,
		IsSymmetric,
		NrofSymmetricSites,
		AtomCount2BranchEnd,
		RelSpan {
			@Override
			public String getLabel() {
				return "RS";
			}
		}, Mol_RotableBondCount, Mol_BondsToEnd, Mol_AtomCount, Mol_TPSA, Mol_TPSAperAtom, 
		Mol_Volume, Mol_HAcount, Mol_HDcount, Mol_PIsystemSize, Branch_RotableBondCount, 
		Branch_BondsToEnd, Branch_AtomCount, Branch_TPSA, Branch_TPSAperAtom, Branch_Volume, 
		Branch_HAcount, Branch_HDcount, Branch_PIsystemSize;

		public String  getLabel()  { return "";};

		public void set(IAtom atom, Number value) {
			atom.setProperty(toString(), value);
		}
		
		public Number get(IAtom atom) {
			Object o = atom.getProperty(toString());
			return (o==null)?null:(Number)o;
		}

		public String atomProperty2String(IAtom atom) {
			return String.format("%s:%s",getLabel(),get(atom));
		}

	}

	// Local variables
	private static final long serialVersionUID = 1L;	
	private int debug = 0;


	// Constructor
	// This constructor also calls the methods that calculate all needed properties
	// This constructor is the only way to create MoleculeKU and Atom objects, -there is no add() method
	public MoleculeKU(IAtomContainer iAtomContainer) throws CloneNotSupportedException
	{
		// Calls the constructor in org.openscience.cdk.AtomContainer
		// Atoms are stored in the array atoms[] and accessed by getAtom() and setAtom()
		super(iAtomContainer);			
		int number = 1;
		for (int atomIndex=0; atomIndex < iAtomContainer.getAtomCount(); atomIndex++) {
			iAtomContainer.getAtom(atomIndex).setID(String.valueOf(number));
			number++;
		}
	}

	public void calculateRelativeSpan() throws CloneNotSupportedException{


		int[][] adjacencyMatrix = AdjacencyMatrix.getMatrix(this);

		// Calculate the maximum topology distance
		// Takes an adjacency matrix and outputs and MinTopDist matrix of the same size
		int[][] minTopDistMatrix = PathTools.computeFloydAPSP(adjacencyMatrix);

		// Find the longest Path of all, "longestMaxTopDistInMolecule"
		double longestMaxTopDistInMolecule = 0;
		double currentMaxTopDist = 0;
		for(int x = 0; x < this.getAtomCount(); x++){
			for(int y = 0; y < this.getAtomCount(); y++){
				currentMaxTopDist =  minTopDistMatrix[x][y];
				if(currentMaxTopDist > longestMaxTopDistInMolecule) longestMaxTopDistInMolecule = currentMaxTopDist;
			}
		}

		// Find the "longest shortestPath" for each atom

		// ITERATE REFERENCE ATOMS
		for (int refAtomNr=0; refAtomNr < this.getAtomCount(); refAtomNr++){

			// ITERATE COMPARISON ATOMS
			double highestMaxTopDistInMatrixRow = 0;
			IAtom refAtom;
			for (int compAtomNr = 0; compAtomNr < this.getAtomCount(); compAtomNr++){
				if(highestMaxTopDistInMatrixRow < minTopDistMatrix[refAtomNr][compAtomNr]) highestMaxTopDistInMatrixRow = minTopDistMatrix[refAtomNr][compAtomNr];
			}	

			refAtom = this.getAtom(refAtomNr);
			// Set the relative span of the Atom
			NEWDESC_PROPERTY.RelSpan.set(refAtom,(highestMaxTopDistInMatrixRow / longestMaxTopDistInMolecule));
		}
	}
	
	public int[] findAtomsatEndOfMolecule() throws CloneNotSupportedException{
		//returns integer vector, if a number is 1 then that number has relative span 1 and is at end of molecule, else it is 0

		int[] endofmolatomindices = new int [this.getAtomCount()];
		IAtom atom;
		double relspan;
		
		// ITERATE ATOMS
		for (int AtomNr=0; AtomNr < this.getAtomCount(); AtomNr++){
			atom = this.getAtom(AtomNr);
			relspan = NEWDESC_PROPERTY.RelSpan.get(atom).doubleValue();
			if (relspan == 1.0) {
				endofmolatomindices[AtomNr] = 1;
				//System.out.println(AtomNr);
			}
			else endofmolatomindices[AtomNr] = 0;
		}
		
		return endofmolatomindices;
	}
	
	public int[] findAtomsatEndOfBranch() throws CloneNotSupportedException{
		//returns integer vector, if a number is 1 then that number has relative span 1 and is at end of molecule, else it is 0

		int[] endofbranchatomindices = new int [this.getAtomCount()];
		IAtom atom;
		double relspan;
		double neighborrelspan;
		boolean endofbranchatom;
		
		// ITERATE ATOMS
		for (int AtomNr=0; AtomNr < this.getAtomCount(); AtomNr++){
			endofbranchatom = false;
			atom = this.getAtom(AtomNr);
			
			relspan = NEWDESC_PROPERTY.RelSpan.get(atom).doubleValue();
			if (relspan == 1.0) {
				//this atom is at the end of the molecule, no need to look further
				endofbranchatom = true;
				//System.out.println(AtomNr);
			}
			else {
				//this atom is not at the end of the molecule, let's check if it has connected atoms that are closer to the end of the molecule
				//any atom with another atom within three bonds that have a higher relspan is either not at a branchend, or it's a single atom branch 
				//and single atom branches should be excluded
				int[][] adjacencyMatrix = AdjacencyMatrix.getMatrix(this);
				int[][] minTopDistMatrix = PathTools.computeFloydAPSP(adjacencyMatrix);
				endofbranchatom = true; //now lets set it to false if neighbor closer to end of mol exists
				for (int i = 0; i < minTopDistMatrix.length; i++){
					if (minTopDistMatrix[AtomNr][i] < 4){
						//System.out.println(AtomNr + ":" + i);
						neighborrelspan = NEWDESC_PROPERTY.RelSpan.get(this.getAtom(i)).doubleValue();
						if (neighborrelspan > relspan) {
							//the neighbor has a higher relspan, this is not a branchend
							endofbranchatom = false;
						}
					}
				}
			}	
			if (endofbranchatom) {
				endofbranchatomindices[AtomNr] = 1;
				//System.out.println(AtomNr);
			}
			else endofbranchatomindices[AtomNr] = 0;
			
		}
		
		return endofbranchatomindices;
	}
	
	public IAtomContainer[] getAtoms2EndOfMolMolecules(int[] EndofMolAtomIndices, String prefix) throws CloneNotSupportedException{
		//returns array of molecules for each atom
		int[] endofmolatomindices = EndofMolAtomIndices;
		IAtomContainer[] AtomMolecules = new IAtomContainer[this.getAtomCount()];
		IAtomContainer Molecule2EndofMol;
		IAtomContainer Alt2EndofMol;
		//get the number of bonds between all atoms to minTopDistMatrix
		int[][] adjacencyMatrix = AdjacencyMatrix.getMatrix(this);
		int[][] minTopDistMatrix = PathTools.computeFloydAPSP(adjacencyMatrix);
		int shortest = 0;
		int distance;
		int ClosestAtomIndex;
		int[] AlternativeAtomIndices;
		int AlternativeAtomCount;
		boolean OnOtherBranch;
		//list of atomindices to keep in molecule
		int[] Atoms2Keep = new int[this.getAtomCount()];
		int[] Atoms2KeepAlt = new int[this.getAtomCount()];
			
		// ITERATE ATOMS and build their Molecule2EndofMol
		for (int AtomNr = 0; AtomNr < this.getAtomCount(); AtomNr++){
			IAtom ThisAtom = this.getAtom(AtomNr);
			AlternativeAtomIndices = new int[9];
			//reset the Atoms2Keep array
			for (int i = 0; i < this.getAtomCount(); i++){
				Atoms2Keep[i] = 0;
			}
			//First let's check that this is not and endofmol atom
			if (endofmolatomindices[AtomNr] == 1){
				//Molecule2EndofMole can be set to null
				AtomMolecules[AtomNr] = null;
				shortest = 0;
				//System.out.println("endofmolatom: " + AtomNr);
			}
			else{
				//This is not an atom at the end of the molecule.
				//let's find the closest endofmol atom
				shortest = 99;
				ClosestAtomIndex = 999;
				AlternativeAtomCount = 0;
				for (int trialendofmolatomindex = 0; trialendofmolatomindex < endofmolatomindices.length; trialendofmolatomindex++){
					if(endofmolatomindices[trialendofmolatomindex] == 1){
						//this is and endofmolatom index
						distance = minTopDistMatrix[AtomNr][trialendofmolatomindex];
						if (distance < shortest) {
							shortest = distance;
							ClosestAtomIndex = trialendofmolatomindex;
						}
						if (distance == shortest){
							//first check symmetry, if symmetric ignore this fragment
							if (NEWDESC_PROPERTY.SymmetryNumber.get(this.getAtom(trialendofmolatomindex)) != NEWDESC_PROPERTY.SymmetryNumber.get(this.getAtom(ClosestAtomIndex))){
								//we got a non-symmetric fragment at same distance
								AlternativeAtomIndices[AlternativeAtomCount] = trialendofmolatomindex;
								AlternativeAtomCount++;
							}
						}
					}
				}
				//now the closest endofmol atom index is ClosestAtomIndex
				//Let's build an IAtomContainer of the atoms in between AtomNr and ClosestAtomIndex
				for (int trialatomindex = 0; trialatomindex < this.getAtomCount(); trialatomindex++){
					OnOtherBranch = false;
					//first, atoms in this IAtomContainer should have a distance to ClosestAtomIndex that is shorter than AtomNr
					//and it should not be AtomNr atom
					if (trialatomindex != AtomNr && minTopDistMatrix[trialatomindex][ClosestAtomIndex] <= minTopDistMatrix[AtomNr][ClosestAtomIndex]){
						//and the trial atom should be closer to AtomNr than the distance between AtomNr and ClosesAtomIndex
						if((minTopDistMatrix[trialatomindex][ClosestAtomIndex] < minTopDistMatrix[AtomNr][ClosestAtomIndex] && (minTopDistMatrix[trialatomindex][AtomNr] <= minTopDistMatrix[AtomNr][ClosestAtomIndex] 
						|| trialatomindex == ClosestAtomIndex))
						|| (minTopDistMatrix[trialatomindex][ClosestAtomIndex] == minTopDistMatrix[AtomNr][ClosestAtomIndex] && NEWDESC_PROPERTY.RelSpan.get(this.getAtom(trialatomindex)).doubleValue() > NEWDESC_PROPERTY.RelSpan.get(this.getAtom(AtomNr)).doubleValue())
						){
							//one final check. The atom must not be closer to another ClosestAtomIndex, 
							//because then it is on another branch than the one we are interested in traversing
							for (int otherendofmolatomindex = 0; otherendofmolatomindex < endofmolatomindices.length; otherendofmolatomindex++){
								if(endofmolatomindices[otherendofmolatomindex] == 1){
									if (otherendofmolatomindex != ClosestAtomIndex 
									&& minTopDistMatrix[trialatomindex][ClosestAtomIndex] > minTopDistMatrix[trialatomindex][otherendofmolatomindex] 
                                    && minTopDistMatrix[trialatomindex][AtomNr] != minTopDistMatrix[ClosestAtomIndex][AtomNr]){
										OnOtherBranch = true;
									}
								}
							}
							if (OnOtherBranch == false){
								//Molecule2EndofMol.addAtom(this.getAtom(trialatomindex));
								Atoms2Keep[trialatomindex] = 1;
								//System.out.println(AtomNr + ":" + trialatomindex);
							}
						}
					}
				}
					
			}
			
			//now let's build the new molecule
			//let's keep the atom of interest in the molecule to get correct bond counts
			Molecule2EndofMol = new AtomContainer(this);
			for (int i = 0; i < this.getAtomCount(); i++){
				if(Atoms2Keep[i] == 0 && i != AtomNr){
					Molecule2EndofMol.removeAtomAndConnectedElectronContainers(this.getAtom(i));
				}
			}
			
			//now let's check if there are any alternative atoms and if so generate their fragments for comparison
			for (int altindex = 0; altindex < 4; altindex++){
				int altatomindex = AlternativeAtomIndices[altindex];
				//reset the Atoms2KeepAlt array
				for (int i = 0; i < this.getAtomCount(); i++){
					Atoms2KeepAlt[i] = 0;
				}
				if(altatomindex > 0){
					//we got an alternative atom which is not symmetric to the first one
					//let's build it's fragment
					for (int trialatomindex = 0; trialatomindex < this.getAtomCount(); trialatomindex++){
						OnOtherBranch = false;
						//first, atoms in this IAtomContainer should have a distance to ClosestAtomIndex that is shorter than AtomNr
						//and it should not be AtomNr atom
						if (trialatomindex != AtomNr && minTopDistMatrix[trialatomindex][altatomindex] <= minTopDistMatrix[AtomNr][altatomindex]){
							//and the trial atom should be closer to AtomNr than the distance between AtomNr and ClosesAtomIndex
							if((minTopDistMatrix[trialatomindex][altatomindex] < minTopDistMatrix[AtomNr][altatomindex] && (minTopDistMatrix[trialatomindex][AtomNr] <= minTopDistMatrix[AtomNr][altatomindex] 
							|| trialatomindex == altatomindex))
							|| (minTopDistMatrix[trialatomindex][altatomindex] == minTopDistMatrix[AtomNr][altatomindex] && NEWDESC_PROPERTY.RelSpan.get(this.getAtom(trialatomindex)).doubleValue() > NEWDESC_PROPERTY.RelSpan.get(this.getAtom(AtomNr)).doubleValue())
							){
								//one final check. The atom must not be closer to another ClosestAtomIndex, 
								//because then it is on another branch than the one we are interested in traversing
								for (int otherendofmolatomindex = 0; otherendofmolatomindex < endofmolatomindices.length; otherendofmolatomindex++){
									if(endofmolatomindices[otherendofmolatomindex] == 1){
										if (otherendofmolatomindex != altatomindex 
										&& minTopDistMatrix[trialatomindex][altatomindex] > minTopDistMatrix[trialatomindex][otherendofmolatomindex] 
	                                    && minTopDistMatrix[trialatomindex][AtomNr] != minTopDistMatrix[altatomindex][AtomNr]){
											OnOtherBranch = true;
										}
									}
								}
								if (OnOtherBranch == false){
									//Molecule2EndofMol.addAtom(this.getAtom(altatomindex));
									Atoms2KeepAlt[trialatomindex] = 1;
									//System.out.println(AtomNr + ":" + altatomindex);
								}
							}
						}
					}
					//now let's build the alternative new molecule
					//let's keep the atom of interest in the molecule to get correct bond counts
					Alt2EndofMol = new AtomContainer(this);
					for (int i = 0; i < this.getAtomCount(); i++){
						if(Atoms2KeepAlt[i] == 0 && i != AtomNr){
							Alt2EndofMol.removeAtomAndConnectedElectronContainers(this.getAtom(i));
						}
					}
					
					//compare the alternative fragment to Molecule2EndofMol
					if(Alt2EndofMol.getAtomCount() < Molecule2EndofMol.getAtomCount()) {
						Molecule2EndofMol = Alt2EndofMol;
					}
					else if (Alt2EndofMol.getAtomCount() == Molecule2EndofMol.getAtomCount()){
						//check rotable bonds count
						int rotbondsalt = FragmentRotableBondsCount(Alt2EndofMol); 
						int rotbonds = FragmentRotableBondsCount(Molecule2EndofMol);
						if(rotbondsalt > rotbonds) {
							Molecule2EndofMol = Alt2EndofMol;
						}
						else if (rotbondsalt == rotbonds){
							//the two fragments are of same atom count and have the same number of rotable bonds
							//now let's check for MW,if still the same let's skip the alternative fragment, the difference will not be significant
							WeightDescriptor mwd = new WeightDescriptor();
							IDescriptorResult mwalt = mwd.calculate(Alt2EndofMol).getValue();
							IDescriptorResult mw = mwd.calculate(Molecule2EndofMol).getValue();
							DoubleResult ifmwalt = (DoubleResult) mwalt; 
							DoubleResult ifmw = (DoubleResult) mw;
							if(ifmwalt.doubleValue() < ifmw.doubleValue()) {
								Molecule2EndofMol = Alt2EndofMol;
							}
						}
					}					
				}				
			}
			
			//before this, keep the atom of interest
			//now compute bond descriptors
			int rotbonds = FragmentRotableBondsCount(Molecule2EndofMol);
			Molecule2EndofMol.setProperty(prefix + "_RotableBondCount", rotbonds);
			Molecule2EndofMol.setProperty(prefix + "_BondsToEnd", shortest);
			
			//remove the atom of interest
			Molecule2EndofMol.removeAtomAndConnectedElectronContainers(this.getAtom(AtomNr));
			//now compute other properties
			//Atom Count
			Molecule2EndofMol.setProperty(prefix + "_AtomCount", Molecule2EndofMol.getAtomCount());
			//Topological polar surface area
			TPSADescriptor tpsa = new TPSADescriptor();
			DescriptorValue tpsavalue = tpsa.calculate(Molecule2EndofMol);
			DoubleResult tpsavaluedouble = (DoubleResult) tpsavalue.getValue();
			Number tpsavaluenumber = (Number) tpsavaluedouble.doubleValue();
			Molecule2EndofMol.setProperty(prefix + "_TPSA", tpsavaluedouble);
			//Topological polar surface area per Atom, a measure of the polar fraction
			double TPSAperAtom = 0;
			if(Molecule2EndofMol.getAtomCount() != 0) TPSAperAtom = tpsavaluedouble.doubleValue()/Molecule2EndofMol.getAtomCount();
			Molecule2EndofMol.setProperty(prefix + "_TPSAperAtom", TPSAperAtom);
			//Volume descriptor
			VABCDescriptor vabc = new VABCDescriptor();
			DescriptorValue vabcvalue = vabc.calculate(Molecule2EndofMol);
			DoubleResult vabcvaluedouble = (DoubleResult) vabcvalue.getValue();
			Number vabcvaluenumber = (Number) vabcvaluedouble.doubleValue();
			Molecule2EndofMol.setProperty(prefix + "_Volume", vabcvaluedouble);
			//H-bond acceptor count descriptor
			HBondAcceptorCountDescriptor hacd = new HBondAcceptorCountDescriptor();
			DescriptorValue hacdvalue = hacd.calculate(Molecule2EndofMol);
			IntegerResult hacdvalueint = (IntegerResult) hacdvalue.getValue();
			Number hacdvaluenumber = (Number) hacdvalueint.intValue();
			Molecule2EndofMol.setProperty(prefix + "_HAcount", hacdvalueint);
			//H-bond donor count descriptor
			HBondDonorCountDescriptor hdcd = new HBondDonorCountDescriptor();
			DescriptorValue hdcdvalue = hdcd.calculate(Molecule2EndofMol);
			IntegerResult hdcdvalueint = (IntegerResult) hdcdvalue.getValue();
			Number hdcdvaluenumber = (Number) hdcdvalueint.intValue();
			Molecule2EndofMol.setProperty(prefix + "_HDcount", hdcdvalueint);
			//Largest PI system descriptor
			LargestPiSystemDescriptor lpsd = new LargestPiSystemDescriptor();
			DescriptorValue lpsdvalue = lpsd.calculate(Molecule2EndofMol);
			IntegerResult lpsdvalueint = (IntegerResult) lpsdvalue.getValue();
			Number lpsdvaluenumber = (Number) lpsdvalueint.intValue();
			Molecule2EndofMol.setProperty(prefix + "_PISystemSize", lpsdvalueint);
			//Mannhold LogP descriptor
			MannholdLogPDescriptor logp = new MannholdLogPDescriptor();
			DescriptorValue logpvalue = logp.calculate(Molecule2EndofMol);
			DoubleResult logpvaluedouble = (DoubleResult) logpvalue.getValue();
			Molecule2EndofMol.setProperty(prefix + "_logP", logpvaluedouble);
			
			
			if (prefix == "Mol"){
				NEWDESC_PROPERTY.Mol_RotableBondCount.set(ThisAtom,rotbonds);
				NEWDESC_PROPERTY.Mol_BondsToEnd.set(ThisAtom, shortest);
				NEWDESC_PROPERTY.Mol_AtomCount.set(ThisAtom, Molecule2EndofMol.getAtomCount());
				NEWDESC_PROPERTY.Mol_TPSA.set(ThisAtom, tpsavaluenumber);
				NEWDESC_PROPERTY.Mol_TPSAperAtom.set(ThisAtom, TPSAperAtom);
				NEWDESC_PROPERTY.Mol_Volume.set(ThisAtom, vabcvaluenumber);
				NEWDESC_PROPERTY.Mol_HAcount.set(ThisAtom, hacdvaluenumber);
				NEWDESC_PROPERTY.Mol_HDcount.set(ThisAtom, hdcdvaluenumber);
				NEWDESC_PROPERTY.Mol_PIsystemSize.set(ThisAtom, lpsdvaluenumber);
			}
			if (prefix == "Branch"){
				NEWDESC_PROPERTY.Branch_RotableBondCount.set(ThisAtom,rotbonds);
				NEWDESC_PROPERTY.Branch_BondsToEnd.set(ThisAtom, shortest);
				NEWDESC_PROPERTY.Branch_AtomCount.set(ThisAtom, Molecule2EndofMol.getAtomCount());
				NEWDESC_PROPERTY.Branch_TPSA.set(ThisAtom, tpsavaluenumber);
				NEWDESC_PROPERTY.Branch_TPSAperAtom.set(ThisAtom, TPSAperAtom);
				NEWDESC_PROPERTY.Branch_Volume.set(ThisAtom, vabcvaluenumber);
				NEWDESC_PROPERTY.Branch_HAcount.set(ThisAtom, hacdvaluenumber);
				NEWDESC_PROPERTY.Branch_HDcount.set(ThisAtom, hdcdvaluenumber);
				NEWDESC_PROPERTY.Branch_PIsystemSize.set(ThisAtom, lpsdvaluenumber);
			}
			
			//assign the "best" fragment to the return array, with all properties set
			AtomMolecules[AtomNr] = Molecule2EndofMol;
		}
		
		return AtomMolecules;
	}
	

	private int FragmentRotableBondsCount(IAtomContainer ac) {
		int rotableBondsCount = 0;
		for (IBond bond : ac.bonds()) {
			IAtom atom0 = bond.getAtom(0);
			IAtom atom1 = bond.getAtom(1);
            if (bond.getOrder() == CDKConstants.BONDORDER_SINGLE) {
				if ((BondManipulator.isLowerOrder(ac.getMaximumBondOrder(atom0), IBond.Order.TRIPLE)) && 
					(BondManipulator.isLowerOrder(ac.getMaximumBondOrder(atom1), IBond.Order.TRIPLE))) {
					if (!bond.getFlag(CDKConstants.ISINRING)) {
						//now we got rotable bonds, including terminal bonds
						//do not count bond if either atom has only one other heavy atom attached
						if((atom0.getFormalNeighbourCount() - atom0.getImplicitHydrogenCount()) > 1 && (atom1.getFormalNeighbourCount() - atom1.getImplicitHydrogenCount()) > 1 ){
							//now exclude amides, thioamides, and sulfonamides
							String atom0type = atom0.getAtomTypeName();
							String atom1type = atom1.getAtomTypeName();
							if ((atom0type.equals("N.amide") && atom1type.equals("C.sp2")) || (atom0type.equals("C.sp2") && atom1type.equals("N.amide"))){
								//an amide/peptide bond is not a rotable bond
							}
							else if((atom0type.equals("N.thioamide") && atom1type.equals("C.sp2")) || (atom0type.equals("C.sp2") && atom1type.equals("N.thioamide"))){
								//a thioamide bond is not a rotable bond
							}
							else if((atom0type.equals("N.sp3") && atom1type.equals("S.onyl")) || (atom0type.equals("S.onyl") && atom1type.equals("N.sp3"))){
								//a sulfonamide bond is not a rotable bond
							} 
							else {
								rotableBondsCount += 1;
							}
						}
                        
					}
				}
			}
		}
		return rotableBondsCount;
	}

	// Symmetric atoms have identical values in the array from getTopoEquivClassbyHuXu
	public void setSymmetryNumbers() throws CDKException{
		Atom atom;
		//set charges so that they are not null
		for(int atomIndex = 0; atomIndex < this.getAtomCount(); atomIndex++){
			atom = (Atom) this.getAtom(atomIndex);
			atom.setCharge((double) atom.getFormalCharge());
		}
		//compute symmetry
		EquivalentClassPartitioner symmtest = new EquivalentClassPartitioner((AtomContainer) this);
		int[] symmetryNumbersArray = symmtest.getTopoEquivClassbyHuXu((AtomContainer) this);
		symmetryNumbersArray[0]=0;//so we can count the number of symmetric sites for each atom without double counting for the ones with the highest symmetrynumber
		int symmsites;
		for(int atomIndex = 0; atomIndex < this.getAtomCount(); atomIndex++){
			symmsites = 0;
			atom = (Atom) this.getAtom(atomIndex);
			NEWDESC_PROPERTY.SymmetryNumber.set(atom,symmetryNumbersArray[atomIndex+1]);
			// Compute how many symmetric sites the atom has, 1=only itself
			symmsites = FindInArray(symmetryNumbersArray,symmetryNumbersArray[atomIndex+1]);
			NEWDESC_PROPERTY.NrofSymmetricSites.set(atom,symmsites);
		}
	}
	
	public static int FindInArray(int[] arr, int numToFind) {
		int occurence=0;
		for (int i = 0; i < arr.length; i++) { 
			if (arr[i] == numToFind) occurence++;
		}
		return occurence;
	}

	
	public void setID(String id){
		super.setID(id);
	}

	public String toString(){
		for(int atomNr=0; atomNr < this.getAtomCount(); atomNr++) System.out.println(this.getAtom(atomNr).toString());
		return "MoleculeKU " + super.toString();
	}
	
}





