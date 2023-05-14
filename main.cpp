// This code is used to compute the interaction fingerprint of a ligand and a protein.

// The first step is to read the protein and ligand from files.
// The protein file is in PDB format, and the ligand file is in MOL2 format.

// Once the protein and ligand have been read, they are converted to RDKit molecules.

// Next, the protein is broken down into residues.

// The residues that are close enough to the ligand are then identified.
// The distance between two atoms is calculated using the Euclidean distance formula.

// An interaction fingerprint is a binary vector that indicates which types of interactions are present between the ligand and the protein.
// The interaction types are defined in the `InteractionType` struct.

// The interaction fingerprint is computed by looping over the residues that are close to the ligand.
// For each residue, the code checks if the residue has the specified interaction type with the ligand.
// If the residue has the interaction type, the corresponding bit in the interaction fingerprint is set to 1.

// The interaction fingerprint is then printed to the console.


#include <iostream>
#include <vector>
#include <cmath>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Geometry/point.h>
#include <Numerics/Vector.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/ROMol.h>
#include <boost/shared_ptr.hpp>

using namespace RDKit;
using namespace RDGeom;
using namespace RDKit;
using namespace RDGeom;
// Define the interaction types and their SMARTS patterns
struct InteractionType {
    std::string name;
    std::string ligandSmarts;
    std::string proteinSmarts;
    double maxDistance;
    double minAngle;
    double maxAngle;
};
std::vector<InteractionType> interactionTypes = {
        {"Anionic", "[-1]", "[+1]", 4.5, 0.0, 0.0},
        {"Cationic", "[+1]", "[-1]", 4.5, 0.0, 0.0},
        {"CationPi", "[+1]", "a1aaaaa1", 4.5, 0.0, 30.0},
        {"PiCation", "a1aaaaa1", "[+1]", 4.5, 0.0, 30.0},
        {"PiStacking", "a1aaaaa1", "a1aaaaa1", 6.0, 0.0, 90.0},
        {"EdgeToFace", "a1aaaaa1", "a1aaaaa1", 6.0, 50.0, 90.0},
        {"FaceToFace", "a1aaaaa1", "a1aaaaa1", 4.5, 0.0, 40.0},
        {"HBAcceptor", "[N,O,F;-1;!+1]", "[#7,#8,#16][H]", 3.5, 130.0, 180.0},
        {"HBDonor", "[#7,#8,#16][H]", "[N,O,F;-1;!+1]", 3.5, 130.0, 180.0},
        {"XBAcceptor", "[#7,#8,P,S,Se,Te;a;!+1][*]", "[#6,#7,Si,F,Cl,Br,I]-[Cl,Br,I,At]", 3.5, 130.0, 180.0},
        {"XBDonor", "[#6,#7,Si,F,Cl,Br,I]-[Cl,Br,I,At]", "[#7,#8,P,S,Se,Te;a;!+1][*]", 3.5, 130.0, 180.0},
        {"MetalAcceptor", "[O,N;-1;!+1]", "[Ca,Cd,Co,Cu,Fe,Mg,Mn,Ni,Zn]", 2.8, 0.0, 0.0},
        {"MetalDonor", "[Ca,Cd,Co,Cu,Fe,Mg,Mn,Ni,Zn]", "[O,N;-1;!+1]", 2.8, 0.0, 0.0},
        {"Hydrophobic", "[#6,#16,F,Cl,Br,I,At;+0]", "[#6,#16,F,Cl,Br,I,At;+0]", 4.5, 0.0, 0.0}
};
bool checkInteraction(const InteractionType& interaction, const ROMol& ligand, const ROMol& residue) {
    // Find substructures in the ligand and residue
    std::vector<MatchVectType> ligandMatches, residueMatches;
    ROMol* ligandQuery = SmartsToMol(interaction.ligandSmarts);
    ROMol* residueQuery = SmartsToMol(interaction.proteinSmarts);
    SubstructMatch(ligand, *ligandQuery, ligandMatches);
    SubstructMatch(residue, *residueQuery, residueMatches);
    // Check geometric constraints for each combination of substructures
    for (const auto& ligandMatch : ligandMatches) {
        for (const auto& residueMatch : residueMatches) {
            Point3D ligandPos = ligand.getConformer().getAtomPos(ligandMatch[0].second);
            Point3D residuePos = residue.getConformer().getAtomPos(residueMatch[0].second);
            double distance = (ligandPos - residuePos).length();
            if (distance <= interaction.maxDistance) {
                if (interaction.minAngle == 0.0 && interaction.maxAngle == 0.0) {
                    return true;
                } else {
                    Point3D ligandDir = ligand.getConformer().getAtomPos(ligandMatch[0].first) - ligandPos;
                    Point3D residueDir = residue.getConformer().getAtomPos(residueMatch[0].first) - residuePos;
                    double angle = ligandDir.angleTo(residueDir) * 180.0 / M_PI;
                    if (angle >= interaction.minAngle && angle <= interaction.maxAngle) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}
ExplicitBitVect computeInteractionFingerprint(const ROMol& protein, const ROMol& ligand) {
    // Break down the protein into residues
    std::vector<boost::shared_ptr<ROMol>> residues = MolOps::getMolFrags(protein); // Change the type here
    // Find residues close enough to the ligand
    std::vector<boost::shared_ptr<ROMol>> closeResidues;
    for (const auto& residue: residues) {
        bool isClose = false;
        for (const auto& proteinAtom : residue->atoms()) {
            for (const auto& ligandAtom : ligand.atoms()) {
                Point3D proteinPos = protein.getConformer().getAtomPos(proteinAtom->getIdx());
                Point3D ligandPos = ligand.getConformer().getAtomPos(ligandAtom->getIdx());
                double distance = (proteinPos - ligandPos).length();
                if (distance <= 7.0) {
                    isClose = true;
                    break;
                }
            }
            if (isClose) {
                closeResidues.push_back(residue);
                break;
            }
        }
    }
    // Initialize the interaction fingerprint
    ExplicitBitVect fp(interactionTypes.size() * closeResidues.size(), false);
    // Compute the interaction fingerprint
    for (size_t i = 0; i < closeResidues.size(); ++i) {
        for (size_t j = 0; j < interactionTypes.size(); ++j) {
            if (checkInteraction(interactionTypes[j], ligand, *closeResidues[i])) {
                if (i < closeResidues.size() && j < interactionTypes.size()) {
                    if (i < closeResidues.size() && j < interactionTypes.size()) {
                        fp.setBit(i * interactionTypes.size() + j);
                    } else {
                        std::cerr << "Error: index out of bounds" << std::endl;
                    }
                }
            }
        }
    }
    return fp;
}

int main() {
    // Read the protein and ligand from files
    std::string proteinFile = "protein.pdb";
    std::string ligandFile = "crystal.mol2";
    boost::shared_ptr<ROMol> protein(PDBFileToMol(proteinFile));
    boost::shared_ptr<ROMol> ligand(Mol2FileToMol(ligandFile));

    if (!protein) {
        std::cerr << "Error reading protein file: " << proteinFile << std::endl;
        return 1;
    }
    if (!ligand) {
        std::cerr << "Error reading ligand file: " << ligandFile << std::endl;
        return 1;
    }

    // Add conformers to the protein and ligand molecules
    protein->addConformer(new Conformer(protein->getNumAtoms()), true);
    ligand->addConformer(new Conformer(ligand->getNumAtoms()), true);

    // Compute the interaction fingerprint
    ExplicitBitVect fp = computeInteractionFingerprint(*protein, *ligand);

    // Print the fingerprint
    std::cout << "Interaction fingerprint: ";
    for (unsigned int i = 0; i < fp.getNumBits(); ++i) {
        std::cout << (fp[i] ? "1" : "0");
    }
    std::cout << std::endl;

    return 0;
}