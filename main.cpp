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
// The parameters for the interactions come from the ProLIF paper https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00548-6 by Cédric Bouysset and Sébastien Fiorucci.


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
#include <omp.h>

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

Point3D computeRingCentroid(const ROMol& mol, const MatchVectType& match) {
    Point3D centroid(0, 0, 0);
    int atomCount = 0;
    for (const auto& pair : match) {
        Point3D pos = mol.getConformer().getAtomPos(pair.second);
        centroid += pos;
        atomCount++;
    }
    centroid /= static_cast<double>(atomCount);
    return centroid;
}

bool checkInteraction(const InteractionType& interaction, const ROMol& ligand, const ROMol& residue) {
    std::vector<MatchVectType> ligandMatches, residueMatches;
    ROMol* ligandQuery = SmartsToMol(interaction.ligandSmarts);
    ROMol* residueQuery = SmartsToMol(interaction.proteinSmarts);
    SubstructMatch(ligand, *ligandQuery, ligandMatches);
    SubstructMatch(residue, *residueQuery, residueMatches);
    for (const auto& ligandMatch : ligandMatches) {
        for (const auto& residueMatch : residueMatches) {
            Point3D ligandPos, residuePos;
            if (interaction.name == "PiStacking" || interaction.name == "CationPi") {
                ligandPos = computeRingCentroid(ligand, ligandMatch);
                residuePos = computeRingCentroid(residue, residueMatch);
            } else {
                ligandPos = ligand.getConformer().getAtomPos(ligandMatch[0].second);
                residuePos = residue.getConformer().getAtomPos(residueMatch[0].second);
            }
                }
            }
    return false;
}

ExplicitBitVect computeInteractionFingerprint(const ROMol& protein, const ROMol& ligand) {
    // Break down the protein into residues
    std::vector<boost::shared_ptr<ROMol>> residues = MolOps::getMolFrags(protein); // Change the type here

    // Initialize the interaction fingerprint
    ExplicitBitVect fp(interactionTypes.size() * residues.size(), false);

    // Find residues close enough to the ligand
    std::vector<boost::shared_ptr<ROMol>> closeResidues;
#pragma omp parallel for
    for (int i = 0; i < residues.size(); ++i) {
        const auto& residue = residues[i];
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
#pragma omp critical
                closeResidues.push_back(residue);
                break;
            }
        }
    }

    // Compute the interaction fingerprint
#pragma omp parallel for
    for (int i = 0; i < closeResidues.size(); ++i) {
        for (size_t j = 0; j < interactionTypes.size(); ++j) {
            if (checkInteraction(interactionTypes[j], ligand, *closeResidues[i])) {
#pragma omp critical
                {
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


std::vector<int> calculateInteractionCounts(const ExplicitBitVect& fp) {
    std::vector<int> interactionCounts(interactionTypes.size(), 0);
    for (unsigned int i = 0; i < fp.getNumBits(); ++i) {
        interactionCounts[i % interactionTypes.size()] += fp[i] ? 1 : 0;
    }
    return interactionCounts;
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

    std::vector<int> interactionCounts = calculateInteractionCounts(fp);

    // Print the interaction counts
    std::cout << "Interaction counts: ";
    for (int count : interactionCounts) {
        std::cout << count << " ";
    }
    std::cout << std::endl;

    return 0;

    return 0;
}
