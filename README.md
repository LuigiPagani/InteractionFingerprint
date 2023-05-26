# Interaction Fingerprint Computation

This code computes the interaction fingerprint between a ligand and a protein. The interaction fingerprint is a binary vector that indicates the presence or absence of specific types of interactions between the ligand and the protein. The code uses the RDKit library to handle molecular representations and calculations.

## Code Overview

The code performs the following steps:

1. Read the protein and ligand from files. The protein file is in PDB format, and the ligand file is in MOL2 format.
2. Convert the protein and ligand into RDKit molecules.
3. Break down the protein into residues.
4. Identify the residues that are close enough to the ligand based on a specified distance threshold.
5. Define interaction types using the `InteractionType` struct, which contains the name, SMARTS patterns for the ligand and protein, maximum distance, and angle constraints for each interaction type.
6. Check for each residue-ligand pair if any interaction type is satisfied.
7. Compute the interaction fingerprint by setting the corresponding bits for each satisfied interaction type.
8. Print the computed interaction fingerprint to the console.

## Dependencies

The code requires the following dependencies:

- RDKit: A collection of cheminformatics and machine learning tools for handling molecular data.

## Usage

1. Ensure that the RDKit library is properly installed and configured.
2. Provide the protein and ligand files in the PDB and MOL2 formats, respectively. Update the file paths in the code if necessary.
3. Compile and execute the code.
4. The interaction fingerprint will be printed to the console as a binary vector indicating the presence or absence of each interaction type.

Note: The code assumes that the ligand and protein files are properly formatted and compatible with the RDKit library.

## References

The parameters for the interaction types used in this code are based on the ProLIF paper by Cédric Bouysset and Sébastien Fiorucci. For more details, refer to the following paper:

- [ProLIF: An extended encoding for molecular interactions with protein-ligand fingerprints](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00548-6) (Cédric Bouysset and Sébastien Fiorucci, Journal of Cheminformatics, 2021)

Please cite the above paper if you use or reference the interaction fingerprint parameters in your work.
