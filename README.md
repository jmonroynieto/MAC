# MAC2.0

Merging assemblies by using adjacency algebraic model and classification.

## About the rewrite

MAC2.0 uses C++ to refactor the main contents of the original MAC, fixing errors in dealing large files.
MAC2.0 removes the step of evaluation, which does't need raw reads as input any more.

## Dependency:


// TODO check MUMmer4 is compatible
Please make sure that [MUMmer](https://github.com/mummer4/mummer) has been added into environment's PATH.

## Usage: 
### compile from source:
`g++ MAC2.0.cpp -o MAC2.0`

### Simple use case
`MAC2.0 <query.fa> <reference.fa>`

This requires the following directory structure.
```
analysisDir
├── input
│   ├── query.fa
│   ├── reference.fa
```
- The input files need to be placed in the ./input folder.
- The output file will be output to the ./output folder.

## Notification:

- Because of some implementation issues, the MAs may be a little more than the original approach.

- If you have multiple contig files (over 3) to merge, please try to merge any two of them first, then merge the result with the other files iteratively.
