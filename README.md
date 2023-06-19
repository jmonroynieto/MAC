# MAC2.0
<sup>Copyright © 2023 bioinfomaticsCSU</sup>

Merging assemblies by using adjacency algebraic model and classification.

## About the rewrite
MAC2.x uses C++ to refactor the main contents of the original MAC, fixing errors in dealing with large files.
MAC2.0 removes the step of evaluation, therefore raw reads are not needed anymore.

## Dependency:
Please make sure that [MUMmer](https://github.com/mummer4/mummer) has been added into your environment's PATH variable.<sub><sup> MUMmer3 output is also supported. As a separate installation, you are responsible to uphold their licensing and you should include [their citation](https://github.com/mummer4/mummer/tree/master#description) too.</sub></sup>

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

### options
For CLI options see `MAC2.0 -h`.

## Notification:

- Because of some implementation issues, the MAs may be a little more than the original approach.

- If you have multiple fasta files (over 3) to merge, please try to merge any two of them first, then merge the result with the other files iteratively.

