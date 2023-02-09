# change shell working directory to ADFR suite executables
cd ADFRsuite-1.0/bin

# set root directory/project directory
ROOT_DIR="/home/romfahrury/vina"

# set working directory/current directory
WORKING_DIR="$ROOT_DIR/ADFRsuite-1.0/bin"

# set mol2 files directory for copying
COPY_FROM_DIR="$ROOT_DIR/dataset/processed/mol2"

# set output directory
OUTPUT_DIR="$ROOT_DIR/dataset/processed/pdbqt/ligand"

# shell parameter variable 1 (input name)
INPUT_NAME=$1

# shell parameter variable 2 (output name)
OUTPUT_NAME=$2

# copy mol2 files into working directory
cp $COPY_FROM_DIR/$INPUT_NAME $WORKING_DIR/

# run the ADFR prepare ligand program
./prepare_ligand -l $INPUT_NAME -o $OUTPUT_DIR/$OUTPUT_NAME

# remove copied files for convenient
rm $INPUT_NAME