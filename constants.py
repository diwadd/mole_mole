from enum import Enum

TRAIN_CVS_FILE_NAME = "../champs-scalar-coupling/train.csv"
STRUCTURE_XYZ_FILES_PATH = "../champs-scalar-coupling/structures/"


class ReadingMode(Enum):
    SHORT = "short"
    LONG = "all"