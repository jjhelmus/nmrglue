""" nmrglue test setup module """
# location of the data directory
import os.path
TESTS_DIR = os.path.dirname(__file__)
NMRGLUE_ROOT = os.path.dirname(TESTS_DIR)
DATA_DIR = os.path.join(NMRGLUE_ROOT, 'data')
