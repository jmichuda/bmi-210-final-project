from owlready2 import *
from src.generate_data import add_civic
import defopt

def main(input_path: str, civic_path:str, output_path: str):
	onto = get_ontology(input_path).load()
	with onto:
		onto = add_civic(onto, civic_path)
		onto.save(file = output_path)

if __name__ == "__main__":
	defopt.run(main)

