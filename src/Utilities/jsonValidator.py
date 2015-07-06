#!/usr/bin/python3.3

import sys
import json
import jsonschema

'''
This script validates the json against its schema
'''


def main(argc, argv):
    if not argc == 3:
        print("The validator must receive two parameters:")
        print("	* first parameter: the schema")
        print("	* second parameter: the json to validate")
        return 0

    with open(argv[1], "r") as f:
        schema = json.load(f)
    with open(argv[2], "r") as f:
        j_son = json.load(f)
    try:
        jsonschema.validate(j_son, schema)
        print("the json is valid")
    except jsonschema.ValidationError as e:
        print(e.message)
    except jsonschema.SchemaError as e:
        print(e)

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
