import json
import yaml

def json_to_yaml(json_data):
    # Parse the JSON data
    packages = json.loads(json_data)
    
    # Prepare the YAML structure
    yaml_data = {
        'name': 'my_environment',
        'channels': ['defaults'],
        'dependencies': []
    }
    
    # Populate dependencies
    for package in packages:
        if 'name' in package and 'version' in package and 'build' in package:
            dependency = f"{package['name']}={package['version']}={package['build']}"
            yaml_data['dependencies'].append(dependency)
    
    # Convert to YAML
    return yaml.dump(yaml_data, default_flow_style=False)

def read_json_file(file_path):
    with open(file_path, 'r') as file:
        json_data = json.load(file)
    return json_data

# Example usage
json_input = read_json_file("env.json")
print(json_input)

yaml_output = json_to_yaml(json_input)
print(yaml_output)