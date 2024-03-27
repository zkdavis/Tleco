import os
import re
import fnmatch
import subprocess
import toml




def remove_existing_functions_section(readme_path):
    with open(readme_path, 'r+', encoding='utf-8') as readme:
        content = readme.readlines()

    start_of_functions = None
    end_of_functions = None
    for i, line in enumerate(content):
        if line.startswith('# Functions'):
            start_of_functions = i
        if start_of_functions is not None and "# " in line and "Functions" not in line and i > start_of_functions:
            end_of_functions = i
            break
    if(end_of_functions is None):
        end_of_functions = len(content)

    if start_of_functions is not None and end_of_functions is not None:
        del content[start_of_functions:end_of_functions]

    with open(readme_path, 'w', encoding='utf-8') as readme:
        readme.writelines(content)


def parse_gitignore(gitignore_path):
    with open(gitignore_path, 'r') as file:
        ignore_patterns = file.read().splitlines()
    return ignore_patterns


def should_ignore(path, ignore_patterns):
    path = path.split("/")
    paths = []
    paths.extend(path)
    ret = False
    for path in paths:
        ret = any(fnmatch.fnmatch(path, pattern) for pattern in ignore_patterns) or any(path in pattern for pattern in ignore_patterns)
        if(ret == True):
            return True
    return ret


def extract_functions_from_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    functions = []
    for i, line in enumerate(lines):
        function_match_rust = re.match(r'^\s*pub fn\s+[^\(]+\([^)]*\)\s*->\s*[^{]+{', line)
        function_match_python = re.match(r'^\s*def\s+[^\(]+\([^)]*\)\s*:', line)
        if function_match_rust or function_match_python:
            function_signature = line.strip()
            comment_block = ""
            comment_start = '/*' if function_match_rust else '"""'
            comment_end = '*/' if function_match_rust else '"""'

            if i + 1 < len(lines) and lines[i + 1].strip().startswith(comment_start):
                in_comment_block = False
                for j in range(i + 1, len(lines)):
                    comment_line = lines[j].strip()
                    if comment_line.startswith(comment_start):
                        in_comment_block = True
                        comment_line = comment_line[len(comment_start):].strip()  # Remove the start token
                    if in_comment_block:
                        comment_block += comment_line + "\n"
                    if comment_line.endswith(comment_end):
                        if in_comment_block:  # Found the end of the comment block
                            comment_block = comment_block[:-(len(comment_end) + 1)].strip()  # Remove the end token
                            break

            if "@func:" in comment_block:
                func_desc_pattern = re.compile(r'@func: (.*?)(?=(@param|@return|$))', re.DOTALL)
                param_desc_pattern = re.compile(r'@param (.*?): (.*?)(?=(@param|@return|$))', re.DOTALL)
                return_desc_pattern = re.compile(r'@return\s+(.*?):\s+(.*?)(?=(@param|@func|$))', re.DOTALL)

                func_desc_match = re.search(func_desc_pattern, comment_block)
                param_desc_matches = re.finditer(param_desc_pattern, comment_block)
                return_desc_match = re.search(return_desc_pattern, comment_block)

                func_description = func_desc_match.group(1).strip() if func_desc_match else ""
                param_descriptions = {m.group(1).strip(): m.group(2).strip() for m in param_desc_matches}
                return_description = f"{return_desc_match.group(1).strip()}: {return_desc_match.group(2).strip()}" if return_desc_match else ""

                func_description = func_description.replace("\n","")
                for keys in param_descriptions.keys():
                    param_descriptions[keys] = param_descriptions[keys].replace("\n","")
                return_description = return_description.replace("\n","")

                description = {"description": func_description, "params": param_descriptions, "return": return_description}

                if function_match_rust:
                    functions.append((i, function_signature, description))
                elif function_match_python:
                    functions.append((i, function_signature, description))

    return functions


def document_functions(repo_path, ignore_patterns=[]):
    gitignore_path = os.path.join(repo_path, '.gitignore')
    if os.path.exists(gitignore_path):
        ignore_patterns.extend(parse_gitignore(gitignore_path))

    rust_functions = []
    python_functions = []

    for root, dirs, files in os.walk(repo_path):
        for file in files:
            file_path = os.path.join(root, file)
            relative_path = os.path.relpath(file_path, repo_path)

            if should_ignore(relative_path, ignore_patterns):
                continue

            if file.endswith('.rs') or file.endswith('.py'):
                functions = extract_functions_from_file(file_path)
                for line_number, function_signature, comment_block in functions:
                    if file.endswith('.rs'):
                        rust_functions.append((relative_path, line_number, function_signature, comment_block))
                    else:
                        python_functions.append((relative_path, line_number, function_signature, comment_block))

    python_functions.sort(key=lambda x: x[2])  # Sort Python functions alphabetically by signature
    return python_functions, rust_functions


def parse_function_signature(signature, language):
    """
    Parses the function signature to extract parameter names, types, optionality,
    and return type.
    """
    params = []
    return_type = ""
    if language == "rust":
        param_section = re.search(r'fn\s+\w+\((.*?)\)\s*->\s*(.*)', signature)
        if param_section:
            param_pattern = re.compile(r'(\w+)\s*:\s*([^,]+)')
            params = param_pattern.findall(param_section.group(1))
            params = [(name, ty.strip(), 'optional' if 'Option<' in ty else '') for name, ty in params]
            return_type = param_section.group(2).strip()  # Capture the return type
    elif language == "python":
        param_section = signature.split("):", 1)
        if len(param_section) > 1 and "->" in param_section[1]:
            params_str, return_type_str = param_section
            param_pattern = re.compile(r'(\w+)\s*:\s*([^=,]+)(\s*=\s*.+)?')
            params = param_pattern.findall(params_str + ")")
            params = [(name, ty.strip(), 'optional' if default else '') for name, ty, default in params]
            return_type = return_type_str.split("->", 1)[1].strip()

    return params, return_type

def write_functions_to_readme(rust_functions, python_functions, readme_path, repo_base_url):
    with open(readme_path, 'a') as readme:
        readme.write('# Functions\n')
        readme.write('This section is not complete and will be updated over time\n')

        if rust_functions:
            readme.write('### Rust Functions\n')
            for file_path, line_number, function_signature, description in rust_functions:
                readme.write('\n')
                function_name = re.search(r'fn\s+(\w+)', function_signature).group(1)
                location_url = f'{repo_base_url}/blob/master/{file_path}#L{line_number}'
                location_text = f'[{file_path}:L{line_number}]({location_url})'
                readme.write(f'- **{function_name}** - {description["description"]} [see in {location_text}]\n')
                readme.write('  - **Parameters:**\n')
                params, return_type = parse_function_signature(function_signature, "rust")
                for param, param_type, optional in params:
                    optional_str = " (optional)" if optional else ""
                    readme.write(f'    - `{param}` (*{param_type}*{optional_str}): {description["params"].get(param, "No description")}\n')
                if description["return"]:
                    return_type = return_type.replace("{","").replace(" ","")
                    return_var, return_desc = description["return"].split(":", 1)
                    readme.write('  - **Returns:**\n')
                    readme.write(f'    - `{return_var}` (*{return_type}*): {return_desc.strip()}\n')
                readme.write('\n')
        if python_functions:
            readme.write('### Python Functions\n')
            for file_path, line_number, function_signature, description in python_functions:
                readme.write('\n')
                function_name = re.search(r'def\s+(\w+)', function_signature).group(1)
                location_url = f'{repo_base_url}/blob/master/{file_path}#L{line_number}'
                location_text = f'[{file_path}:L{line_number}]({location_url})'
                readme.write(f'- **{function_name}** - {description["description"]} [see in {location_text}]\n')
                readme.write('  - **Parameters:**\n')
                params, return_type = parse_function_signature(function_signature, "python")
                for param, param_type, optional in params:
                    optional_str = "(optional)" if optional else ""
                    readme.write(f'    - `{param}` (*{param_type}*{optional_str}): {description["params"].get(param, "No description")}\n')
                if description["return"]:
                    return_var, return_desc = description["return"].split(":", 1)
                    readme.write('  - **Returns:**\n')
                    readme.write(f'    - `{return_var}` (*{return_type}*): {return_desc.strip()}\n')
                readme.write('\n')


def generate_requirements_txt(project_path, exclude_dirs=None,  additional_deps=None):
    """
    Uses pipreqs to generate requirements.txt for the given project path.
    """
    command = ['pipreqs', '--force', project_path]

    if exclude_dirs:
        exclude_dirs_str = ','.join(exclude_dirs)
        command += ['--ignore', exclude_dirs_str]

    subprocess.run(command, check=True)

    if additional_deps:
        requirements_path = f"{project_path}/requirements.txt"
        with open(requirements_path, "a") as req_file:
            for dep in additional_deps:
                req_file.write(f"{dep}\n")

def update_pyproject_toml(project_path):
    """
    Update pyproject.toml using poetry.
    """
    with open(f"{project_path}/requirements.txt", "r") as req_file:
        for line in req_file:
            package_name, version = line.strip().split('==')
            subprocess.run(['poetry', 'add', f'{package_name}=={version}'], check=True)


def get_python_version_requirement(pyproject_path):
    """
    Reads the Python version requirement from pyproject.toml.
    """
    with open(pyproject_path, 'r', encoding='utf-8') as file:
        data = toml.load(file)
        python_version = data['tool']['poetry']['dependencies']['python']
    return python_version


def update_readme_requirements(project_path):
    """
    Updates README.md with a Requirements section and dynamically adds the Python version requirement.
    """
    pyproject_path = f"{project_path}/pyproject.toml"
    python_version_requirement = get_python_version_requirement(pyproject_path)

    readme_path = f"{project_path}/README.md"
    requirements_text = "# Requirements\n\nThe following Python packages are used in this project:\n\n"

    with open(f"{project_path}/requirements.txt", "r", encoding='utf-8') as req_file:
        for line in req_file:
            package_name, version = line.strip().split('==')
            requirements_text += f"- {package_name} {version}\n"
    requirements_text += f"#### Python Version\n\nThis project requires Python {python_version_requirement}.\n"
    with open(readme_path, "r", encoding='utf-8') as readme_file:
        readme_content = readme_file.read()
    pattern = re.compile(r'(# Requirements\n.*?)(#### Python Version\n.*?)(?=# |\Z)', flags=re.DOTALL)
    new_readme_content = re.sub(pattern, '', readme_content, 1)  # Replace only the first match
    new_readme_content += "\n" + requirements_text
    with open(readme_path, "w", encoding='utf-8') as readme_file:
        readme_file.write(new_readme_content)


def update_project_dependencies_and_docs(project_path,ignore_dir=None,additional_deps=None):
    """
    Main function to update the project dependencies and documentation.
    """
    print("Generating requirements.txt...")
    generate_requirements_txt(project_path,ignore_dir,additional_deps=additional_deps)

    print("Updating pyproject.toml...")
    update_pyproject_toml(project_path)

    print("Updating README.md...")
    update_readme_requirements(project_path)

    print("Update complete.")

def read_version_from_file(version_file_path='VERSION.txt'):
    """
    Reads the project version from the specified VERSION.txt file.
    """
    with open(version_file_path, 'r', encoding='utf-8') as file:
        version = file.readline().strip()
        return version

def update_pyproject_version_from_file(version_file_path='VERSION.txt', pyproject_path='pyproject.toml'):
    """
    Update the project version in pyproject.toml based on the version found in VERSION.txt.
    """
    version = read_version_from_file(version_file_path)

    with open(pyproject_path,'r', encoding='utf-8') as f:
        data = toml.load(f)
    data['tool']['poetry']['version'] = version
    data['project']['version'] = version
    with open(pyproject_path, 'w', encoding='utf-8') as f:
        toml.dump(data, f)

    print(f"Updated version in {pyproject_path} to {version}")

def update_cargo_version_from_file(version_file_path='VERSION.txt', cargo_path='Cargo.toml'):
    """
    Update the project version in Cargo.toml based on the version found in VERSION.txt.
    """
    version = read_version_from_file(version_file_path)

    with open(cargo_path,'r', encoding='utf-8') as f:
        data = toml.load(f)
    data['package']['version'] = version
    with open(cargo_path, 'w', encoding='utf-8') as f:
        toml.dump(data, f)

    print(f"Updated version in {cargo_path} to {version}")


if __name__ == "__main__":
    repo_path = './'
    readme_path = os.path.join(repo_path, 'README.md')
    repo_base_url = 'https://github.com/zkdavis/Tleco'
    script_filename = os.path.basename(__file__)
    extra_ignores = ['.git', '.gitignore', '.github', 'README.md',
                       'Cargo.lock', 'Cargo.toml', 'pyproject.toml',
                       'requirements.txt', 'distribs_unfinished', script_filename]
    remove_existing_functions_section(readme_path)
    pyf, rsf = document_functions(repo_path, extra_ignores)
    write_functions_to_readme(rsf,pyf, readme_path,repo_base_url)
    gitignore_path = os.path.join(repo_path, '.gitignore')
    if os.path.exists(gitignore_path):
        extra_ignores.extend(parse_gitignore(gitignore_path))
    add_deps=['maturin==1.5']
    update_project_dependencies_and_docs('./',ignore_dir=extra_ignores,additional_deps=add_deps)
    update_pyproject_version_from_file()
    update_cargo_version_from_file()
