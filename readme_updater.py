import os
import re
import fnmatch


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
        # Remove the section in place
        del content[start_of_functions:end_of_functions]

    # Rewrite the file with the modified content
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
            # Determine the start and end of the comment block based on the file type
            comment_start = '/*' if function_match_rust else '"""'
            comment_end = '*/' if function_match_rust else '"""'

            # Check for comment block after the function
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

            # Extract @func:, @param, and @return information
            if "@func:" in comment_block:
                func_desc_pattern = re.compile(r'@func: (.*?)(?=(@param|@return|$))', re.DOTALL)
                param_desc_pattern = re.compile(r'@param (.*?): (.*?)(?=(@param|@return|$))', re.DOTALL)
                return_desc_pattern = re.compile(r'@return\s+(.*?):\s+(.*?)(?=(@param|@func|$))', re.DOTALL)

                # Extracting descriptions
                func_desc_match = re.search(func_desc_pattern, comment_block)
                param_desc_matches = re.finditer(param_desc_pattern, comment_block)
                return_desc_match = re.search(return_desc_pattern, comment_block)

                # Formatting descriptions
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
# def parse_function_signature(signature, language):
#     """
#     Parses the function signature to extract parameter names, types, and optionality,
#     ensuring the return type is not included with parameters.
#     """
#     params = []
#     if language == "rust":
#         # Rust: Adjust to stop before "-> ReturnType"
#         param_section = re.search(r'fn\s+\w+\((.*?)\)\s*->', signature)
#         if param_section:
#             param_pattern = re.compile(r'(\w+)\s*:\s*([^,]+)')
#             params = param_pattern.findall(param_section.group(1))
#             params = [(name, ty.strip(), 'optional' if 'Option<' in ty else '') for name, ty in params]
#     elif language == "python":
#         # Python: Adjust to correctly parse default values and type hints.
#         # Split signature at the first colon (assumes no colons in default values for simplicity).
#         param_section = signature.split("):", 1)[0] + ")"
#         param_pattern = re.compile(r'(\w+)\s*:\s*([^=,]+)(\s*=\s*.+)?')
#         matches = param_pattern.findall(param_section)
#         params = [(name, ty.strip(), 'optional' if default else '') for name, ty, default in matches]
#
#     return params

def parse_function_signature(signature, language):
    """
    Parses the function signature to extract parameter names, types, optionality,
    and return type.
    """
    params = []
    return_type = ""
    if language == "rust":
        # Extracting parameters and return type for Rust
        param_section = re.search(r'fn\s+\w+\((.*?)\)\s*->\s*(.*)', signature)
        if param_section:
            param_pattern = re.compile(r'(\w+)\s*:\s*([^,]+)')
            params = param_pattern.findall(param_section.group(1))
            params = [(name, ty.strip(), 'optional' if 'Option<' in ty else '') for name, ty in params]
            return_type = param_section.group(2).strip()  # Capture the return type
    elif language == "python":
        # Adjust for Python if using type hints for return types
        # This is a simplistic approach; adjust as necessary for your codebase
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

                # Adjust return section to match parameters format
                if description["return"]:
                    # Assuming the return description is in the format "variable: description"
                    return_type = return_type.replace("{","").replace(" ","")
                    return_var, return_desc = description["return"].split(":", 1)
                    readme.write('  - **Returns:**\n')
                    readme.write(f'    - `{return_var}` (*{return_type}*): {return_desc.strip()}\n')
                readme.write('\n')

        # Similar adjustment for Python functions
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

                # Adjust return section to match parameters format for Python functions
                if description["return"]:
                    return_var, return_desc = description["return"].split(":", 1)
                    readme.write('  - **Returns:**\n')
                    readme.write(f'    - `{return_var}` (*{return_type}*): {return_desc.strip()}\n')
                readme.write('\n')


if __name__ == "__main__":
    # Example usage
    repo_path = './'
    readme_path = os.path.join(repo_path, 'README.md')
    repo_base_url = 'https://github.com/zkdavis/PyParamo'
    script_filename = os.path.basename(__file__)
    extra_ignores = ['.git', '.gitignore', '.github', 'README.md',
                       'Cargo.lock', 'Cargo.toml', 'pyproject.toml',
                       'requirements.txt', 'distribs_unfinished', script_filename]
    remove_existing_functions_section(readme_path)
    pyf, rsf = document_functions(repo_path, extra_ignores)
    write_functions_to_readme(rsf,pyf, readme_path,repo_base_url)

