import re
import os
import fnmatch

# Regex patterns to identify different parts of the comment
func_regex = re.compile(r'/\*\s*@func: (.*?)\*/', re.DOTALL)
param_regex = re.compile(r'//\s*@param (.*?): (.*?)$', re.MULTILINE)
return_regex = re.compile(r'//\s*@return (.*?): (.*?)$', re.MULTILINE)
comment_block_regex = re.compile(
    r'/\*\s*@func: (.*?)\*/(.*?)//\s*@return (.*?): (.*?)$',
    re.DOTALL | re.MULTILINE
)

def load_gitignore_patterns(repo_path):
    """
    Loads and returns a list of patterns from the .gitignore file.
    """
    ignore_patterns = ['.git', '.gitignore', '.github', 'README.md', 'Cargo.lock', 'Cargo.toml', 'pyproject.toml', 'requirements.txt']
    gitignore_path = os.path.join(repo_path, '.gitignore')

    try:
        with open(gitignore_path, 'r') as file:
            file_patterns = file.readlines()
        file_patterns = [pattern.strip() for pattern in file_patterns if pattern.strip() != '' and not pattern.startswith('#')]
        ignore_patterns.extend(file_patterns)
    except FileNotFoundError:
        print(f"No .gitignore file found in {repo_path}. Continuing with default ignore patterns.")

    return ignore_patterns

def should_ignore_path(path, ignore_patterns, repo_path):
    """
    Determines if the given path matches any of the .gitignore patterns.
    """
    normalized_path = os.path.relpath(path, start=repo_path)
    for pattern in ignore_patterns:
        if fnmatch.fnmatch(normalized_path, pattern) or fnmatch.fnmatch(os.path.basename(path), pattern)  or normalized_path in pattern:
            return True
    return False

def detect_language(file_path):
    """
    Detects the programming language of a file based on its file extension.

    """
    _, file_extension = os.path.splitext(file_path)
    if file_extension == '.py':
        return "python"
    elif file_extension == '.rs':
        return "rust"
    else:
        return None

def extract_function_signatures_and_comments(file_content, language):
    """
    Extracts function signatures and their associated following comment blocks for Python and Rust.
    """
    if language == "python":
        pattern = re.compile(
            r'(def \w+\(.*?\))\s*:\s*(?:(?:#.*\n)|(?:\"\"\".*?\"\"\"|\'\'\'.*?\'\'\'))*\s*(/\*\s*@func:.*?\*/)',
            re.DOTALL
        )
    elif language == "rust":
        pattern = re.compile(
            r'(fn \w+\(.*?\))\s*->\s*.*?{\s*.*?(/\*\s*@func:.*?\*/)',
            re.DOTALL
        )
    else:
        raise ValueError(f"Unsupported language: {language}")

    matches = pattern.findall(file_content)
    formatted_matches = [(match[0], match[-1]) for match in matches]
    return formatted_matches

def format_comment_block(function_signature, comment_block):
    """
    Formats a comment block for README inclusion, fixing issues with parameter details and styling.
    """

    func_desc_match = re.search(r'/\*\s*@func: (.*?)\*/', comment_block, re.DOTALL)
    func_description = func_desc_match.group(1).strip() if func_desc_match else "No description available."


    if "->" in function_signature:
        func_name_params, return_type = function_signature.rsplit("->", 1)
        return_type = return_type.replace("{", "").strip()
        return_type = return_type[:-1].strip() if return_type.endswith(":") else return_type
    else:
        func_name_params = function_signature
        return_type = "None"

    param_start = func_name_params.find("(") + 1
    param_end = func_name_params.rfind(")")
    params_raw = func_name_params[param_start:param_end]
    params = params_raw.split(",")

    formatted_signature = f"### {func_name_params} -> {return_type} :\n"

    formatted_details = f"#### Description of function: {func_description}\n\n#### Parameters:\n"
    if params_raw.strip():
        for param in params:
            param = param.strip()
            formatted_details += f"- **{param}**\n"
    else:
        formatted_details += "None\n"

    if return_type:
        formatted_details += f"\n#### Returns:\n- **{return_type}**\n"

    return formatted_signature + formatted_details


def find_and_format_comments_in_file(file_path):
    """
    Finds and formats comments within a given file.
    """
    language = detect_language(file_path)
    if language is None:
        print(f"Unsupported file type: {file_path}")
        return ""

    with open(file_path, 'r') as file:
        print(f"processing file {file_path}")
        content = file.read()
        formatted_comments = ''
        function_blocks = extract_function_signatures_and_comments(content,language)

        for signature, comment in function_blocks:
            formatted_comments += format_comment_block(signature, comment)

        return formatted_comments

def scan_repository_for_comments(repo_path, ignore_patterns):
    """
    Scans the repository for comments to format and compile.
    """
    comments = ''
    for subdir, dirs, files in os.walk(repo_path):
        dirs[:] = [d for d in dirs if not should_ignore_path(os.path.join(subdir, d), ignore_patterns, repo_path)]
        for filename in files:
            file_path = os.path.join(subdir, filename)
            if not should_ignore_path(file_path, ignore_patterns, repo_path):
                comments += find_and_format_comments_in_file(file_path)
    return comments

def append_comments_to_readme(comments, readme_path):
    """
    Appends comments to the README.md file.
    """
    if not comments.strip():
        print("No new comments to append.")
        return

    try:
        with open(readme_path, 'r+', encoding='utf-8') as readme:
            readme_content = readme.read()
            if comments.strip() not in readme_content:
                readme.write('\n' + comments.strip())
    except FileNotFoundError:
        print(f"README.md file not found at path: {readme_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    repo_path = './'
    ignore_patterns = load_gitignore_patterns(repo_path)
    comments = scan_repository_for_comments(repo_path, ignore_patterns)
    readme_path = os.path.join(repo_path, 'README.md')
    append_comments_to_readme(comments, readme_path)
