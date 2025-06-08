import sys

def edit_cmake_file(file_name):
    """
    Comment the cuda enable command at `cuda.cmake`.
    """
    with open(file_name) as f:
        code_lines = f.readlines()

    to_comment_line = 'enable_language(CUDA)\n'
    comment_line = '# enable_language(CUDA)\n'

    if code_lines.count(to_comment_line) > 0:
        idx = code_lines.index(to_comment_line)
        code_lines[idx] = comment_line

        with open(file_name, 'w') as f:
            f.writelines(code_lines)


if __name__ == "__main__":
    file_name = sys.argv[1]
    edit_cmake_file(file_name)
