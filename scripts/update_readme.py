with open('version.txt', 'r') as file:
    version = file.read().strip()

with open('README.md', 'r') as file:
    readme_content = file.read()

updated_readme = readme_content.replace('{VERSION_NUMBER}', version)

with open('README.md', 'w') as file:
    file.write(updated_readme)
