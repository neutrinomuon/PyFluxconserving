import hashlib
import requests

def calculate_sha256(url):
    response = requests.get(url)
    content = response.content
    sha256 = hashlib.sha256(content).hexdigest()
    return sha256

if __name__ == "__main__":
    import sys
    url = sys.argv[1]
    sha256 = calculate_sha256(url)
    print(sha256)
