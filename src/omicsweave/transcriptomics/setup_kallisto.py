"""Download kallisto Windows binary and mouse transcriptome index."""
import os
import sys
import ssl
import zipfile
import urllib.request
import time

from omicsweave.core.project import project_dir

KALLISTO_URL = "https://github.com/pachterlab/kallisto/releases/download/v0.51.1/kallisto_windows-v0.51.1.zip"


def download_file(url, dest, label=""):
    # 保持默认的证书与主机名校验；网络环境无法校验时请手动下载
    ctx = ssl.create_default_context()
    
    req = urllib.request.Request(url, headers={
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'
    })
    
    print(f"  Downloading {label}...", flush=True)
    print(f"  URL: {url}", flush=True)
    
    start = time.time()
    response = urllib.request.urlopen(req, timeout=120, context=ctx)
    total = int(response.headers.get('Content-Length', 0))
    
    downloaded = 0
    chunk_size = 256 * 1024
    last_report = 0
    
    with open(dest, 'wb') as f:
        while True:
            chunk = response.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            mb = downloaded / (1024*1024)
            if mb - last_report >= 1:
                if total:
                    pct = downloaded / total * 100
                    print(f"    {mb:.1f}/{total/(1024*1024):.1f} MB ({pct:.0f}%)", flush=True)
                else:
                    print(f"    {mb:.1f} MB", flush=True)
                last_report = mb
    
    elapsed = time.time() - start
    print(f"  Done: {downloaded/(1024*1024):.1f} MB in {elapsed:.0f}s", flush=True)


def main():
    # 外部程序安装到项目目录（运行时解析，以便 --project-dir 生效），而不是安装后的包目录
    tools_dir = os.path.join(str(project_dir()), "tools")
    kallisto_exe = os.path.join(tools_dir, "kallisto", "kallisto.exe")
    print("=" * 60, flush=True)
    print("  Setup kallisto for Windows", flush=True)
    print("=" * 60, flush=True)
    
    os.makedirs(tools_dir, exist_ok=True)
    
    # 1. Download kallisto
    if os.path.exists(kallisto_exe):
        print(f"kallisto already exists: {kallisto_exe}", flush=True)
    else:
        zip_path = os.path.join(tools_dir, "kallisto.zip")
        try:
            download_file(KALLISTO_URL, zip_path, "kallisto")
        except Exception as e:
            print(f"  ERROR: {e}", flush=True)
            # Try alternate version
            alt_url = "https://github.com/pachterlab/kallisto/releases/download/v0.50.1/kallisto_windows-v0.50.1.zip"
            print(f"  Trying alternate version...", flush=True)
            download_file(alt_url, zip_path, "kallisto v0.50.1")
        
        print("  Extracting...", flush=True)
        with zipfile.ZipFile(zip_path, 'r') as z:
            z.extractall(tools_dir)
        os.remove(zip_path)
        
        # Find the exe
        for root, dirs, files in os.walk(tools_dir):
            for f in files:
                if f == "kallisto.exe":
                    found = os.path.join(root, f)
                    print(f"  kallisto found at: {found}", flush=True)
                    break
    
    # 2. Verify
    import subprocess
    try:
        result = subprocess.run([kallisto_exe, "version"], capture_output=True, text=True)
        print(f"\nkallisto version: {result.stdout.strip()}", flush=True)
        print("Setup complete!", flush=True)
    except Exception as e:
        # Try to find it
        for root, dirs, files in os.walk(tools_dir):
            for f in files:
                if f == "kallisto.exe":
                    exe = os.path.join(root, f)
                    print(f"  Found exe at: {exe}", flush=True)
                    result = subprocess.run([exe, "version"], capture_output=True, text=True)
                    print(f"  Version: {result.stdout.strip()}", flush=True)


if __name__ == "__main__":
    main()
