"""
Download mouse transcriptome from GENCODE/Ensembl with multiple mirrors.
Tries faster sources first (GENCODE > Ensembl mirrors).
"""
import os
import sys
import time
import urllib.request
import ssl

from omicsweave.core.project import project_dir

OUTPUT_NAME = "Mus_musculus.GRCm39.cdna.all.fa.gz"

# Multiple download sources, ordered by reliability
URLS = [
    # GENCODE (hosted on EBI, often faster)
    ("GENCODE M35", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.transcripts.fa.gz"),
    # GENCODE M34
    ("GENCODE M34", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.transcripts.fa.gz"),
    # Ensembl (main)
    ("Ensembl 113", "https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"),
    # Ensembl (US mirror)
    ("Ensembl US", "https://useast.ensembl.org/info/data/ftp/index.html"),
]


def download_with_progress(url, output_path, timeout=600):
    """Download with chunked progress display."""
    # 保持默认的证书与主机名校验；网络环境无法校验时请手动下载
    ctx = ssl.create_default_context()

    req = urllib.request.Request(url, headers={
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'
    })

    start_time = time.time()

    try:
        response = urllib.request.urlopen(req, timeout=30, context=ctx)
    except Exception as e:
        print(f"    Connection failed: {e}", flush=True)
        return False

    total_size = int(response.headers.get('Content-Length', 0))
    total_mb = total_size / (1024 * 1024) if total_size else 0
    print(f"    File size: {total_mb:.1f} MB", flush=True)

    downloaded = 0
    chunk_size = 1024 * 256  # 256 KB chunks
    last_report = 0

    with open(output_path, 'wb') as f:
        while True:
            elapsed = time.time() - start_time
            if elapsed > timeout:
                print(f"\n    Timeout after {timeout}s", flush=True)
                return False

            try:
                chunk = response.read(chunk_size)
            except Exception as e:
                print(f"\n    Read error: {e}", flush=True)
                return False

            if not chunk:
                break

            f.write(chunk)
            downloaded += len(chunk)
            mb_done = downloaded / (1024 * 1024)

            # Report every 2 MB
            if mb_done - last_report >= 2:
                speed = mb_done / elapsed if elapsed > 0 else 0
                if total_size:
                    pct = downloaded / total_size * 100
                    eta = (total_size - downloaded) / (downloaded / elapsed) if elapsed > 0 and downloaded > 0 else 0
                    print(f"    {mb_done:.1f}/{total_mb:.1f} MB ({pct:.0f}%) - {speed:.1f} MB/s - ETA: {eta:.0f}s", flush=True)
                else:
                    print(f"    {mb_done:.1f} MB - {speed:.1f} MB/s", flush=True)
                last_report = mb_done

    elapsed = time.time() - start_time
    final_mb = downloaded / (1024 * 1024)
    speed = final_mb / elapsed if elapsed > 0 else 0
    print(f"    Complete: {final_mb:.1f} MB in {elapsed:.0f}s ({speed:.1f} MB/s)", flush=True)
    return True


def main():
    # 参考数据下载到项目目录（运行时解析，以便 --project-dir 生效），而不是安装后的包目录
    output_dir = os.path.join(str(project_dir()), "reference")
    output_file = os.path.join(output_dir, OUTPUT_NAME)
    print("=" * 60, flush=True)
    print("  Download Mouse Transcriptome (GRCm39)", flush=True)
    print("=" * 60, flush=True)

    os.makedirs(output_dir, exist_ok=True)

    # Check if already exists
    if os.path.exists(output_file):
        size_mb = os.path.getsize(output_file) / (1024 * 1024)
        if size_mb > 30:
            print(f"File already exists: {output_file} ({size_mb:.1f} MB)", flush=True)
            return output_file

    for name, url in URLS:
        print(f"\nTrying {name}:", flush=True)
        print(f"  URL: {url}", flush=True)

        # For GENCODE, save with a different name then rename
        temp_file = output_file + ".downloading"
        if os.path.exists(temp_file):
            os.remove(temp_file)

        success = download_with_progress(url, temp_file, timeout=600)

        if success and os.path.exists(temp_file):
            size_mb = os.path.getsize(temp_file) / (1024 * 1024)
            if size_mb > 30:
                # Rename to final name
                if os.path.exists(output_file):
                    os.remove(output_file)
                os.rename(temp_file, output_file)
                print(f"\nSaved to: {output_file} ({size_mb:.1f} MB)", flush=True)
                return output_file
            else:
                print(f"    File too small ({size_mb:.1f} MB), trying next source...", flush=True)
                os.remove(temp_file)
        else:
            if os.path.exists(temp_file):
                os.remove(temp_file)
            print(f"    Failed, trying next source...", flush=True)

    print("\nAll sources failed. Please download manually:", flush=True)
    print("  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/", flush=True)
    print(f"  Save to: {output_file}", flush=True)
    return None


if __name__ == "__main__":
    result = main()
    if result:
        print(f"\nTranscriptome ready: {result}", flush=True)
